import json
import re
from pathlib import Path

import citations
import provider
import retrieval
import verifier

REPO_ROOT = Path(__file__).resolve().parent.parent
C1_INDEX_DIR = REPO_ROOT / "data" / "processed" / "rag" / "c1_index"
PARAMETERS_PATH = REPO_ROOT / "docs" / "plan" / "PARAMETERS.md"
MANIFEST_PATH = REPO_ROOT / "data" / "processed" / "build_manifest.json"

MAX_TOOL_TURNS = 6

SYSTEM_PROMPT = """You are the methodology Q&A agent for an explainable cell-line selection \
tool. Researchers ask you why a parameter has the value it does, what its status is (cited, \
derived, EDA-inherited, or invented), what paper backs it, and what limitation is disclosed \
about it.

You have four tools: search_docs, get_parameter, get_citation, get_manifest_value. Use them --
never answer from memory. If a tool returns nothing (an empty result), say plainly that the \
project's documentation does not contain that information -- do not guess, estimate, or invent \
a value, a status, or a citation.

When you give your final answer (no more tool calls needed), write it as plain text, not JSON. \
Every number you state must come from a tool result, at the same precision the tool gave it. \
Cite any paper as [Author Year] exactly as a get_citation/search_docs result gives it."""

TOOL_SCHEMAS = [
    {
        "type": "function",
        "function": {
            "name": "search_docs",
            "description": "Search the project documentation (C1 corpus: docs/**.md, research summaries, module docstrings) for passages relevant to a query.",
            "parameters": {
                "type": "object",
                "properties": {"query": {"type": "string"}, "top_k": {"type": "integer", "default": 5}},
                "required": ["query"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "get_parameter",
            "description": "Look up a parameter by name (e.g. 'R_EXPONENT', 'min_lineage_n', 'LAYER_WEIGHTS') in docs/plan/PARAMETERS.md. Returns every matching row: value, where it's defined, its status, and notes.",
            "parameters": {
                "type": "object",
                "properties": {"name": {"type": "string"}},
                "required": ["name"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "get_citation",
            "description": "Resolve a citation like 'Lewis 2020' or a data source name against the registered docs/research/*/PAPERS.md and docs/reference/DATA_SOURCES.md entries.",
            "parameters": {
                "type": "object",
                "properties": {"key": {"type": "string"}},
                "required": ["key"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "get_manifest_value",
            "description": "Look up a key in data/processed/build_manifest.json (row counts, coverage stats, join-loss reports).",
            "parameters": {
                "type": "object",
                "properties": {"key": {"type": "string"}},
                "required": ["key"],
            },
        },
    },
]

_PARAMETER_ROW_RE = re.compile(r"^\|(.+)\|\s*$")


def get_parameter(name):
    if not PARAMETERS_PATH.exists():
        return []
    matches = []
    for line in PARAMETERS_PATH.read_text(encoding="utf-8").splitlines():
        row_match = _PARAMETER_ROW_RE.match(line)
        if not row_match:
            continue
        cells = [c.strip() for c in row_match.group(1).split("|")]
        if all(re.fullmatch(r":?-+:?", c) for c in cells if c):
            continue  # separator row
        row_text = " | ".join(cells)
        if name.lower() in row_text.lower():
            matches.append(row_text)
    return matches


def get_citation(key):
    registry = citations.build_registry()
    entry = citations.resolve(key, registry)
    return entry if entry else None


def get_manifest_value(key):
    if not MANIFEST_PATH.exists():
        return None
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    if key in manifest:
        return manifest[key]
    # dotted-path lookup, e.g. "coverage.measured" -- build_manifest.json nests per-table stats
    node = manifest
    for part in key.split("."):
        if isinstance(node, dict) and part in node:
            node = node[part]
        else:
            return None
    return node


_c1_index_cache = None


def search_docs(query, top_k=5):
    global _c1_index_cache
    if _c1_index_cache is None:
        if not (C1_INDEX_DIR / "chunks.jsonl").exists():
            return []
        _c1_index_cache = retrieval.HybridIndex.load(C1_INDEX_DIR)
    return _c1_index_cache.search(query, top_k=top_k)


TOOL_FUNCTIONS = {
    "search_docs": lambda args: search_docs(args["query"], args.get("top_k", 5)),
    "get_parameter": lambda args: get_parameter(args["name"]),
    "get_citation": lambda args: get_citation(args["key"]),
    "get_manifest_value": lambda args: get_manifest_value(args["key"]),
}


def _evidence_record_from_tool_results(tool_results_text):
    values = []
    for text in tool_results_text:
        for token in verifier.extract_numbers_raw(text):
            cleaned = token.rstrip("%").replace(",", "")
            try:
                value = float(cleaned)
            except ValueError:
                continue
            values.append(value / 100.0 if token.endswith("%") else value)
    return {"values": values}


def answer(question, registry=None, max_turns=MAX_TOOL_TURNS):
    if registry is None:
        registry = citations.build_registry()

    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": question},
    ]
    tool_results_text = []

    for _ in range(max_turns):
        message = provider.chat(messages, tools=TOOL_SCHEMAS)
        tool_calls = getattr(message, "tool_calls", None)
        if not tool_calls:
            final_text = message.content or ""
            context_text = "\n\n".join(tool_results_text)
            evidence_record = _evidence_record_from_tool_results(tool_results_text)
            check = verifier.verify(final_text, evidence_record, registry)
            abstained = _looks_like_abstention(final_text)
            return {
                "answer": final_text,
                "context": context_text,
                "abstained": abstained,
                "verification": check,
            }

        messages.append({
            "role": "assistant", "content": message.content or "",
            "tool_calls": [tc.model_dump() for tc in tool_calls],
        })
        for tool_call in tool_calls:
            name = tool_call.function.name
            args = json.loads(tool_call.function.arguments or "{}")
            try:
                result = TOOL_FUNCTIONS.get(name, lambda a: None)(args)
            except Exception as exc:
                # A tool failure (a bad manifest key, a malformed search) must not crash the
                # whole answer -- the model sees the error and can adapt or abstain honestly.
                result = {"error": f"{type(exc).__name__}: {exc}"}
            result_text = json.dumps(result, default=str)
            tool_results_text.append(result_text)
            messages.append({
                "role": "tool", "tool_call_id": tool_call.id, "name": name, "content": result_text,
            })

    return {
        "answer": "Could not reach a final answer within the tool-call turn limit.",
        "context": "\n\n".join(tool_results_text),
        "abstained": True,
        "verification": None,
    }


_ABSTENTION_KEYWORDS = (
    "no data", "not found", "no evidence", "does not exist", "doesn't exist", "no such",
    "cannot find", "can't find", "not a registered", "not available", "no record",
    "not in this project", "does not contain", "not documented", "i do not have", "i don't have",
)


def _looks_like_abstention(text):
    lowered = text.lower()
    return any(keyword in lowered for keyword in _ABSTENTION_KEYWORDS)
