import json

import citations as citations_module
import kg_tools
import provider
import verifier

MAX_TOOL_TURNS = 10
# Past this many turns, tool schemas are withheld from the request entirely -- a model call
# with no `tools` argument physically cannot return a tool_call, which forces a text answer
# instead of relying on the model's own judgment to stop exploring. Observed live: even with an
# explicit "make at most 3-4 tool calls" instruction, meta/llama-3.3-70b-instruct sometimes kept
# calling tools well past that budget and never produced a final answer within 10 turns -- this
# is the deterministic backstop for that failure mode, not a second polite request.
FORCE_FINAL_AFTER_TURN = 4

TOOL_SCHEMAS = [
    {"type": "function", "function": {
        "name": "kg_neighbors",
        "description": "The k-hop knowledge-graph neighbourhood of one gene (Reactome pathways, STRING/BioGRID partners, fusion partners, co-dependent genes), each with source/edge-type/score provenance.",
        "parameters": {"type": "object", "properties": {
            "gene": {"type": "string"}, "k": {"type": "integer", "default": 1},
        }, "required": ["gene"]},
    }},
    {"type": "function", "function": {
        "name": "shared_pathways",
        "description": "Reactome pathways that ALL the given genes are 1-hop neighbours of -- strong co-membership evidence.",
        "parameters": {"type": "object", "properties": {
            "genes": {"type": "array", "items": {"type": "string"}},
        }, "required": ["genes"]},
    }},
    {"type": "function", "function": {
        "name": "interaction_evidence",
        "description": "Direct graph edges between two specific genes (STRING score, BioGRID call, shared fusion partner, co-dependency).",
        "parameters": {"type": "object", "properties": {
            "gene_a": {"type": "string"}, "gene_b": {"type": "string"},
        }, "required": ["gene_a", "gene_b"]},
    }},
    {"type": "function", "function": {
        "name": "resolve_symbol",
        "description": "Resolve a gene symbol or Ensembl ID against gene_reference.csv before looking it up in the graph.",
        "parameters": {"type": "object", "properties": {"token": {"type": "string"}}, "required": ["token"]},
    }},
]

TOOL_FUNCTIONS = {
    "kg_neighbors": lambda args: kg_tools.kg_neighbors(args["gene"], args.get("k", 1)),
    "shared_pathways": lambda args: kg_tools.shared_pathways(args["genes"]),
    "interaction_evidence": lambda args: kg_tools.interaction_evidence(args["gene_a"], args["gene_b"]),
    "resolve_symbol": lambda args: kg_tools.resolve_symbol(args["token"]),
}

SYSTEM_PROMPT = """You are the query-expansion agent for a cell-line ranking tool. A researcher \
has supplied one or more inclusion genes. Your job is to suggest OTHER genes worth considering, \
using the kg_neighbors/shared_pathways/interaction_evidence/resolve_symbol tools to find real \
graph evidence -- never suggest a gene from memory or general biology knowledge.

You may call a tool more than once and chain calls (e.g. check a 1-hop neighbour's own \
connections before suggesting it), but budget your calls: **make at most 3-4 tool calls in total**, \
then stop and give your final answer using the strongest evidence you have already found -- do \
not keep exploring indefinitely looking for a "best" answer. When you are ready, give your FINAL \
answer as a JSON array, and nothing else, of objects shaped exactly like this:
[{"gene": "SYMBOL", "reason": "one sentence citing the evidence", "source_db": "string_v12|biogrid|reactome|fusions|dependency", "edge_type": "...", "score": <number or null>}]

These are SUGGESTIONS ONLY, presented for the researcher to accept or reject -- never state that \
a gene has been added to the query, because nothing here changes the query. If no tool call \
finds any related gene, return an empty JSON array []."""


class QueryExpansionFailedError(RuntimeError):
    pass


def expand(inclusion_genes, registry=None, max_turns=MAX_TOOL_TURNS):
    registry = registry or citations_module.build_registry()

    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": f"Inclusion genes: {', '.join(inclusion_genes)}"},
    ]
    tool_results_text = []

    for turn in range(max_turns):
        if turn == FORCE_FINAL_AFTER_TURN:
            messages.append({"role": "user", "content": (
                "You have used your tool-call budget. Give your final JSON array answer now, "
                "using only the evidence already gathered above -- no more tool calls."
            )})
        force_final = turn >= FORCE_FINAL_AFTER_TURN
        message = provider.chat(messages, tools=None if force_final else TOOL_SCHEMAS)
        tool_calls = getattr(message, "tool_calls", None)
        if not tool_calls and not message.content:
            # Some models occasionally return an empty turn (no tool call, no content) rather
            # than a real final answer -- a live, observed failure mode, not a hypothetical
            # one. Nudge once instead of treating silence as "done"; still bounded by
            # max_turns, so this cannot loop forever.
            messages.append({"role": "user", "content": (
                "Please give your final answer now: a JSON array of suggestions (or [] if none), "
                "and nothing else."
            )})
            continue
        if not tool_calls:
            try:
                suggestions = json.loads(message.content)
            except (json.JSONDecodeError, TypeError) as exc:
                raise QueryExpansionFailedError(f"final response was not valid JSON: {exc}")

            # Every suggestion is checked against R2's provenance requirement structurally
            # (source_db/edge_type present) and against R1/R3 via the verifier over the raw
            # tool-result text this conversation actually retrieved.
            missing_provenance = [
                s for s in suggestions if not s.get("source_db") or not s.get("edge_type")
            ]
            evidence_values = []
            for text in tool_results_text:
                for token in verifier.extract_numbers_raw(text):
                    try:
                        evidence_values.append(float(token.rstrip("%").replace(",", "")))
                    except ValueError:
                        continue
            reasons_text = "\n".join(s.get("reason", "") for s in suggestions)
            check = verifier.verify(reasons_text, {"values": evidence_values}, registry)

            return {
                "suggestions": suggestions,
                "missing_provenance": missing_provenance,
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
                result = TOOL_FUNCTIONS.get(name, lambda a: {"error": f"unknown tool {name!r}"})(args)
            except Exception as exc:
                # A graph-lookup failure (a malformed gene token, etc.) must not crash the
                # whole expansion -- the model sees the error and can retry or move on.
                result = {"error": f"{type(exc).__name__}: {exc}"}
            result_text = json.dumps(result, default=str)
            tool_results_text.append(result_text)
            messages.append({
                "role": "tool", "tool_call_id": tool_call.id, "name": name, "content": result_text,
            })

    raise QueryExpansionFailedError(f"no final answer within {max_turns} tool-call turns")
