import json

import literature
import provider

MAX_TOOL_TURNS = 6

SYSTEM_PROMPT = """You are the literature agent for an explainable cell-line selection tool. \
Given a question about supporting or contradicting evidence for a gene/biology claim, use the \
pubmed_search, pmc_fetch, cache_lookup and rerank tools to find real PubMed abstracts.

Hard rule, absolute: for every finding, you must give a VERBATIM QUOTE copied exactly from an \
abstract you actually fetched, plus that abstract's PMID. If you cannot find a passage you can \
quote exactly, do not report a finding for that claim at all -- do not paraphrase, summarise, or \
describe what an abstract "shows" in your own words. A finding without an exact quote is worse \
than no finding.

When you are done searching, respond with ONLY a JSON object matching this shape:
{"findings": [{"pmid": "...", "quote": "the exact substring, copied verbatim from the abstract", "relevance": "one sentence on why this bears on the question"}]}
If you found nothing quotable, respond with {"findings": []}."""

TOOL_SCHEMAS = [
    {
        "type": "function",
        "function": {
            "name": "pubmed_search",
            "description": "Search PubMed and return matching PMIDs.",
            "parameters": {
                "type": "object",
                "properties": {"terms": {"type": "string"}, "retmax": {"type": "integer", "default": 10}},
                "required": ["terms"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "pmc_fetch",
            "description": "Fetch the title and abstract text for one PMID.",
            "parameters": {
                "type": "object",
                "properties": {"pmid": {"type": "string"}},
                "required": ["pmid"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "cache_lookup",
            "description": "Check the local cache for a previously fetched search or abstract by its cache key.",
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
            "name": "rerank",
            "description": "Rerank a list of already-fetched {pmid, text} passages by relevance to a query.",
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {"type": "string"},
                    "passages": {"type": "array", "items": {"type": "object"}},
                },
                "required": ["query", "passages"],
            },
        },
    },
]

TOOL_FUNCTIONS = {
    "pubmed_search": lambda args: literature.pubmed_search(args["terms"], args.get("retmax", 10)),
    "pmc_fetch": lambda args: literature.pmc_fetch(args["pmid"]),
    "cache_lookup": lambda args: literature.cache_lookup(args["key"]),
    "rerank": lambda args: literature.rerank(args["query"], args["passages"]),
}


class LiteratureAgentFailedError(RuntimeError):
    pass


def _fetched_texts_by_pmid(tool_results):
    texts = {}
    for result in tool_results:
        if isinstance(result, dict) and result.get("pmid") and result.get("text"):
            texts[str(result["pmid"])] = result["text"]
    return texts


def find_evidence(question, max_turns=MAX_TOOL_TURNS):
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": question},
    ]
    tool_results = []
    json_repairs_left = 1

    for _ in range(max_turns):
        message = provider.chat(messages, tools=TOOL_SCHEMAS)
        tool_calls = getattr(message, "tool_calls", None)
        if not tool_calls and not message.content:
            # Some models occasionally return an empty turn (no tool call, no content) rather
            # than a real final answer -- same live, observed failure mode
            # query_expansion_agent.py already guards against. Nudge once instead of burning the
            # JSON-repair budget on it; still bounded by max_turns.
            messages.append({
                "role": "user",
                "content": "Please give your final answer now: the JSON object described in the "
                           "system prompt (or {\"findings\": []} if you found nothing), and "
                           "nothing else.",
            })
            continue
        if not tool_calls:
            try:
                parsed = json.loads(message.content)
                findings = parsed.get("findings", [])
            except (json.JSONDecodeError, TypeError, AttributeError) as exc:
                if json_repairs_left > 0:
                    # One repair chance, same convention as narrator.py's repair turn -- a
                    # malformed final answer is usually a format slip after a search that already
                    # succeeded, not a sign the search itself failed.
                    json_repairs_left -= 1
                    messages.append({"role": "assistant", "content": message.content or ""})
                    messages.append({
                        "role": "user",
                        "content": f"That response was not valid JSON ({exc}). Respond again with "
                                   "ONLY the JSON object described in the system prompt.",
                    })
                    continue
                raise LiteratureAgentFailedError(f"final response was not valid JSON: {exc}")

            fetched = _fetched_texts_by_pmid(tool_results)
            confirmed = [
                f for f in findings
                if str(f.get("pmid")) in fetched and f.get("quote", "") in fetched[str(f["pmid"])]
            ]
            return {
                "findings": confirmed,
                "n_reported_by_model": len(findings),
                "n_quote_confirmed": len(confirmed),
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
                # A flaky external call (NCBI rate limit, a bad PMID the model guessed) must
                # not crash the whole request -- the model sees the error and can retry with a
                # different argument or a different tool, same as a human researcher would.
                result = {"error": f"{type(exc).__name__}: {exc}"}
            tool_results.append(result)
            messages.append({
                "role": "tool", "tool_call_id": tool_call.id, "name": name,
                "content": json.dumps(result, default=str),
            })

    raise LiteratureAgentFailedError(f"no final answer within {max_turns} tool-call turns")
