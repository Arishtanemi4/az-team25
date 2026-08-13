import json
from pathlib import Path

import citations as citations_module
import docs_lookup
import provider
import verifier

REPO_ROOT = Path(__file__).resolve().parent.parent
GOLD_SET_PATH = Path(__file__).resolve().parent / "resources" / "gold_set.json"
REPORT_DIR = Path(__file__).resolve().parent / "resources" / "eval_reports"

REFUSAL_PROBES = [
    {
        "id": "refusal-unknown-gene",
        "question": "What is the RNA evidence for the gene ZZZFAKE1 in this project's data?",
        "probe_type": "unknown_gene",
    },
    {
        "id": "refusal-out-of-panel-cell-line",
        "question": "What is the confidence tier for cell line ACH-999999?",
        "probe_type": "out_of_panel_cell_line",
    },
    {
        "id": "refusal-nonexistent-parameter",
        "question": "What is the value of the FOO_BAR_THRESHOLD parameter?",
        "probe_type": "nonexistent_parameter",
    },
]

# A keyword heuristic, not a model judgement -- deliberately simple and auditable. A real
# abstention says one of these things; a fabricated answer essentially never does.
_ABSTENTION_KEYWORDS = (
    "no data", "not found", "no evidence", "does not exist", "doesn't exist", "no such",
    "cannot find", "can't find", "unknown parameter", "not a registered", "not available",
    "no record", "not assayed", "abstain", "not in the panel", "not in this project",
    "i do not have", "i don't have", "no information",
)


def load_gold_set(path=GOLD_SET_PATH):
    return json.loads(path.read_text(encoding="utf-8"))


def gate_numeric_fidelity(items_with_text, evidence_record):
    results = []
    for item_id, text in items_with_text:
        failures = verifier.verify_numeric_fidelity(text, evidence_record)
        results.append({"id": item_id, "passed": not failures, "failures": failures})
    n_passed = sum(1 for r in results if r["passed"])
    return {
        "gate": "numeric_fidelity",
        "n_items": len(results),
        "n_passed": n_passed,
        "fraction_passed": (n_passed / len(results)) if results else None,
        "results": results,
    }


def gate_citation_resolution(items_with_text, registry=None):
    if registry is None:
        registry = citations_module.build_registry()
    results = []
    for item_id, text in items_with_text:
        failures = verifier.verify_citations(text, registry)
        results.append({"id": item_id, "passed": not failures, "failures": failures})
    n_passed = sum(1 for r in results if r["passed"])
    return {
        "gate": "citation_resolution",
        "n_items": len(results),
        "n_passed": n_passed,
        "fraction_passed": (n_passed / len(results)) if results else None,
        "results": results,
    }


def _looks_like_abstention(text):
    lowered = text.lower()
    return any(keyword in lowered for keyword in _ABSTENTION_KEYWORDS)


def gate_refusal_behaviour(answer_fn):
    results = []
    for probe in REFUSAL_PROBES:
        text = answer_fn(probe["question"])
        results.append({
            "id": probe["id"],
            "probe_type": probe["probe_type"],
            "abstained": _looks_like_abstention(text),
            "answer": text,
        })
    n_passed = sum(1 for r in results if r["abstained"])
    return {
        "gate": "refusal_behaviour",
        "n_items": len(results),
        "n_passed": n_passed,
        "fraction_passed": n_passed / len(results),
        "results": results,
    }


def _llm_judge_claims(prompt_text):
    messages = [
        {"role": "system", "content": "You are a strict fact-checking judge. Respond with only valid JSON, no markdown fences, no prose."},
        {"role": "user", "content": prompt_text},
    ]
    reply = provider.chat(messages)
    try:
        return json.loads(reply.content)
    except (json.JSONDecodeError, TypeError):
        messages.append({"role": "user", "content": "That was not valid JSON. Reply again with ONLY the JSON object, nothing else."})
        retry = provider.chat(messages)
        return json.loads(retry.content)


def faithfulness(answer_text, context_text):
    prompt = (
        "Decompose the ANSWER into short atomic factual claims. For each claim, decide whether "
        "it is directly supported by the CONTEXT.\n\n"
        'Respond with exactly this JSON shape: {"claims": [{"claim": "...", "supported": true|false}]}\n\n'
        f"CONTEXT:\n{context_text}\n\nANSWER:\n{answer_text}"
    )
    result = _llm_judge_claims(prompt)
    claims = result.get("claims", [])
    if not claims:
        return None
    return sum(1 for c in claims if c.get("supported")) / len(claims)


def context_recall(gold_answer, context_text):
    prompt = (
        "Decompose the GOLD ANSWER into short atomic factual claims. For each claim, decide "
        "whether it can be attributed to (found in, or directly inferable from) the CONTEXT.\n\n"
        'Respond with exactly this JSON shape: {"claims": [{"claim": "...", "attributable": true|false}]}\n\n'
        f"CONTEXT:\n{context_text}\n\nGOLD ANSWER:\n{gold_answer}"
    )
    result = _llm_judge_claims(prompt)
    claims = result.get("claims", [])
    if not claims:
        return None
    return sum(1 for c in claims if c.get("attributable")) / len(claims)


def gate_ragas_metrics(gold_set, answer_fn, context_fn, sample_size=None):
    items = gold_set if sample_size is None else gold_set[:sample_size]
    results = []
    for item in items:
        answer = answer_fn(item["question"])
        context = context_fn(item["question"])
        results.append({
            "id": item["id"],
            "faithfulness": faithfulness(answer, context),
            "context_recall": context_recall(item["answer"], context),
        })
    faith_scores = [r["faithfulness"] for r in results if r["faithfulness"] is not None]
    recall_scores = [r["context_recall"] for r in results if r["context_recall"] is not None]
    return {
        "gate": "ragas_metrics",
        "n_items": len(results),
        "mean_faithfulness": sum(faith_scores) / len(faith_scores) if faith_scores else None,
        "mean_context_recall": sum(recall_scores) / len(recall_scores) if recall_scores else None,
        "results": results,
    }


def _dry_run_answer_fn(question, gold_by_question):
    return gold_by_question.get(question, "No data available for this question.")


def run(answer_fn=None, context_fn=None, registry=None, ragas_sample_size=10):
    gold_set = load_gold_set()
    registry = registry or citations_module.build_registry()

    dry_run = answer_fn is None
    if dry_run:
        gold_by_question = {item["question"]: item["answer"] for item in gold_set}
        answer_fn = lambda q: _dry_run_answer_fn(q, gold_by_question)

    numeric_items = [
        (item["id"], answer_fn(item["question"]))
        for item in gold_set if item.get("contains_number")
    ]
    # The gold set draws from all three of PARAMETERS.md/ALGORITHM_SPEC.md/METHOD_DECISION.md
    # (RAG_ARCHITECTURE.md SS7), so a "parameter"-category answer's number and a
    # "method_decision"-category answer's number are graded against the same combined fact-set.
    # docs_lookup's default source_files already covers all three.
    known_facts = docs_lookup.flatten_all_parameter_values()
    numeric_report = gate_numeric_fidelity(numeric_items, known_facts)

    citation_items = [(item["id"], answer_fn(item["question"])) for item in gold_set]
    citation_report = gate_citation_resolution(citation_items, registry)

    refusal_report = gate_refusal_behaviour(answer_fn)

    report = {
        "dry_run": dry_run,
        "numeric_fidelity": numeric_report,
        "citation_resolution": citation_report,
        "refusal_behaviour": refusal_report,
    }

    if not dry_run and context_fn is not None:
        report["ragas_metrics"] = gate_ragas_metrics(gold_set, answer_fn, context_fn, ragas_sample_size)

    return report


def write_report(report, out_dir=REPORT_DIR):
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "latest_report.json"
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    lines = ["# RAG evaluation report", ""]
    lines.append(f"Dry run: {report['dry_run']}")
    for gate_name in ("numeric_fidelity", "citation_resolution", "refusal_behaviour"):
        gate = report[gate_name]
        lines.append(f"- **{gate_name}**: {gate['n_passed']}/{gate['n_items']} passed")
    if "ragas_metrics" in report:
        rm = report["ragas_metrics"]
        lines.append(f"- **mean faithfulness**: {rm['mean_faithfulness']}")
        lines.append(f"- **mean context recall**: {rm['mean_context_recall']}")
    md_path = out_dir / "latest_report.md"
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, md_path


if __name__ == "__main__":
    result = run()
    paths = write_report(result)
    print(f"Numeric fidelity: {result['numeric_fidelity']['n_passed']}/{result['numeric_fidelity']['n_items']}")
    print(f"Citation resolution: {result['citation_resolution']['n_passed']}/{result['citation_resolution']['n_items']}")
    print(f"Refusal behaviour: {result['refusal_behaviour']['n_passed']}/{result['refusal_behaviour']['n_items']}")
    print(f"Report written to {paths[0]} / {paths[1]}")
