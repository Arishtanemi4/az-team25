"""Bounded, reproducible evaluation summaries for the approved research extension.

This module consumes records captured from complete, already-run query workflows.  It never invokes
the scorer, changes similarity weights, opens processed data, or promotes an observed top result
to a biological claim.  The frozen query battery lives beside it and is loaded before any result
is inspected.  Real runs are written only under ``backend/research_service/runtime/evaluation``.
"""

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


EVALUATION_SCHEMA_VERSION = "bounded-research-evaluation-v1"
BATTERY_PATH = Path(__file__).with_name("evaluation_queries.json")
PARTITION_BUCKETS = (
    "ranked_cell_lines",
    "ranked_beyond_top_n",
    "low_confidence_lines",
    "insufficient_evidence_lines",
    "disqualified_lines",
)
MISSING_LAYER_STATES = {"not_assayed", "non_detected", "unavailable", "not_available", "error"}


def _canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _content_hash(value):
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def load_query_battery(path=None):
    """Load the frozen S13 battery and reject malformed/duplicate query definitions."""
    with open(Path(path) if path is not None else BATTERY_PATH, encoding="utf-8") as handle:
        battery = json.load(handle)
    if battery.get("schema_version") != "research-evaluation-battery-v1":
        raise ValueError("Unexpected evaluation query battery schema_version")
    cases = battery.get("cases")
    repeat = battery.get("lineage_repeat")
    if not isinstance(cases, list) or not isinstance(repeat, dict):
        raise ValueError("Evaluation query battery requires cases and lineage_repeat")
    all_cases = cases + [repeat]
    ids = [case.get("case_id") for case in all_cases]
    if any(not isinstance(case_id, str) or not case_id for case_id in ids) or len(ids) != len(set(ids)):
        raise ValueError("Evaluation query case_id values must be non-empty and unique")
    for case in all_cases:
        for field in ("inclusion_symbols", "exclusion_symbols"):
            if not isinstance(case.get(field), list) or not all(isinstance(value, str) and value for value in case[field]):
                raise ValueError(f"Evaluation case {case['case_id']!r} has invalid {field}")
    if battery.get("cohort_scope") != "unfiltered human cohort":
        raise ValueError("S13 battery must retain the unfiltered human cohort scope")
    return battery


def select_most_populated_lineage(metadata_rows):
    """Select the lineage-repeat filter from pinned metadata, with no name/alias inference."""
    counts = Counter()
    for row in metadata_rows or []:
        if not isinstance(row, dict):
            continue
        species = str(row.get("species", "human")).strip().lower()
        if species not in {"human", "homo sapiens"}:
            continue
        lineage = row.get("lineage")
        if lineage is None:
            continue
        lineage = str(lineage).strip()
        if not lineage or lineage.lower() in {"unknown", "none", "nan"}:
            continue
        counts[lineage] += 1
    if not counts:
        return {"status": "unavailable", "reason": "No explicit human, nonmissing lineage is available in metadata."}
    maximum = max(counts.values())
    selected = sorted(lineage for lineage, count in counts.items() if count == maximum)[0]
    return {
        "status": "ok",
        "lineage": selected,
        "model_count": counts[selected],
        "all_lineage_counts": dict(sorted(counts.items())),
    }


def resolved_battery_cases(battery, metadata_rows):
    """Return the fixed six unfiltered cases and the deterministic EGFR lineage repeat."""
    resolved = []
    for case in battery["cases"]:
        item = dict(case)
        item["filters"] = {}
        item["cohort_scope"] = battery["cohort_scope"]
        resolved.append(item)
    repeat = dict(battery["lineage_repeat"])
    lineage = select_most_populated_lineage(metadata_rows)
    repeat["cohort_scope"] = battery["cohort_scope"]
    repeat["lineage_selection"] = lineage
    repeat["filters"] = {"lineage": lineage["lineage"]} if lineage["status"] == "ok" else {}
    resolved.append(repeat)
    return resolved


def _iter_native_lines(native_result):
    for bucket in PARTITION_BUCKETS:
        for line in native_result.get(bucket, []) or []:
            if isinstance(line, dict):
                yield bucket, line


def native_signature(native_result):
    """The precise D/tier/eligible-order signature used for on/off regression checks."""
    eligible = []
    all_rows = []
    for bucket, line in _iter_native_lines(native_result):
        row = (bucket, line.get("model_id"), line.get("D"), line.get("confidence_tier"), line.get("veto"))
        all_rows.append(row)
        if bucket in {"ranked_cell_lines", "ranked_beyond_top_n"}:
            eligible.append(row)
    return {"all_partitions": all_rows, "eligible_order": eligible}


def layer_missingness(native_result):
    """Count source states without recoding measured absence or vetoed evidence as missing."""
    totals = {}
    for _, line in _iter_native_lines(native_result):
        for gene in line.get("per_gene", []) or []:
            for layer_id, layer in (gene.get("layers") or {}).items():
                state = layer.get("state") if isinstance(layer, dict) else None
                key = state if isinstance(state, str) and state else "unknown_state"
                entry = totals.setdefault(layer_id, {"total": 0, "missing": 0, "state_counts": Counter()})
                entry["total"] += 1
                entry["state_counts"][key] += 1
                if key in MISSING_LAYER_STATES:
                    entry["missing"] += 1
    result = {}
    for layer_id, entry in sorted(totals.items()):
        total = entry["total"]
        result[layer_id] = {
            "total": total,
            "missing": entry["missing"],
            "coverage_fraction": (total - entry["missing"]) / total if total else None,
            "state_counts": dict(sorted(entry["state_counts"].items())),
        }
    return result


def context_source_coverage(contexts):
    """Summarise packets actually captured; absence is not transformed into a negative result."""
    measured = 0
    graph_modes = Counter()
    for packet in (contexts or {}).values():
        if not isinstance(packet, dict):
            continue
        measured_packet = packet.get("measured_context", packet.get("measured"))
        if isinstance(measured_packet, dict):
            measured += 1
        graph = packet.get("graph_context", packet.get("graph"))
        source = graph.get("graph_source") if isinstance(graph, dict) else None
        if isinstance(source, dict):
            mode = source.get("mode") or source.get("status") or "unknown"
        else:
            mode = "not_recorded"
        graph_modes[str(mode)] += 1
    return {
        "context_models_recorded": measured,
        "graph_mode_counts": dict(sorted(graph_modes.items())),
    }


def _available_views(candidate, omitted=None, equal_weights=False):
    views = []
    for view_id, view in (candidate.get("view_scores") or {}).items():
        if view_id == omitted or not isinstance(view, dict) or view.get("status") != "available":
            continue
        score, weight = view.get("score"), view.get("weight")
        if not isinstance(score, (int, float)) or not math.isfinite(score):
            continue
        if not isinstance(weight, (int, float)) or weight <= 0:
            continue
        views.append((view_id, float(score), 1.0 if equal_weights else float(weight)))
    return views


def _rank_candidates(candidates, omitted=None, equal_weights=False, minimum_available_views=2):
    ranked = []
    for candidate in candidates or []:
        model_id = candidate.get("model_id") if isinstance(candidate, dict) else None
        views = _available_views(candidate, omitted=omitted, equal_weights=equal_weights)
        if not isinstance(model_id, str) or len(views) < minimum_available_views:
            continue
        total_weight = sum(view[2] for view in views)
        score = sum(view[1] * view[2] for view in views) / total_weight
        ranked.append((model_id, score))
    return sorted(ranked, key=lambda item: (-item[1], item[0]))


def _top_overlap(baseline, alternative):
    base_ids = [model_id for model_id, _ in baseline[:3]]
    alt_ids = [model_id for model_id, _ in alternative[:3]]
    if not base_ids or not alt_ids:
        return {
            "status": "undefined",
            "reason": "No eligible alternatives were available in the baseline or sensitivity scenario.",
            "baseline_top3": base_ids,
            "scenario_top3": alt_ids,
            "overlap_count": None,
            "jaccard": None,
        }
    intersection = set(base_ids) & set(alt_ids)
    union = set(base_ids) | set(alt_ids)
    return {
        "status": "ok",
        "baseline_top3": base_ids,
        "scenario_top3": alt_ids,
        "overlap_count": len(intersection),
        "jaccard": len(intersection) / len(union),
    }


def alternative_sensitivity(candidates, minimum_available_views=2):
    """Compare pre-truncation C5 candidates under fixed equal-weight and leave-view-out variants.

    Callers must supply the complete eligible alternative set from a run, not only the displayed
    top three.  This keeps the reported top-three overlap honest.
    """
    baseline = _rank_candidates(candidates, minimum_available_views=minimum_available_views)
    all_views = sorted({view_id for candidate in candidates or [] for view_id in (candidate.get("view_scores") or {})})
    scenarios = {
        "equal_view_weights": _top_overlap(
            baseline,
            _rank_candidates(candidates, equal_weights=True, minimum_available_views=minimum_available_views),
        )
    }
    for view_id in all_views:
        scenarios[f"omit_{view_id}"] = _top_overlap(
            baseline,
            _rank_candidates(candidates, omitted=view_id, minimum_available_views=minimum_available_views),
        )
    return {
        "status": "ok" if baseline else "undefined",
        "reason": None if baseline else "No eligible alternatives were supplied for this anchor.",
        "baseline_top3": [model_id for model_id, _ in baseline[:3]],
        "scenarios": scenarios,
    }


def top_k_invariance(observations):
    """Check the S13 top_k=1/10/20 invariant from pre-recorded full-query observations."""
    expected = {1, 10, 20}
    by_top_k = {item.get("top_k"): item for item in observations or [] if isinstance(item, dict)}
    if not expected.issubset(by_top_k):
        return {"status": "not_run", "reason": "top_k observations for 1, 10, and 20 were not all supplied."}
    fingerprints = [
        (by_top_k[value].get("profile_hash"), tuple(by_top_k[value].get("full_backup_model_ids") or []))
        for value in sorted(expected)
    ]
    return {
        "status": "pass" if len(set(fingerprints)) == 1 else "fail",
        "observations": {str(value): {"profile_hash": by_top_k[value].get("profile_hash"), "full_backup_model_ids": by_top_k[value].get("full_backup_model_ids", [])} for value in sorted(expected)},
    }


def _regression_check(native_result, variants):
    baseline = native_signature(native_result)
    if not variants:
        return {"status": "not_run", "reason": "No extension on/off native exports were supplied."}
    mismatches = [name for name, value in variants.items() if native_signature(value) != baseline]
    return {"status": "pass" if not mismatches else "fail", "mismatched_variants": mismatches}


def _alternative_summary(result):
    """Keep alternative similarity, reference and joint evidence separate from native scores."""
    top_three = []
    for item in result.get("top_alternatives", []) or []:
        joint = {
            view_id: {
                "joint_gene_count": view.get("joint_gene_count"),
                "joint_gene_fraction": view.get("joint_gene_fraction"),
                "status": view.get("status"),
            }
            for view_id, view in sorted((item.get("view_scores") or {}).items())
            if isinstance(view, dict)
        }
        top_three.append({
            "model_id": item.get("model_id"),
            "S": item.get("S"),
            "available_weight_fraction": item.get("available_weight_fraction"),
            "joint_by_view": joint,
        })
    return {
        "status": result.get("status"),
        "candidate_counts": result.get("candidate_counts"),
        "reference_count": (result.get("reference_cohort") or {}).get("count"),
        "top3": top_three,
        "legacy_rna_spearman_available": bool(result.get("legacy_native_similar_lines")),
    }


def evaluate_record(record):
    """Summarise one already-executed case without re-running or interpreting its score."""
    native_result = record.get("native_result") or {}
    alternatives = record.get("alternatives_by_anchor") or {}
    sensitivity_inputs = record.get("sensitivity_candidates_by_anchor") or {}
    return {
        "case_id": record["case_id"],
        "status": record.get("status", "completed"),
        "runtime_seconds": record.get("runtime_seconds"),
        "candidate_counts": native_result.get("candidate_counts"),
        "layer_missingness": layer_missingness(native_result),
        "context_source_coverage": context_source_coverage(record.get("context_by_model")),
        "alternatives": {
            anchor_id: _alternative_summary(result)
            for anchor_id, result in sorted(alternatives.items()) if isinstance(result, dict)
        },
        "sensitivity": {
            anchor_id: alternative_sensitivity(candidates)
            for anchor_id, candidates in sorted(sensitivity_inputs.items())
        },
        "top_k_invariance": top_k_invariance(record.get("top_k_observations")),
        "extension_regression": _regression_check(native_result, record.get("extension_variants")),
    }


def evaluate_records(records, battery, metadata_rows=None):
    """Produce an auditable result with explicit not-run entries for every fixed battery case."""
    resolved = resolved_battery_cases(battery, metadata_rows or [])
    case_ids = {case["case_id"] for case in resolved}
    supplied = {}
    for record in records or []:
        case_id = record.get("case_id") if isinstance(record, dict) else None
        if case_id not in case_ids:
            raise ValueError(f"Evaluation result has an unknown case_id: {case_id!r}")
        if case_id in supplied:
            raise ValueError(f"Evaluation result is duplicated for case_id: {case_id!r}")
        supplied[case_id] = evaluate_record(record)
    outcomes = []
    for case in resolved:
        outcome = supplied.get(case["case_id"])
        if outcome is None:
            outcome = {"case_id": case["case_id"], "status": "not_run", "reason": "No complete-pipeline record was supplied for this fixed case."}
        outcomes.append({"query": case, "outcome": outcome})
    complete = all(item["outcome"].get("status") == "completed" for item in outcomes)
    return {
        "schema_version": EVALUATION_SCHEMA_VERSION,
        "query_battery_hash": _content_hash(battery),
        "status": "complete" if complete else "incomplete",
        "cases": outcomes,
        "limitations": [
            "This evaluation reports observed software outputs and coverage only; it does not establish laboratory success or biological superiority.",
            "Sensitivity rankings require complete pre-truncation eligible alternatives per anchor; displayed top-three values alone are insufficient.",
        ],
    }


def write_runtime_evaluation(result, manifest, output_dir):
    """Write rebuildable local outputs only; committed reports remain human-reviewed records."""
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    for filename, payload in (("evaluation.json", result), ("experiment_manifest.json", manifest)):
        with open(directory / filename, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Summarise a predeclared S13 research evaluation record.")
    parser.add_argument("--input", required=True, help="JSON with metadata_rows, records, and optional manifest fields.")
    parser.add_argument("--output-dir", required=True, help="Ignored backend/research_service/runtime/evaluation destination.")
    parser.add_argument("--battery", default=str(BATTERY_PATH), help="Frozen evaluation query battery JSON.")
    args = parser.parse_args(argv)
    with open(args.input, encoding="utf-8") as handle:
        source = json.load(handle)
    battery = load_query_battery(args.battery)
    result = evaluate_records(source.get("records"), battery, source.get("metadata_rows"))
    manifest = dict(source.get("manifest") or {})
    manifest.update({
        "schema_version": "research-experiment-manifest-v1",
        "query_battery_hash": result["query_battery_hash"],
        "evaluation_schema_version": EVALUATION_SCHEMA_VERSION,
        "result_status": result["status"],
    })
    write_runtime_evaluation(result, manifest, args.output_dir)


if __name__ == "__main__":
    main()
