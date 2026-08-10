import json
import math

import numpy as np

BOUNDARY_STATEMENT = (
    "High confidence means the available evidence is consistent. It does not mean the "
    "result has been experimentally proven. Calibrated probabilities require held-out "
    "validation that does not yet exist."
)

SCHEMA_VERSION = "rank-result-v2"


def _export_diagnostic(diagnostic):
    if diagnostic is None:
        return None

    exported = dict(diagnostic)
    veto_counts = diagnostic.get("veto_counts", {})
    veto_counts_by_gene = []
    veto_counts_by_ensembl = {}

    for key, count in veto_counts.items():
        if isinstance(key, tuple) and len(key) == 2:
            ensembl_id, symbol = key
        else:
            # Accept an already JSON-safe mapping from a caller without guessing a symbol.
            ensembl_id, symbol = str(key), None
        veto_counts_by_gene.append(
            {"ensembl_id": ensembl_id, "symbol": symbol, "count": count}
        )
        veto_counts_by_ensembl[ensembl_id] = veto_counts_by_ensembl.get(ensembl_id, 0) + count

    veto_counts_by_gene.sort(key=lambda item: (item["ensembl_id"], item["symbol"] or ""))
    exported["veto_counts"] = veto_counts_by_ensembl
    exported["veto_counts_by_gene"] = veto_counts_by_gene
    return exported


def _candidate_counts(ranked, ranked_beyond_top_n, low_confidence, insufficient, disqualified):
    eligible = len(ranked) + len(ranked_beyond_top_n)
    return {
        "evaluated": eligible + len(low_confidence) + len(insufficient) + len(disqualified),
        "eligible": eligible,
        "displayed": len(ranked),
        "ranked_beyond_top_n": len(ranked_beyond_top_n),
        "low_confidence": len(low_confidence),
        "insufficient": len(insufficient),
        "disqualified": len(disqualified),
    }


def export_query_result(query_result, inclusion_tokens, exclusion_tokens, filters=None):
    ranked = query_result["ranked"]
    insufficient = query_result["insufficient"]
    low_confidence = query_result.get("low_confidence", [])
    ranked_beyond_top_n = query_result.get("ranked_beyond_top_n", [])
    disqualified = query_result.get("disqualified", [])
    counts = _candidate_counts(
        ranked, ranked_beyond_top_n, low_confidence, insufficient, disqualified
    )

    return {
        "schema_version": SCHEMA_VERSION,
        "query": {
            "inclusion_genes": list(inclusion_tokens),
            "exclusion_genes": list(exclusion_tokens),
            "filters": dict(filters or {}),
        },
        "ranked_cell_lines": ranked,
        "ranked_beyond_top_n": ranked_beyond_top_n,
        "low_confidence_lines": low_confidence,
        "insufficient_evidence_lines": insufficient,
        "disqualified_lines": disqualified,
        "total_ranked": query_result.get("total_ranked", counts["eligible"]),
        "candidate_counts": counts,
        "diagnostic": _export_diagnostic(query_result["diagnostic"]),
        "boundary_statement": BOUNDARY_STATEMENT,
    }


def sanitize_for_json(obj):
    if isinstance(obj, dict):
        return {key: sanitize_for_json(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [sanitize_for_json(value) for value in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, float):
        return None if not math.isfinite(obj) else obj
    return obj


_sanitize = sanitize_for_json


def write_export_json(export_dict, path):
    with open(path, "w") as f:
        json.dump(sanitize_for_json(export_dict), f, indent=2, allow_nan=False)
