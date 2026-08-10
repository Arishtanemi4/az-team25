import json
import random

import numpy as np

import score
import sensitivity

DATA_DIR = "data/processed"
RESOURCES_DIR = "scoring/resources"

KNOWN_NEGATIVE_PAIR_QUERY_NAME = "stress_egfr_kras_joint_inclusion"

KNOWN_NEGATIVE_PAIR_QUERY_LUNG = {
    "name": "stress_egfr_kras_joint_inclusion_lung",
    "category": "known_negative_pair_lung_restricted",
    "inclusion": ["EGFR", "KRAS"],
    "exclusion": [],
    "lineage": "lung",
}
CONFIDENT_TIER_CEILING = 0.10  # disclosed heuristic, PARAMETERS.md
VETO_MAJORITY_SYMBOLS = {"EGFR", "KRAS"}


def run_known_negative_pair_control(query, gene_reference_df, cell_lines_df, coverage_df,
                                     data_dir=DATA_DIR, layer_frames=None):
    if layer_frames is None:
        battery_genes = sensitivity.resolve_battery_genes([query], gene_reference_df)
        layer_frames = sensitivity.load_layer_frames(battery_genes, data_dir)
    query_cache = sensitivity.build_query_cache(
        query, gene_reference_df, cell_lines_df, coverage_df, layer_frames, data_dir
    )

    weights, thresholds = sensitivity.default_sample()
    tier_params = sensitivity.tier_params_from_thresholds(thresholds)
    tables = dict(query_cache["tables"])

    results = [
        score.score_one_line(
            model_id, query_cache["inclusion_genes"], query_cache["exclusion_genes"],
            tables, query_cache["correlation_weights"], query_cache["cell_lines_by_id"][model_id],
            layer_weights=weights, tier_params=tier_params, build_narrative=False,
        )
        for model_id in query_cache["model_ids"]
    ]

    total = len(results)
    tier_counts = {}
    for r in results:
        tier_counts[r["confidence_tier"]] = tier_counts.get(r["confidence_tier"], 0) + 1
    n_confident = tier_counts.get("High", 0) + tier_counts.get("Moderate", 0)
    confident_fraction = n_confident / total if total else 0.0

    non_insufficient = [r for r in results if r["confidence_tier"] != "Insufficient"]
    veto_named = [
        r for r in non_insufficient
        if r["veto"] is not None and r["veto"]["symbol"] in VETO_MAJORITY_SYMBOLS
    ]
    veto_fraction = len(veto_named) / len(non_insufficient) if non_insufficient else 0.0

    veto_summary = {}
    for r in results:
        if r["veto"] is not None:
            key = r["veto"]["symbol"]
            veto_summary[key] = veto_summary.get(key, 0) + 1

    passed = confident_fraction <= CONFIDENT_TIER_CEILING and veto_fraction > 0.5

    return {
        "query_name": query["name"],
        "total_candidates": total,
        "tier_counts": tier_counts,
        "n_confident": n_confident,
        "confident_fraction": confident_fraction,
        "non_insufficient_count": len(non_insufficient),
        "veto_summary": veto_summary,
        "veto_fraction_of_non_insufficient": veto_fraction,
        "pass": passed,
    }
