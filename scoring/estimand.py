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


RANDOM_BATTERY_SEED = 20260812  # disclosed, fixed -- PARAMETERS.md
RANDOM_BATTERY_REPS = 5
D_HIGH_THRESHOLD = 0.7

RANDOM_MEDIAN_D_CEILING = 0.2  # disclosed heuristic: shadow median D must sit near the floor
RANDOM_MEAN_D_RATIO_CEILING = 0.5  # shadow mean D must be < half the real battery's own mean D


def build_shadow_battery(real_battery, gene_reference_df, exclude=None,
                          n_reps=RANDOM_BATTERY_REPS, seed=RANDOM_BATTERY_SEED):
    exclude = set(exclude or [])
    symbol_counts = gene_reference_df["symbol"].value_counts()
    unique_symbols = set(symbol_counts[symbol_counts == 1].index)
    pool_df = gene_reference_df[
        gene_reference_df["symbol"].isin(unique_symbols)
        & ~gene_reference_df["ensembl_id"].isin(exclude)
    ]
    pool = pool_df["symbol"].dropna().tolist()

    rng = random.Random(seed)
    shadow_queries = []
    for query in real_battery:
        n_inclusion = len(query["inclusion"])
        n_exclusion = len(query["exclusion"])
        n_needed = n_inclusion + n_exclusion
        for rep in range(n_reps):
            draw = rng.sample(pool, n_needed)
            shadow_queries.append({
                "name": f"shadow__{query['name']}__rep{rep}",
                "category": "random_gene_list_control",
                "inclusion": draw[:n_inclusion],
                "exclusion": draw[n_inclusion:],
                "lineage": query.get("lineage"),
                "source_query": query["name"],
                "rep": rep,
            })
    return shadow_queries


def _d_summary(values):
    if not values:
        return {"n": 0, "mean": None, "median": None, "fraction_d_ge_0_7": None}
    arr = np.array(values)
    return {
        "n": len(arr),
        "mean": float(arr.mean()),
        "median": float(np.median(arr)),
        "fraction_d_ge_0_7": float((arr >= D_HIGH_THRESHOLD).mean()),
    }


def run_random_control(shadow_battery, real_baseline_results, gene_reference_df, cell_lines_df,
                        coverage_df, data_dir=DATA_DIR):
    battery_genes = sensitivity.resolve_battery_genes(shadow_battery, gene_reference_df)
    layer_frames = sensitivity.load_layer_frames(battery_genes, data_dir)
    weights, thresholds = sensitivity.default_sample()
    tier_params = sensitivity.tier_params_from_thresholds(thresholds)

    tier_counts = {}
    d_values = []
    for query in shadow_battery:
        query_cache = sensitivity.build_query_cache(
            query, gene_reference_df, cell_lines_df, coverage_df, layer_frames, data_dir
        )
        tables = dict(query_cache["tables"])
        for model_id in query_cache["model_ids"]:
            result = score.score_one_line(
                model_id, query_cache["inclusion_genes"], query_cache["exclusion_genes"],
                tables, query_cache["correlation_weights"], query_cache["cell_lines_by_id"][model_id],
                layer_weights=weights, tier_params=tier_params, build_narrative=False,
            )
            tier_counts[result["confidence_tier"]] = tier_counts.get(result["confidence_tier"], 0) + 1
            if result["D"] is not None:
                d_values.append(result["D"])

    total = sum(tier_counts.values())
    collapse_fraction = (
        (tier_counts.get("Low", 0) + tier_counts.get("Insufficient", 0)) / total if total else 0.0
    )

    real_tier_counts = {}
    real_d_values = []
    for r in real_baseline_results:
        real_tier_counts[r["confidence_tier"]] = real_tier_counts.get(r["confidence_tier"], 0) + 1
        if r["D"] is not None:
            real_d_values.append(r["D"])

    d_summary = _d_summary(d_values)
    real_d_summary = _d_summary(real_d_values)
    passed = (
        d_summary["median"] is not None and d_summary["median"] <= RANDOM_MEDIAN_D_CEILING
        and real_d_summary["mean"] is not None
        and d_summary["mean"] < RANDOM_MEAN_D_RATIO_CEILING * real_d_summary["mean"]
    )

    return {
        "n_shadow_queries": len(shadow_battery),
        "tier_counts": tier_counts,
        "total_scored": total,
        "low_insufficient_collapse_fraction": collapse_fraction,
        "D_summary": d_summary,
        "real_tier_counts": real_tier_counts,
        "real_total_scored": sum(real_tier_counts.values()),
        "real_D_summary": real_d_summary,
        "pass": passed,
    }
