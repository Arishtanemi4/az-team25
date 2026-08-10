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


SHUFFLE_DRAWS = 20
SHUFFLE_SEED_BASE = 20260813  # disclosed, fixed -- PARAMETERS.md

SHUFFLE_RBO_CEILING = 0.90


def find_override_informative_queries(battery, rna_constants, gene_reference_df):
    per_lineage = rna_constants.get("per_lineage", {})
    informative = []
    for query in battery:
        lineage = query.get("lineage")
        if lineage is None or lineage not in per_lineage:
            continue
        override_genes = per_lineage[lineage]
        resolved, _ambiguous, _unresolved = score.resolve_genes(
            query["inclusion"] + query["exclusion"], gene_reference_df
        )
        if any(eid in override_genes for eid in resolved.values()):
            informative.append(query["name"])
    return informative


def shuffle_lineage_labels(cell_lines_by_id, seed):
    rng = random.Random(seed)
    model_ids = list(cell_lines_by_id.keys())
    lineages = [cell_lines_by_id[m].get("lineage") for m in model_ids]
    rng.shuffle(lineages)
    shuffled = {}
    for model_id, lineage in zip(model_ids, lineages):
        row = dict(cell_lines_by_id[model_id])
        row["lineage"] = lineage
        shuffled[model_id] = row
    return shuffled


def run_lineage_shuffle_control(query_cache, full_cell_lines_by_id, K=SHUFFLE_DRAWS,
                                 seed_base=SHUFFLE_SEED_BASE, top_n=10):
    weights, thresholds = sensitivity.default_sample()
    tier_params = sensitivity.tier_params_from_thresholds(thresholds)
    rna_constants = query_cache["tables"]["rna_constants"]
    extended_constants = query_cache["tables"]["extended_constants"]
    essentiality_constants = query_cache["tables"]["essentiality_constants"]

    real_top10, _real_full = sensitivity.run_query_under_sample(
        query_cache, weights, tier_params, rna_constants, extended_constants,
        essentiality_constants, top_n=top_n,
    )
    real_rho_bar = query_cache["correlation_weights"]["rho_bar"]
    real_m_eff = query_cache["correlation_weights"]["m_eff"]

    rbo_draws, membership_draws = [], []
    rho_bar_invariant, m_eff_invariant = True, True
    for i in range(K):
        shuffled_full = shuffle_lineage_labels(full_cell_lines_by_id, seed_base + i)
        shuffled = {m: shuffled_full[m] for m in query_cache["model_ids"]}
        top10, _full = sensitivity.run_query_under_sample(
            query_cache, weights, tier_params, rna_constants, extended_constants,
            essentiality_constants, top_n=top_n, cell_lines_by_id_override=shuffled,
        )
        rbo_draws.append(sensitivity.rbo_at_10(real_top10, top10))
        _changed, n_diff = sensitivity.top10_membership_change(real_top10, top10)
        membership_draws.append(n_diff)
        if query_cache["correlation_weights"]["rho_bar"] != real_rho_bar:
            rho_bar_invariant = False
        if query_cache["correlation_weights"]["m_eff"] != real_m_eff:
            m_eff_invariant = False

    valid_rbo = [r for r in rbo_draws if r is not None]
    mean_rbo = float(np.mean(valid_rbo)) if valid_rbo else None
    passed = mean_rbo is not None and mean_rbo <= SHUFFLE_RBO_CEILING

    return {
        "query_name": query_cache["name"],
        "k_draws": K,
        "rbo_at_10_per_draw": rbo_draws,
        "mean_rbo_at_10": mean_rbo,
        "membership_change_per_draw": membership_draws,
        "rho_bar_invariant": rho_bar_invariant,
        "m_eff_invariant": m_eff_invariant,
        "pass": passed,
    }


def run_shuffled_lineage_null(battery, gene_reference_df, cell_lines_df, coverage_df,
                               data_dir=DATA_DIR, K=SHUFFLE_DRAWS, layer_frames=None):
    if layer_frames is None:
        battery_genes = sensitivity.resolve_battery_genes(battery, gene_reference_df)
        layer_frames = sensitivity.load_layer_frames(battery_genes, data_dir)

    sample_cache = sensitivity.build_query_cache(
        battery[0], gene_reference_df, cell_lines_df, coverage_df, layer_frames, data_dir
    )
    rna_constants = sample_cache["tables"]["rna_constants"]
    informative_names = set(
        find_override_informative_queries(battery, rna_constants, gene_reference_df)
    )
    full_cell_lines_by_id = cell_lines_df.set_index("ModelID").to_dict(orient="index")

    per_query_results = []
    for query in battery:
        if query["name"] not in informative_names:
            continue
        query_cache = sensitivity.build_query_cache(
            query, gene_reference_df, cell_lines_df, coverage_df, layer_frames, data_dir
        )
        per_query_results.append(
            run_lineage_shuffle_control(query_cache, full_cell_lines_by_id, K=K)
        )

    n_informative = len(per_query_results)
    n_passed = sum(1 for r in per_query_results if r["pass"])
    passed = n_informative > 0 and n_passed >= (n_informative / 2)

    return {
        "informative_queries": sorted(informative_names),
        "per_query": per_query_results,
        "n_informative": n_informative,
        "n_passed": n_passed,
        "pass": passed,
    }

LOO_TOP_K_LINEAGES = 5
LOO_MIN_N = 100
LOO_MEMBERSHIP_CEILING = 2  # of 10, disclosed heuristic


def select_holdout_lineages(cell_lines_df, top_k=LOO_TOP_K_LINEAGES, min_n=LOO_MIN_N):
    counts = cell_lines_df["lineage"].value_counts()
    eligible = counts[counts >= min_n]
    return eligible.head(top_k).index.tolist()
