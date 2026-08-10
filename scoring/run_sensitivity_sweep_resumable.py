import json
import os
import time

import numpy as np
import pandas as pd

import sensitivity

SEED = 20260810
N_SOBOL = 8
N_DIRICHLET = 32
CHECKPOINT_PATH = "scoring/resources/.sensitivity_sweep_checkpoint.json"  # gitignored, ephemeral
FINAL_OUT_PATH = f"{sensitivity.RESOURCES_DIR}/sensitivity_sweep_results.json"


def load_checkpoint():
    if os.path.exists(CHECKPOINT_PATH):
        with open(CHECKPOINT_PATH) as f:
            return json.load(f)
    return {"baseline": None, "dirichlet": [], "sobol": [], "grid": []}


def save_checkpoint(ckpt):

    tmp = CHECKPOINT_PATH + ".tmp"
    with open(tmp, "w") as f:
        json.dump(ckpt, f)
    os.replace(tmp, CHECKPOINT_PATH)


def main():
    t0 = time.time()
    gene_reference = pd.read_csv(f"{sensitivity.DATA_DIR}/gene_reference.csv")
    cell_lines = pd.read_csv(f"{sensitivity.DATA_DIR}/cell_lines.csv")
    coverage = pd.read_csv(f"{sensitivity.DATA_DIR}/coverage.csv")
    battery = sensitivity.load_battery()
    print(f"[{time.time()-t0:.1f}s] reference tables + battery loaded", flush=True)

    battery_genes = sensitivity.resolve_battery_genes(battery, gene_reference)
    layer_frames = sensitivity.load_layer_frames(battery_genes)
    query_caches = [
        sensitivity.build_query_cache(q, gene_reference, cell_lines, coverage, layer_frames)
        for q in battery
    ]
    print(f"[{time.time()-t0:.1f}s] Tier-1 caches built ({len(query_caches)} queries)", flush=True)

    rna_floor = sensitivity.restrict_rna_floor_to_genes(
        sensitivity._load_json(f"{sensitivity.RESOURCES_DIR}/desirability_constants_rna_sweep_floor.json"),
        battery_genes,
    )
    extended_floor = sensitivity.restrict_extended_floor_to_genes(
        sensitivity._load_json(f"{sensitivity.RESOURCES_DIR}/desirability_constants_extended_sweep_floor.json"),
        battery_genes,
    )
    extended_production = sensitivity._load_json(f"{sensitivity.RESOURCES_DIR}/desirability_constants_extended.json")
    protein_global = extended_production["protein"]["global"]
    dep_gene_codes, dep_scores, dep_gene_categories = sensitivity.load_dependency_raw(
        battery_genes, f"{sensitivity.DATA_DIR}/dependency.parquet"
    )

    def constants_for(thresholds):
        return sensitivity.build_constants_for_thresholds(
            thresholds, rna_floor, extended_floor, protein_global,
            dep_gene_codes, dep_scores, dep_gene_categories,
        )

    default_weights, default_thresholds = sensitivity.default_sample()
    default_rna, default_extended, default_essentiality = constants_for(default_thresholds)
    default_tier_params = sensitivity.tier_params_from_thresholds(default_thresholds)

    ckpt = load_checkpoint()

    # -- baseline (once) --
    if ckpt["baseline"] is None:
        baseline = {}
        for qc in query_caches:
            top10, full_ranked = sensitivity.run_query_under_sample(
                qc, default_weights, default_tier_params, default_rna, default_extended, default_essentiality,
            )
            baseline[qc["name"]] = {"top10": top10, "full_ranked": full_ranked}
        ckpt["baseline"] = baseline
        save_checkpoint(ckpt)
        print(f"[{time.time()-t0:.1f}s] baseline done, checkpointed", flush=True)
    baseline = ckpt["baseline"]

    # -- Dirichlet weight-stability pass --
    dirichlet_samples = sensitivity.dirichlet_weight_samples(N_DIRICHLET, SEED)
    for i in range(len(ckpt["dirichlet"]), len(dirichlet_samples)):
        weights = dirichlet_samples[i]
        per_query = []
        for qc in query_caches:
            top10, full_ranked = sensitivity.run_query_under_sample(
                qc, weights, default_tier_params, default_rna, default_extended, default_essentiality,
            )
            per_query.append(sensitivity._compare_to_baseline(baseline[qc["name"]], top10, full_ranked))
        ckpt["dirichlet"].append({"weights": weights, "per_query": per_query})
        save_checkpoint(ckpt)
        print(f"[{time.time()-t0:.1f}s] dirichlet {i+1}/{len(dirichlet_samples)} done", flush=True)

    # -- joint 12-factor Sobol pass (weights + T1-T6) --
    problem = sensitivity.build_sobol_problem()
    sobol_design = sensitivity.sample_sobol(problem, N_SOBOL)
    for i in range(len(ckpt["sobol"]), len(sobol_design)):
        row = sobol_design[i]
        weights, thresholds = sensitivity.row_to_sample(row, problem)
        tier_params = sensitivity.tier_params_from_thresholds(thresholds)
        rna_c, extended_c, essentiality_c = constants_for(thresholds)
        per_query_rbo = []
        for qc in query_caches:
            top10, full_ranked = sensitivity.run_query_under_sample(
                qc, weights, tier_params, rna_c, extended_c, essentiality_c,
            )
            comparison = sensitivity._compare_to_baseline(baseline[qc["name"]], top10, full_ranked)
            per_query_rbo.append(comparison["rbo_at_10"] if comparison["rbo_at_10"] is not None else 0.0)
        mean_rbo = float(np.mean(per_query_rbo))
        ckpt["sobol"].append({"weights": weights, "thresholds": thresholds, "mean_rbo_at_10": mean_rbo})
        save_checkpoint(ckpt)
        print(f"[{time.time()-t0:.1f}s] sobol {i+1}/{len(sobol_design)} done", flush=True)

    # -- T7/T8 discrete layer-count grid --
    for i in range(len(ckpt["grid"]), len(sensitivity.LAYER_COUNT_GRID)):
        high_min, moderate_min = sensitivity.LAYER_COUNT_GRID[i]
        tier_params = sensitivity.tier_params_from_thresholds(
            default_thresholds, high_min_layers=high_min, moderate_min_layers=moderate_min,
        )
        per_query = []
        for qc in query_caches:
            top10, full_ranked = sensitivity.run_query_under_sample(
                qc, default_weights, tier_params, default_rna, default_extended, default_essentiality,
            )
            per_query.append(sensitivity._compare_to_baseline(baseline[qc["name"]], top10, full_ranked))
        ckpt["grid"].append({
            "high_min_layers": high_min, "moderate_min_layers": moderate_min, "per_query": per_query,
        })
        save_checkpoint(ckpt)
        print(f"[{time.time()-t0:.1f}s] grid {i+1}/{len(sensitivity.LAYER_COUNT_GRID)} done", flush=True)

    print(f"[{time.time()-t0:.1f}s] all samples done -- finalizing", flush=True)

    from SALib.analyze import sobol as sobol_analyze
    sobol_mean_rbo = np.array([r["mean_rbo_at_10"] for r in ckpt["sobol"]])
    Si = sobol_analyze.analyze(problem, sobol_mean_rbo)
    sobol_indices = {"names": problem["names"], "S1": [float(v) for v in Si["S1"]], "ST": [float(v) for v in Si["ST"]]}

    results = {
        "seed": SEED, "n_sobol": N_SOBOL, "n_dirichlet": N_DIRICHLET, "battery_size": len(battery),
        "seconds_per_battery_pass": None,  # not separately timed across resumed invocations
        "dirichlet": ckpt["dirichlet"],
        "sobol": {"design_rows": len(sobol_design), "records": ckpt["sobol"], "indices": sobol_indices},
        "layer_count_grid": ckpt["grid"],
    }
    with open(FINAL_OUT_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[{time.time()-t0:.1f}s] wrote {FINAL_OUT_PATH}", flush=True)
    os.remove(CHECKPOINT_PATH)  # done -- don't leave a stale checkpoint for the next run to misread


if __name__ == "__main__":
    main()
