import json
import time

import numpy as np
import pandas as pd
from scipy.stats import kendalltau
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze

import correlation
import evidence
import score

DATA_DIR = "data/processed"
RESOURCES_DIR = "scoring/resources"

LAYER_ORDER = ["rna", "protein", "fusion", "copy_number", "mutation", "dependency"]

LAYER_FRAME_NAMES = [
    "expression_rna", "expression_rna_hpa", "expression_rna_geo", "protein",
    "dependency", "copy_number", "mutations", "fusions",
]

THRESHOLD_BOUNDS = {
    "min_lineage_n": [5, 50],                        # T1, default 15
    "min_calibration_n": [10, 100],                  # T2, default 30 (shared with essentiality's MIN_N)
    "tier_insufficient_fraction": [0.2, 0.8],         # T3, default 0.5
    "tier_fraction_threshold": [0.5, 0.95],           # T4, default 0.8 (shared High/Moderate literal)
    "strong_dependency_cutoff": [-1.5, -0.5],         # T5, default -1.0
    "pan_essential_fraction_threshold": [0.75, 0.99], # T6, default 0.90
}

LAYER_COUNT_GRID = [(2, 1), (3, 1), (3, 2), (4, 2), (4, 3)]  # (high_min_layers, moderate_min_layers)

RBO_P = 0.9 


def rbo_at_10(ranked_a, ranked_b, p=RBO_P):
    k = min(10, len(ranked_a), len(ranked_b))
    if k == 0:
        return None
    seen_a, seen_b, overlaps = set(), set(), []
    for d in range(1, k + 1):
        seen_a.add(ranked_a[d - 1])
        seen_b.add(ranked_b[d - 1])
        overlaps.append(len(seen_a & seen_b))
    weighted_sum = sum((overlaps[d - 1] / d) * (p ** d) for d in range(1, k + 1))
    extrapolation_term = (overlaps[-1] / k) * (p ** k)
    return extrapolation_term + ((1 - p) / p) * weighted_sum


def kendall_tau_full(ranked_a, ranked_b):
    common = set(ranked_a) & set(ranked_b)
    if len(common) < 2:
        return None
    rank_a = {m: i for i, m in enumerate(ranked_a)}
    rank_b = {m: i for i, m in enumerate(ranked_b)}
    xs = [rank_a[m] for m in ranked_a if m in common]
    ys = [rank_b[m] for m in ranked_a if m in common]
    tau, _p_value = kendalltau(xs, ys)
    return float(tau)


def top10_membership_change(ranked_a, ranked_b):
    set_a, set_b = set(ranked_a[:10]), set(ranked_b[:10])
    return set_a != set_b, len(set_a ^ set_b)


def load_dependency_raw(battery_genes, path=f"{DATA_DIR}/dependency.parquet"):
    df = pd.read_parquet(
        path, engine="pyarrow", columns=["ModelID", "ensembl_id", "dependency_score"],
        filters=[("ensembl_id", "in", list(battery_genes))],
    )
    gene_codes_series = df["ensembl_id"].astype("category")
    gene_codes = gene_codes_series.cat.codes.to_numpy()
    gene_categories = gene_codes_series.cat.categories.to_numpy()
    scores = df["dependency_score"].to_numpy(dtype="float32")
    return gene_codes, scores, gene_categories


def essentiality_for_sample(gene_codes, scores, cutoff, fraction_threshold, min_n, gene_categories):
    n_genes = len(gene_categories)
    strong = (scores < cutoff).astype(np.float64)
    n_strong = np.bincount(gene_codes, weights=strong, minlength=n_genes)
    n_total = np.bincount(gene_codes, minlength=n_genes).astype(np.float64)
    frac = np.divide(n_strong, n_total, out=np.zeros(n_genes), where=n_total > 0)
    flagged = np.nonzero((n_total >= min_n) & (frac >= fraction_threshold))[0]
    return {"pan_essential": {str(gene_categories[i]): float(frac[i]) for i in flagged}}


def _load_json(path):
    with open(path) as f:
        return json.load(f)


def restrict_rna_floor_to_genes(rna_floor, genes):
    genes = set(genes)
    per_lineage = {
        lineage: {g: lt for g, lt in gene_map.items() if g in genes}
        for lineage, gene_map in rna_floor["per_lineage"].items()
    }
    per_lineage = {lineage: gm for lineage, gm in per_lineage.items() if gm}
    lineage_gene_counts = {
        lineage: {g: n for g, n in count_map.items() if g in genes}
        for lineage, count_map in rna_floor["lineage_gene_counts"].items()
    }
    return {
        "global": {g: lt for g, lt in rna_floor["global"].items() if g in genes},
        "per_lineage": per_lineage,
        "lineage_gene_counts": lineage_gene_counts,
    }


def restrict_extended_floor_to_genes(extended_floor, genes):
    genes = set(genes)
    return {
        "dependency": {"global": {g: lt for g, lt in extended_floor["dependency"]["global"].items() if g in genes}},
        "copy_number": {"global": {g: lt for g, lt in extended_floor["copy_number"]["global"].items() if g in genes}},
        "dependency_gene_counts": {g: n for g, n in extended_floor["dependency_gene_counts"].items() if g in genes},
        "copy_number_gene_counts": {g: n for g, n in extended_floor["copy_number_gene_counts"].items() if g in genes},
    }


def build_rna_constants_at_threshold(rna_floor, min_lineage_n):
    per_lineage = {}
    counts = rna_floor["lineage_gene_counts"]
    for lineage, genes in rna_floor["per_lineage"].items():
        lineage_counts = counts.get(lineage, {})
        kept = {g: lt for g, lt in genes.items() if lineage_counts.get(g, 0) >= min_lineage_n}
        if kept:
            per_lineage[lineage] = kept
    return {"global": rna_floor["global"], "per_lineage": per_lineage}


def build_extended_constants_at_threshold(extended_floor, min_calibration_n, protein_global):
    dep_counts = extended_floor["dependency_gene_counts"]
    dep_kept = {
        g: lt for g, lt in extended_floor["dependency"]["global"].items()
        if dep_counts.get(g, 0) >= min_calibration_n
    }
    cn_counts = extended_floor["copy_number_gene_counts"]
    cn_kept = {
        g: lt for g, lt in extended_floor["copy_number"]["global"].items()
        if cn_counts.get(g, 0) >= min_calibration_n
    }
    return {
        "protein": {"global": protein_global},
        "dependency": {"global": dep_kept},
        "copy_number": {"global": cn_kept},
    }


def build_constants_for_thresholds(thresholds, rna_floor, extended_floor, protein_global,
                                    dep_gene_codes, dep_scores, dep_gene_categories):
    rna_constants = build_rna_constants_at_threshold(rna_floor, thresholds["min_lineage_n"])
    extended_constants = build_extended_constants_at_threshold(
        extended_floor, thresholds["min_calibration_n"], protein_global,
    )
    essentiality_constants = essentiality_for_sample(
        dep_gene_codes, dep_scores, thresholds["strong_dependency_cutoff"],
        thresholds["pan_essential_fraction_threshold"], thresholds["min_calibration_n"],
        dep_gene_categories,
    )
    return rna_constants, extended_constants, essentiality_constants


def dirichlet_weight_samples(n_samples, seed, alpha=None):
    rng = np.random.default_rng(seed)
    if alpha is None:
        alpha = np.ones(len(LAYER_ORDER))
    draws = rng.dirichlet(alpha, size=n_samples)
    return [dict(zip(LAYER_ORDER, row)) for row in draws]


def build_sobol_problem():
    names = [f"weight_{layer}" for layer in LAYER_ORDER] + list(THRESHOLD_BOUNDS.keys())
    bounds = [[0.0, 1.0]] * len(LAYER_ORDER) + list(THRESHOLD_BOUNDS.values())
    return {"num_vars": len(names), "names": names, "bounds": bounds}


def sample_sobol(problem, n):
    return sobol_sample.sample(problem, n)


def row_to_sample(row, problem):
    raw = dict(zip(problem["names"], row))
    raw_weights = np.array([raw[f"weight_{layer}"] for layer in LAYER_ORDER])
    weight_sum = raw_weights.sum()
    if weight_sum < 1e-9:
        normalized = np.ones(len(LAYER_ORDER)) / len(LAYER_ORDER)
    else:
        normalized = raw_weights / weight_sum
    layer_weights = dict(zip(LAYER_ORDER, normalized))
    thresholds = {name: raw[name] for name in THRESHOLD_BOUNDS}
    thresholds["min_lineage_n"] = int(round(thresholds["min_lineage_n"]))
    thresholds["min_calibration_n"] = int(round(thresholds["min_calibration_n"]))
    return layer_weights, thresholds


def default_sample():
    layer_weights = dict(evidence.LAYER_WEIGHTS)
    thresholds = {
        "min_lineage_n": 15, "min_calibration_n": 30,
        "tier_insufficient_fraction": 0.5, "tier_fraction_threshold": 0.8,
        "strong_dependency_cutoff": -1.0, "pan_essential_fraction_threshold": 0.90,
    }
    return layer_weights, thresholds
