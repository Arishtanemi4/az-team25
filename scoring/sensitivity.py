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
