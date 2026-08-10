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
