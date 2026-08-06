import json
from datetime import datetime, timezone

import pandas as pd

DATA_DIR = "data/processed"
RESOURCES_DIR = "scoring/resources"

MIN_LINEAGE_N = 15  
MIN_LINEAGE_N_FLOOR = 5  

LOW_PERCENTILE = 0.1   # L, the floor
HIGH_PERCENTILE = 0.9  # T, the target


def load_expression_with_lineage(data_dir: str = DATA_DIR) -> pd.DataFrame:
    expr = pd.read_csv(
        f"{data_dir}/expression_rna.csv",
        dtype={"ModelID": "category", "ensembl_id": "category", "log2tpm1": "float32"},
    )
    cell_lines = pd.read_csv(f"{data_dir}/cell_lines.csv", usecols=["ModelID", "lineage"])
    lineage_map = cell_lines.set_index("ModelID")["lineage"]
    expr["lineage"] = expr["ModelID"].map(lineage_map).astype("category")
    return expr


def compute_global_lt(expr: pd.DataFrame) -> pd.DataFrame:
    global_q = expr.groupby("ensembl_id", observed=True)["log2tpm1"].quantile(
        [LOW_PERCENTILE, HIGH_PERCENTILE]
    ).unstack()
    global_q.columns = ["L", "T"]
    return global_q


def compute_lineage_lt(expr: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    grouped = expr.groupby(["lineage", "ensembl_id"], observed=True)["log2tpm1"]
    lineage_counts = grouped.count()
    lineage_q = grouped.quantile([LOW_PERCENTILE, HIGH_PERCENTILE]).unstack()
    lineage_q.columns = ["L", "T"]
    return lineage_q, lineage_counts


def apply_discrimination_gate(quantiles: pd.DataFrame) -> pd.DataFrame:
    return quantiles[quantiles["T"] > quantiles["L"]]


def to_global_dict(global_kept: pd.DataFrame) -> dict:
    return {
        gene: {"L": round(float(row["L"]), 4), "T": round(float(row["T"]), 4)}
        for gene, row in global_kept.iterrows()
    }


def to_per_lineage_dict(lineage_kept: pd.DataFrame, lineage_counts: pd.Series,
                        min_n: int) -> tuple[dict, dict]:
    per_lineage = {}
    counts = {}
    # one pass over the kept pairs, nesting by lineage and dropping anything under the floor
    for (lineage, gene), row in lineage_kept.iterrows():
        n = int(lineage_counts.get((lineage, gene), 0))
        if n < min_n:
            continue
        per_lineage.setdefault(lineage, {})[gene] = {
            "L": round(float(row["L"]), 4), "T": round(float(row["T"]), 4),
        }
        counts.setdefault(lineage, {})[gene] = n
    return per_lineage, counts


def build_production_constants(global_q: pd.DataFrame, global_kept: pd.DataFrame,
                               per_lineage_dict: dict) -> dict:
    return {
        "version": "1.0.0",
        "generated": datetime.now(timezone.utc).date().isoformat(),
        "source_script": "scoring/build_desirability_constants.py",
        "source_data": "data/processed/expression_rna.csv, data/processed/cell_lines.csv",
        "method": (
            "Global: L = 10th percentile, T = 90th percentile of log2(TPM+1) pooled across every "
            "measured line, per gene. Per-lineage override: the same 10th/90th percentile computed "
            "within cell_lines.lineage, for any (gene, lineage) pair with >=15 measured lines "
            "(CONSTRAINTS.md SS3's lineage-size threshold), falling back to the global value below "
            "that count. P5 redesign: computed in-sample from DepMap's own expression_rna.csv, not "
            "from Jin et al. (2023) as in the original v4 design -- that file is reserved for "
            "external validation and is never read for calibration."
        ),
        "min_lineage_n": MIN_LINEAGE_N,
        "gene_filter": "kept genes where T > L (any discrimination)",
        "n_genes_total": len(global_q),
        "n_genes_kept": len(global_kept),
        "n_lineages_with_overrides": len(per_lineage_dict),
        "global": to_global_dict(global_kept),
        "per_lineage": per_lineage_dict,
    }


def build_sweep_floor_constants(global_q: pd.DataFrame, global_kept: pd.DataFrame,
                                per_lineage_floor: dict, lineage_gene_counts: dict) -> dict:
    return {
        "version": "1.0.0",
        "generated": datetime.now(timezone.utc).date().isoformat(),
        "source_script": "scoring/build_desirability_constants.py",
        "source_data": "data/processed/expression_rna.csv, data/processed/cell_lines.csv",
        "method": (
            "Same 10th/90th percentile procedure as the production file, re-cut at a lower "
            "MIN_LINEAGE_N_FLOOR so V6-3's sensitivity sweep can honestly sweep min_lineage_n below "
            "the production threshold of 15. Adds lineage_gene_counts (absent from the production "
            "file) so the sweep can re-gate per sampled threshold via a dict filter, not a "
            "recomputation. docs/plan/PARAMETERS.md SS11, docs/plan/SENSITIVITY.md SS3-4."
        ),
        "min_lineage_n_floor": MIN_LINEAGE_N_FLOOR,
        "gene_filter": "kept genes where T > L (any discrimination)",
        "n_genes_total": len(global_q),
        "n_genes_kept": len(global_kept),
        "n_lineages_with_overrides": len(per_lineage_floor),
        "global": to_global_dict(global_kept),
        "per_lineage": per_lineage_floor,
        "lineage_gene_counts": lineage_gene_counts,
    }


def main(data_dir: str = DATA_DIR, resources_dir: str = RESOURCES_DIR) -> None:
    expr = load_expression_with_lineage(data_dir)
    print(f"{len(expr):,} RNA rows loaded")

    global_q = compute_global_lt(expr)
    lineage_q, lineage_counts = compute_lineage_lt(expr)
    print(f"{len(global_q):,} genes have a global RNA measurement")

    global_kept = apply_discrimination_gate(global_q)
    lineage_kept = apply_discrimination_gate(lineage_q)
    print(f"Global: kept {len(global_kept):,} of {len(global_q):,} genes")
    print(f"Per-lineage: {len(lineage_kept):,} of {len(lineage_q):,} pairs pass the discrimination gate")

    per_lineage_dict, _ = to_per_lineage_dict(lineage_kept, lineage_counts, MIN_LINEAGE_N)
    per_lineage_floor, lineage_gene_counts = to_per_lineage_dict(
        lineage_kept, lineage_counts, MIN_LINEAGE_N_FLOOR
    )

    production_path = f"{resources_dir}/desirability_constants.json"
    with open(production_path, "w") as f:
        json.dump(build_production_constants(global_q, global_kept, per_lineage_dict), f, indent=2)
    print(f"Wrote {production_path}")

    floor_path = f"{resources_dir}/desirability_constants_rna_sweep_floor.json"
    with open(floor_path, "w") as f:
        json.dump(
            build_sweep_floor_constants(global_q, global_kept, per_lineage_floor, lineage_gene_counts),
            f, indent=2,
        )
    print(f"Wrote {floor_path}")

    n_production = sum(len(v) for v in per_lineage_dict.values())
    n_floor = sum(len(v) for v in per_lineage_floor.values())
    print(f"{n_floor:,} (lineage, gene) overrides at floor={MIN_LINEAGE_N_FLOOR} "
          f"vs. {n_production:,} at production MIN_LINEAGE_N={MIN_LINEAGE_N}")


if __name__ == "__main__":
    main()
