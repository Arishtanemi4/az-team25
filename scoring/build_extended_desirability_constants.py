import json
from datetime import datetime, timezone

import pandas as pd

from build_desirability_constants import apply_discrimination_gate, to_global_dict

DATA_DIR = "data/processed"
RESOURCES_DIR = "scoring/resources"

MIN_N = 30       # production gate -- same >=30-line calibration threshold every continuous
                 # layer here uses; below this a percentile is noise, not a real estimate.
MIN_N_FLOOR = 10 

LOW_PERCENTILE = 0.1
HIGH_PERCENTILE = 0.9


def compute_protein_lt(data_dir: str = DATA_DIR) -> tuple[float, float, int]:
    protein = pd.read_csv(
        f"{data_dir}/protein.csv",
        dtype={"ModelID": "category", "ensembl_id": "category", "zscore": "float32", "detected": "bool"},
    )
    detected = protein.loc[protein["detected"], "zscore"]
    L, T = float(detected.quantile(LOW_PERCENTILE)), float(detected.quantile(HIGH_PERCENTILE))
    return L, T, len(detected)


def compute_gene_lt(df: pd.DataFrame, value_col: str) -> tuple[pd.DataFrame, pd.Series]:
    grouped = df.groupby("ensembl_id", observed=True)[value_col]
    counts = grouped.count()
    quantiles = grouped.quantile([LOW_PERCENTILE, HIGH_PERCENTILE]).unstack()
    quantiles.columns = ["L", "T"]
    return quantiles, counts


def load_dependency_strength(data_dir: str = DATA_DIR) -> pd.DataFrame:
    dependency = pd.read_csv(
        f"{data_dir}/dependency.csv",
        dtype={"ModelID": "category", "ensembl_id": "category", "dependency_score": "float32"},
    )
    dependency["strength"] = -dependency["dependency_score"]
    return dependency


def load_copy_number(data_dir: str = DATA_DIR) -> pd.DataFrame:
    return pd.read_csv(
        f"{data_dir}/copy_number.csv",
        usecols=["ModelID", "ensembl_id", "copy_number"],
        dtype={"ModelID": "category", "ensembl_id": "category", "copy_number": "float32"},
    )


def apply_count_floor(quantiles: pd.DataFrame, counts: pd.Series, min_n: int) -> pd.DataFrame:
    eligible = counts[counts >= min_n].index
    return quantiles.loc[quantiles.index.isin(eligible)]


def to_gene_count_dict(kept: pd.DataFrame, counts: pd.Series) -> dict:
    return {gene: int(counts.get(gene, 0)) for gene in kept.index}


def build_production_constants(protein_L: float, protein_T: float, n_detected: int,
                               dep_q: pd.DataFrame, dep_kept: pd.DataFrame,
                               cn_q: pd.DataFrame, cn_kept: pd.DataFrame) -> dict:
    return {
        "version": "1.0.0",
        "generated": datetime.now(timezone.utc).date().isoformat(),
        "source_script": "scoring/build_extended_desirability_constants.py",
        "source_data": "data/processed/protein.csv, data/processed/dependency.csv, data/processed/copy_number.csv",
        "method": (
            "Protein: L = 10th percentile, T = 90th percentile of zscore pooled across all "
            "detected=True rows (one global pair, not per-gene). Dependency: L = 10th percentile, "
            "T = 90th percentile of -dependency_score ('dependency strength'), per gene, gated at "
            "N>=30 lines/gene. Copy number: L = 10th percentile, T = 90th percentile of the "
            "continuous copy_number value, per gene, same N>=30 gate. r_default = 1 (Derringer & "
            "Suich linear desirability, matching the RNA calibration). P5 note: protein's zscore "
            "is now detection-floor imputed upstream in preprocessing/impute.py -- this "
            "calibration still filters on detected==True exactly as it always did, so the "
            "imputed floor value is correctly excluded from the percentile."
        ),
        "min_calibration_n": MIN_N,
        "gene_filter": "kept genes where T > L (any discrimination), same rule as the RNA calibration",
        "no_per_lineage_override": (
            "Unlike RNA, none of these three layers has a per-lineage override -- RNA's own "
            "in-sample per-lineage override is a deliberate one-layer investment, not yet "
            "extended here. Documented simplification, not a silent gap."
        ),
        "protein": {
            "n_detected_measurements": n_detected,
            "global": {"L": round(protein_L, 4), "T": round(protein_T, 4)},
        },
        "dependency": {
            "n_genes_total": len(dep_q), "n_genes_kept": len(dep_kept),
            "global": to_global_dict(dep_kept),
        },
        "copy_number": {
            "n_genes_total": len(cn_q), "n_genes_kept": len(cn_kept),
            "global": to_global_dict(cn_kept),
        },
    }


def build_sweep_floor_constants(dep_q: pd.DataFrame, dep_kept_floor: pd.DataFrame, dep_counts: pd.Series,
                                cn_q: pd.DataFrame, cn_kept_floor: pd.DataFrame, cn_counts: pd.Series) -> dict:
    return {
        "version": "1.0.0",
        "generated": datetime.now(timezone.utc).date().isoformat(),
        "source_script": "scoring/build_extended_desirability_constants.py",
        "source_data": "data/processed/dependency.csv, data/processed/copy_number.csv",
        "method": (
            "Same per-gene 10th/90th percentile procedure as the production dependency/copy_number "
            "calibration, re-cut at a lower MIN_N_FLOOR so V6-3's sensitivity sweep can honestly "
            "sweep min_calibration_n below the production threshold of 30. Adds "
            "dependency_gene_counts/copy_number_gene_counts (absent from the production file) so the "
            "sweep can re-gate per sampled threshold via a dict filter, not a recomputation. Protein is "
            "not included -- it has no per-gene N-gate (one pooled global pair). "
            "docs/plan/PARAMETERS.md SS11, docs/plan/SENSITIVITY.md SS3-4."
        ),
        "min_calibration_n_floor": MIN_N_FLOOR,
        "gene_filter": "kept genes where T > L (any discrimination)",
        "dependency": {
            "n_genes_total": len(dep_q), "n_genes_kept": len(dep_kept_floor),
            "global": to_global_dict(dep_kept_floor),
        },
        "copy_number": {
            "n_genes_total": len(cn_q), "n_genes_kept": len(cn_kept_floor),
            "global": to_global_dict(cn_kept_floor),
        },
        "dependency_gene_counts": to_gene_count_dict(dep_kept_floor, dep_counts),
        "copy_number_gene_counts": to_gene_count_dict(cn_kept_floor, cn_counts),
    }


def main(data_dir: str = DATA_DIR, resources_dir: str = RESOURCES_DIR) -> None:
    protein_L, protein_T, n_detected = compute_protein_lt(data_dir)
    print(f"protein: n_detected={n_detected:,}, L={protein_L:.4f}, T={protein_T:.4f}")

    dependency = load_dependency_strength(data_dir)
    dep_q, dep_counts = compute_gene_lt(dependency, "strength")
    del dependency

    copy_number = load_copy_number(data_dir)
    cn_q, cn_counts = compute_gene_lt(copy_number, "copy_number")
    del copy_number

    dep_production = apply_count_floor(apply_discrimination_gate(dep_q), dep_counts, MIN_N)
    cn_production = apply_count_floor(apply_discrimination_gate(cn_q), cn_counts, MIN_N)
    print(f"dependency: kept {len(dep_production):,} of {len(dep_q):,} genes at MIN_N={MIN_N}")
    print(f"copy_number: kept {len(cn_production):,} of {len(cn_q):,} genes at MIN_N={MIN_N}")

    production_path = f"{resources_dir}/desirability_constants_extended.json"
    with open(production_path, "w") as f:
        json.dump(
            build_production_constants(protein_L, protein_T, n_detected, dep_q, dep_production, cn_q, cn_production),
            f, indent=2,
        )
    print(f"Wrote {production_path}")

    dep_floor = apply_count_floor(apply_discrimination_gate(dep_q), dep_counts, MIN_N_FLOOR)
    cn_floor = apply_count_floor(apply_discrimination_gate(cn_q), cn_counts, MIN_N_FLOOR)
    print(f"dependency: {len(dep_floor):,} genes at floor={MIN_N_FLOOR} vs. {len(dep_production):,} at production MIN_N={MIN_N}")
    print(f"copy_number: {len(cn_floor):,} genes at floor={MIN_N_FLOOR} vs. {len(cn_production):,} at production MIN_N={MIN_N}")

    floor_path = f"{resources_dir}/desirability_constants_extended_sweep_floor.json"
    with open(floor_path, "w") as f:
        json.dump(
            build_sweep_floor_constants(dep_q, dep_floor, dep_counts, cn_q, cn_floor, cn_counts),
            f, indent=2,
        )
    print(f"Wrote {floor_path}")


if __name__ == "__main__":
    main()
