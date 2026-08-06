import json
from datetime import datetime, timezone

import pandas as pd

DATA_DIR = "data/processed"
RESOURCES_DIR = "scoring/resources"

STRONG_DEPENDENCY_CUTOFF = -1.0          # matches eda/single_layer/15's own "strong dependency"
                                          # bucket (5.12% of all gene x model cells fall below it)
PAN_ESSENTIAL_FRACTION_THRESHOLD = 0.90
MIN_N = 30  # same >=30-line calibration gate the other continuous layers use -- below this a
            # fraction estimated from a handful of lines is noise, not a real proxy.


def load_dependency(data_dir: str = DATA_DIR) -> pd.DataFrame:
    return pd.read_csv(
        f"{data_dir}/dependency.csv",
        dtype={"ModelID": "category", "ensembl_id": "category", "dependency_score": "float32"},
    )


def compute_fraction_strong(dependency: pd.DataFrame, cutoff: float = STRONG_DEPENDENCY_CUTOFF,
                            min_n: int = MIN_N) -> pd.DataFrame:
    dependency = dependency.copy()
    dependency["strong"] = dependency["dependency_score"] < cutoff
    grouped = dependency.groupby("ensembl_id", observed=True)["strong"].agg(["sum", "count"])
    grouped = grouped[grouped["count"] >= min_n]
    grouped["fraction_strong"] = grouped["sum"] / grouped["count"]
    return grouped


def flag_pan_essential(grouped: pd.DataFrame,
                       threshold: float = PAN_ESSENTIAL_FRACTION_THRESHOLD) -> pd.DataFrame:
    return grouped[grouped["fraction_strong"] >= threshold]


def sanity_check_top_hits(pan_essential: pd.DataFrame, data_dir: str = DATA_DIR, n: int = 15) -> pd.DataFrame:
    gene_reference = pd.read_csv(f"{data_dir}/gene_reference.csv", usecols=["ensembl_id", "symbol"])
    top_hits = pan_essential.sort_values("fraction_strong", ascending=False).head(n)
    return top_hits.merge(gene_reference, left_index=True, right_on="ensembl_id")


def build_constants(grouped: pd.DataFrame, pan_essential: pd.DataFrame) -> dict:
    pan_essential_dict = {
        gene: round(float(row["fraction_strong"]), 4) for gene, row in pan_essential.iterrows()
    }
    return {
        "version": "1.0.0",
        "generated": datetime.now(timezone.utc).date().isoformat(),
        "source_script": "scoring/build_essentiality_constants.py",
        "source_data": "data/processed/dependency.csv",
        "method": {
            "description": (
                "Per gene, the fraction of measured lines (N>=30) with a Chronos dependency_score "
                "below strong_dependency_cutoff. Flagged pan_essential if that fraction is >= "
                "pan_essential_fraction_threshold. strong_dependency_cutoff=-1.0 matches "
                "eda/single_layer/15's own 'strong dependency' bucket (5.12% of all gene x model "
                "cells fall below it). pan_essential_fraction_threshold=0.90 and the cutoff itself "
                "are EDA-inherited, not yet cited to a specific published numeric threshold -- "
                "Hart et al. (2014/2017) and Dempster et al. (2021) are the closest available "
                "references (core-essential gold standard; the Chronos method itself), not the "
                "source of these exact numbers. Covered by the V6-3 sensitivity sweep "
                "(docs/plan/PARAMETERS.md row 40)."
            ),
            "strong_dependency_cutoff": STRONG_DEPENDENCY_CUTOFF,
            "pan_essential_fraction_threshold": PAN_ESSENTIAL_FRACTION_THRESHOLD,
            "min_calibration_n": MIN_N,
        },
        "n_genes_total": int(len(grouped)),
        "n_genes_flagged": int(len(pan_essential)),
        "pan_essential": pan_essential_dict,
    }


def main(data_dir: str = DATA_DIR, resources_dir: str = RESOURCES_DIR) -> None:
    dependency = load_dependency(data_dir)
    print(f"loaded {len(dependency):,} dependency rows")

    grouped = compute_fraction_strong(dependency)
    print(f"n_genes_total (>= {MIN_N} measured lines) = {len(grouped):,}")

    pan_essential = flag_pan_essential(grouped)
    print(f"n_genes_flagged pan_essential = {len(pan_essential):,}")

    top_hits = sanity_check_top_hits(pan_essential, data_dir)
    print(top_hits[["ensembl_id", "symbol", "fraction_strong"]].to_string(index=False))

    out_path = f"{resources_dir}/common_essential_genes.json"
    with open(out_path, "w") as f:
        json.dump(build_constants(grouped, pan_essential), f, indent=2)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
