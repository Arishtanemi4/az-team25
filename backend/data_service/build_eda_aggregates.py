"""One-time build: RNA cross-source concordance and RNA<->protein concordance, written to
`data/processed/eda_aggregates.json`.

Run BY HAND (`python backend/data_service/build_eda_aggregates.py`), never called from an
endpoint or a test -- same convention as `preprocessing/build_gene_role_reference.py`. The
`/data/eda/concordance` endpoint (`eda_service.concordance`) serves this file's content as-is
and returns a **named** unavailable response when it is absent -- it never scans the underlying
79.2M/22.2M/45.9M/4.5M-row tables per request (PRODUCT_SURFACE.md SS3.3).

**Method, ported from `eda/cross_layer/03_rna_cross_source_concordance_eda.ipynb` SS6 and
`04_rna_protein_concordance_eda.ipynb` SS6** -- per-model (not per-gene) Spearman rho: for each
cell line, rank its expression profile across the shared gene set within each source and
correlate the two rank orders. A model is skipped if fewer than `MIN_GENES` genes have a
non-missing value in **both** sources -- too few points for a stable per-model rho, the same
floor those notebooks used. Ported onto `data/processed/` (already ModelID/ensembl_id-resolved)
rather than `data/raw/`, so this script never imports from `eda/`
(`docs/plan/EXECUTE.md`'s standing rule) -- it only reuses the *method*, cited above, not the code.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from data_app.config import DATA_DIR, EDA_AGGREGATES_JSON  # noqa: E402
from data_app.services import table_io  # noqa: E402

MIN_GENES = 30  # genes required in common per model, reused verbatim from the ported notebooks


def _wide(table: str, value_column: str, genes: set) -> pd.DataFrame:
    """One row per ModelID, one column per gene, restricted to `genes` -- read via predicate
    pushdown so only the needed columns' worth of rows are pulled off disk. Some sources (GEO's
    microarray platform in particular) carry more than one probe per (ModelID, ensembl_id) --
    averaged before pivoting, the same collapse-by-mean convention V6-6's
    `joins.average_duplicate_rna_profiles` already applies elsewhere in this pipeline, rather
    than picking an arbitrary probe or letting pivot() raise on the duplicate index."""
    long = table_io.read_gene_filtered(table, genes)
    long = long.groupby(["ModelID", "ensembl_id"], observed=True)[value_column].mean().reset_index()
    wide = long.pivot(index="ModelID", columns="ensembl_id", values=value_column)
    return wide


def _per_model_spearman(a: pd.DataFrame, b: pd.DataFrame, models: list, min_genes: int = MIN_GENES):
    """For each model, correlates its two rank orders across the shared gene columns of `a`
    and `b`. Returns (rho per model, n_genes_used per model) -- one entry per scored model."""
    genes = sorted(set(a.columns) & set(b.columns))
    rhos, n_genes = {}, {}
    for model in models:
        x = a.loc[model, genes]
        y = b.loc[model, genes]
        mask = x.notna() & y.notna()
        if mask.sum() < min_genes:
            continue
        rho, _ = spearmanr(x[mask].values, y[mask].values)
        rhos[model] = rho
        n_genes[model] = int(mask.sum())
    return pd.Series(rhos, name="spearman_rho"), pd.Series(n_genes, name="n_genes_used")


def _summarize(rho: pd.Series, n_genes: pd.Series, models_in_intersection: int) -> dict:
    rho = rho.dropna()
    if rho.empty:
        return {
            "n_models_scored": 0, "n_models_in_intersection": models_in_intersection,
            "median_rho": None, "iqr_low": None, "iqr_high": None,
            "median_n_genes_used": None, "min_genes_floor": MIN_GENES,
        }
    return {
        "n_models_scored": int(len(rho)),
        "n_models_in_intersection": models_in_intersection,
        "median_rho": round(float(rho.median()), 4),
        "iqr_low": round(float(np.percentile(rho, 25)), 4),
        "iqr_high": round(float(np.percentile(rho, 75)), 4),
        "median_n_genes_used": float(n_genes.median()) if not n_genes.empty else None,
        "min_genes_floor": MIN_GENES,
    }


def _shared_genes(table_a: str, table_b: str) -> set:
    genes_a = set(table_io.read_full_column(table_a, "ensembl_id"))
    genes_b = set(table_io.read_full_column(table_b, "ensembl_id"))
    return genes_a & genes_b


def build_rna_cross_source() -> dict:
    print("[build_eda_aggregates] RNA cross-source: resolving shared gene sets...")
    genes_dep_hpa = _shared_genes("expression_rna", "expression_rna_hpa")
    genes_dep_geo = _shared_genes("expression_rna", "expression_rna_geo")
    print(f"  DepMap<->HPA shared genes: {len(genes_dep_hpa):,}")
    print(f"  DepMap<->GEO shared genes: {len(genes_dep_geo):,}")

    print("[build_eda_aggregates] RNA cross-source: reading and pivoting (DepMap vs HPA)...")
    dep_for_hpa = _wide("expression_rna", "log2tpm1", genes_dep_hpa)
    hpa_wide = _wide("expression_rna_hpa", "log2ntpm1", genes_dep_hpa)
    models_dep_hpa = sorted(set(dep_for_hpa.index) & set(hpa_wide.index))
    rho_dep_hpa, n_dep_hpa = _per_model_spearman(dep_for_hpa, hpa_wide, models_dep_hpa)

    print("[build_eda_aggregates] RNA cross-source: reading and pivoting (DepMap vs GEO)...")
    dep_for_geo = _wide("expression_rna", "log2tpm1", genes_dep_geo)
    geo_wide = _wide("expression_rna_geo", "log1p_expr", genes_dep_geo)
    models_dep_geo = sorted(set(dep_for_geo.index) & set(geo_wide.index))
    rho_dep_geo, n_dep_geo = _per_model_spearman(dep_for_geo, geo_wide, models_dep_geo)

    return {
        "depmap_vs_hpa": _summarize(rho_dep_hpa, n_dep_hpa, len(models_dep_hpa)),
        "depmap_vs_geo": _summarize(rho_dep_geo, n_dep_geo, len(models_dep_geo)),
    }


def build_rna_protein() -> dict:
    print("[build_eda_aggregates] RNA<->protein: resolving shared gene set...")
    genes = _shared_genes("expression_rna", "protein")
    print(f"  DepMap RNA<->protein shared genes: {len(genes):,}")

    rna_wide = _wide("expression_rna", "log2tpm1", genes)
    protein_wide = _wide("protein", "zscore", genes)
    models = sorted(set(rna_wide.index) & set(protein_wide.index))
    print(f"  RNA<->protein model intersection: {len(models)}")
    rho, n_genes = _per_model_spearman(rna_wide, protein_wide, models)

    return {"depmap_rna_vs_protein": _summarize(rho, n_genes, len(models))}


def main():
    aggregates = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "rna_cross_source": build_rna_cross_source(),
        "rna_protein": build_rna_protein(),
    }
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(EDA_AGGREGATES_JSON, "w", encoding="utf-8") as fh:
        json.dump(aggregates, fh, indent=2)
    print(f"[build_eda_aggregates] wrote {EDA_AGGREGATES_JSON}")


if __name__ == "__main__":
    main()
