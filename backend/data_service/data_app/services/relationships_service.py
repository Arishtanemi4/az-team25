"""Backs GET /data/relationships/*. Every read here is gene-keyed via Parquet predicate
pushdown (table_io.read_gene_filtered) -- a one-gene query touches ~2,000 rows, never a full
table scan (PRODUCT_SURFACE.md SS3.2)."""

import numpy as np
import pandas as pd
from fastapi import HTTPException

from data_app.services import table_io

# layer -> (processed table name, its value column, a human-readable unit label). The six
# per-(ModelID, ensembl_id) measurement tables this project has -- mutations/fusions are
# categorical evidence, not a single continuous value, so they are not offered as a "layer" here.
LAYER_TABLES = {
    "rna": {"table": "expression_rna", "value_column": "log2tpm1", "unit": "log2(TPM+1)"},
    "rna_hpa": {"table": "expression_rna_hpa", "value_column": "log2ntpm1", "unit": "log2(nTPM+1)"},
    "rna_geo": {"table": "expression_rna_geo", "value_column": "log1p_expr", "unit": "log1p(expr)"},
    "protein": {"table": "protein", "value_column": "zscore", "unit": "z-score"},
    "dependency": {"table": "dependency", "value_column": "dependency_score", "unit": "Chronos dependency score"},
    "copy_number": {"table": "copy_number", "value_column": "copy_number", "unit": "relative copy number"},
}


def _resolve_layer(layer: str) -> dict:
    if layer not in LAYER_TABLES:
        raise HTTPException(
            status_code=422,
            detail=f"unknown layer {layer!r}; choose one of {sorted(LAYER_TABLES)}",
        )
    return LAYER_TABLES[layer]


def _resolve_symbol(ensembl_id: str) -> str | None:
    matches = table_io.read_gene_filtered("gene_reference", [ensembl_id])
    if matches.empty:
        return None
    return str(matches.iloc[0]["symbol"])


def gene_cell_line(ensembl_id: str, layer: str, limit: int) -> dict:
    spec = _resolve_layer(layer)
    values = table_io.read_gene_filtered(spec["table"], [ensembl_id])
    values = values[values["ensembl_id"] == ensembl_id]

    cell_lines = table_io.read_full("cell_lines")[["ModelID", "cell_line_name", "lineage"]]
    merged = values.merge(cell_lines, on="ModelID", how="left")

    points = [
        {
            "model_id": row["ModelID"],
            "cell_line_name": row["cell_line_name"] if not pd.isna(row["cell_line_name"]) else None,
            "lineage": row["lineage"] if not pd.isna(row["lineage"]) else None,
            "value": float(row[spec["value_column"]]),
        }
        for _, row in merged.head(limit).iterrows()
        if not pd.isna(row[spec["value_column"]])
    ]

    return {
        "ensembl_id": ensembl_id,
        "symbol": _resolve_symbol(ensembl_id),
        "layer": layer,
        "unit": spec["unit"],
        "n_measured": len(merged),
        "n_models_total": len(table_io.read_full("cell_lines")),
        "points": points,
    }


def cell_line_disease(lineage: str | None) -> dict:
    cell_lines = table_io.read_full("cell_lines")
    n_models_total = len(cell_lines)

    scoped = cell_lines if lineage is None else cell_lines[cell_lines["lineage"] == lineage]
    grouped = (
        scoped.dropna(subset=["lineage", "primary_disease"])
        .groupby(["lineage", "primary_disease"], observed=True)
        .size()
        .reset_index(name="n_models")
        .sort_values(["lineage", "primary_disease"])
    )
    cells = [
        {"lineage": row["lineage"], "primary_disease": row["primary_disease"], "n_models": int(row["n_models"])}
        for _, row in grouped.iterrows()
    ]
    return {"n_models_total": n_models_total, "cells": cells}


def gene_disease(ensembl_id: str, layer: str, min_n: int) -> dict:
    spec = _resolve_layer(layer)
    values = table_io.read_gene_filtered(spec["table"], [ensembl_id])
    values = values[values["ensembl_id"] == ensembl_id]

    cell_lines = table_io.read_full("cell_lines")[["ModelID", "primary_disease"]].dropna(subset=["primary_disease"])
    n_models_by_disease = cell_lines.groupby("primary_disease", observed=True).size()

    measured = values.merge(cell_lines, on="ModelID", how="inner")

    groups = []
    for disease, n_models in n_models_by_disease.items():
        disease_values = measured.loc[measured["primary_disease"] == disease, spec["value_column"]].dropna()
        n_measured = len(disease_values)
        if n_measured == 0:
            median = iqr_low = iqr_high = None
        else:
            median = float(np.median(disease_values))
            iqr_low = float(np.percentile(disease_values, 25))
            iqr_high = float(np.percentile(disease_values, 75))
        groups.append({
            "primary_disease": disease,
            "n_measured": int(n_measured),
            "n_models": int(n_models),
            "median": median,
            "iqr_low": iqr_low,
            "iqr_high": iqr_high,
            "flagged": n_measured < min_n,
        })
    groups.sort(key=lambda g: g["n_measured"], reverse=True)

    return {
        "ensembl_id": ensembl_id,
        "symbol": _resolve_symbol(ensembl_id),
        "layer": layer,
        "min_n": min_n,
        "groups": groups,
    }
