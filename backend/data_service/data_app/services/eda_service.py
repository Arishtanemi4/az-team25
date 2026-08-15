"""Backs GET /data/eda/*."""

import json
import os

from data_app.config import EDA_AGGREGATES_JSON
from data_app.services import table_io
from data_app.services.table_registry import COVERAGE_LAYER_COLUMNS

# Reused from docs/plan/PARAMETERS.md SS9 row 4 (`min_lineage_n`) -- the one group-size floor
# this project has already cited (RNA lineage-calibration fallback), rather than inventing a
# second, unrelated "small group" threshold just for this chart.
MIN_LINEAGE_N = 15

_LAYER_COLUMNS = COVERAGE_LAYER_COLUMNS


def coverage() -> dict:
    """Served as-is from build_manifest.json -- no table read, matching the "free" cost row in
    PRODUCT_SURFACE.md SS3.3."""
    manifest = _read_manifest()
    return {"coverage_state_counts": manifest.get("coverage_state_counts", {})}


def _read_manifest() -> dict:
    from data_app.config import BUILD_MANIFEST_JSON
    if not os.path.exists(BUILD_MANIFEST_JSON):
        return {}
    with open(BUILD_MANIFEST_JSON, encoding="utf-8") as fh:
        return json.load(fh)


def lineage_coverage() -> dict:
    coverage_df = table_io.read_full("coverage")
    cell_lines = table_io.read_full("cell_lines")[["ModelID", "lineage"]]
    merged = coverage_df.merge(cell_lines, on="ModelID", how="left").dropna(subset=["lineage"])

    rows = []
    for lineage, group in merged.groupby("lineage", observed=True):
        n = len(group)
        for column in _LAYER_COLUMNS:
            if column.endswith("_available"):
                n_measured = int((group[column] == True).sum())  # noqa: E712 -- explicit bool compare, not truthiness
            else:
                n_measured = int((group[column] != "not_assayed").sum())
            rows.append({
                "lineage": lineage,
                "layer": column,
                "n": n,
                "n_measured": n_measured,
                "measured_fraction": round(n_measured / n, 4) if n else 0.0,
                "flagged": n < MIN_LINEAGE_N,
            })

    return {"min_lineage_n": MIN_LINEAGE_N, "rows": rows}


def concordance() -> dict:
    if not os.path.exists(EDA_AGGREGATES_JSON):
        return {
            "available": False,
            "reason": "data/processed/eda_aggregates.json has not been built on this machine",
            "build_command": "python backend/data_service/build_eda_aggregates.py",
        }
    with open(EDA_AGGREGATES_JSON, encoding="utf-8") as fh:
        aggregates = json.load(fh)
    return {"available": True, **aggregates}
