"""Backs GET /data/validation/*. Reads `validation/*.csv` as-is and reproduces
`validation/validation_report.ipynb`'s verdicts -- never recomputes them (PRODUCT_SURFACE.md
SS3.4)."""

import pandas as pd
from fastapi import HTTPException

from data_app.config import VALIDATION_DIR

STUDY_COMPARABILITY_CSV = "study_comparability.csv"

# Allow-listed per-study files -- never used to build a path from user input.
STUDY_FILES = {
    "jin2023": "jin2023_comparison.csv",
    "tumorcomparer": "tumorcomparer_comparison.csv",
    "celligner": "celligner_coverage.csv",
    "netcellmatch": "netcellmatch_coverage.csv",
    "tumorcomparer_tier2": "tumorcomparer_tier2_status.csv",
}

_PREVIEW_CAP = 500


def _read_csv(filename: str) -> pd.DataFrame:
    path = VALIDATION_DIR / filename
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"{filename} not found under {VALIDATION_DIR}")
    return pd.read_csv(path)


def list_studies() -> dict:
    df = _read_csv(STUDY_COMPARABILITY_CSV)
    return {"studies": df.astype(object).where(pd.notnull(df), None).to_dict(orient="records")}


def get_study(study: str) -> dict:
    if study not in STUDY_FILES:
        raise HTTPException(status_code=404, detail=f"unknown study {study!r}; choose one of {sorted(STUDY_FILES)}")
    df = _read_csv(STUDY_FILES[study])
    total_rows = len(df)
    windowed = df.head(_PREVIEW_CAP)
    return {
        "study": study,
        "total_rows": total_rows,
        "returned_rows": len(windowed),
        "columns": list(windowed.columns),
        "rows": windowed.astype(object).where(pd.notnull(windowed), None).to_dict(orient="records"),
    }
