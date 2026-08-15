"""Paths and settings for data_service, all overridable by environment variables so the same
code runs unchanged locally (repo-root CWD) and inside a Docker container (data mounted at
/app/data/processed, validation results at /app/validation) -- same convention as
gene_service/gene_app/config.py.
"""

import os
from pathlib import Path

# Resolve the repo root relative to this file (backend/data_service/data_app/config.py -> repo
# root is four levels up) so paths are correct no matter what directory the process was launched
# from.
_REPO_ROOT = Path(__file__).resolve().parents[3]

DATA_DIR = Path(os.environ.get("DATA_SERVICE_DATA_DIR", _REPO_ROOT / "data" / "processed"))
BUILD_MANIFEST_JSON = DATA_DIR / "build_manifest.json"
EDA_AGGREGATES_JSON = DATA_DIR / "eda_aggregates.json"

VALIDATION_DIR = Path(os.environ.get("DATA_SERVICE_VALIDATION_DIR", _REPO_ROOT / "validation"))

# The frontend's Vite dev server origin, plus a comma-separated override for other setups --
# only used when this service runs standalone (its own tests / a future split-out deployment);
# backend/main.py sets CORS once for the combined app.
CORS_ORIGINS = os.environ.get(
    "DATA_SERVICE_CORS_ORIGINS", "http://localhost:5173,http://127.0.0.1:5173"
).split(",")
