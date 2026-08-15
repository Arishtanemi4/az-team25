"""Paths and settings for gene_service, all overridable by environment variables so the
same code runs unchanged locally (repo-root CWD) and inside a Docker container (data mounted
at /app/data/processed).
"""

import os
from pathlib import Path

# Resolve the repo root relative to this file (backend/gene_service/app/config.py -> repo root
# is four levels up) so paths are correct no matter what directory the process was launched from.
_REPO_ROOT = Path(__file__).resolve().parents[3]

DATA_DIR = Path(os.environ.get("GENE_SERVICE_DATA_DIR", _REPO_ROOT / "data" / "processed"))
GENE_REFERENCE_CSV = DATA_DIR / "gene_reference.csv"
CELL_LINES_CSV = DATA_DIR / "cell_lines.csv"

# The frontend's Vite dev server origin, plus a comma-separated override for other setups.
CORS_ORIGINS = os.environ.get(
    "GENE_SERVICE_CORS_ORIGINS", "http://localhost:5173,http://127.0.0.1:5173"
).split(",")
