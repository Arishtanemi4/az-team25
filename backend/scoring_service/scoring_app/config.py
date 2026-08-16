"""Paths and settings for scoring_service. scoring/score.py's own default paths
(data/processed/, preprocessing/resources/..., scoring/resources/...) are repo-root-relative,
which only holds if the process is launched with the repo root as its working directory. This
module resolves the same paths as absolutes instead, so scoring works the same whether launched
locally (`uvicorn app.main:app` from backend/scoring_service/) or inside a Docker container.
"""

import os
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[3]

DATA_DIR = Path(os.environ.get("SCORING_SERVICE_DATA_DIR", _REPO_ROOT / "data" / "processed"))
GENE_REFERENCE_CSV = DATA_DIR / "gene_reference.csv"
CELL_LINES_CSV = DATA_DIR / "cell_lines.csv"

RNA_CONSTANTS_PATH = Path(
    os.environ.get(
        "SCORING_SERVICE_RNA_CONSTANTS_PATH",
        _REPO_ROOT / "scoring" / "resources" / "desirability_constants.json",
    )
)
EXTENDED_CONSTANTS_PATH = Path(
    os.environ.get(
        "SCORING_SERVICE_EXTENDED_CONSTANTS_PATH",
        _REPO_ROOT / "scoring" / "resources" / "desirability_constants_extended.json",
    )
)
ESSENTIALITY_CONSTANTS_PATH = Path(
    os.environ.get(
        "SCORING_SERVICE_ESSENTIALITY_CONSTANTS_PATH",
        _REPO_ROOT / "scoring" / "resources" / "common_essential_genes.json",
    )
)
GENE_ROLE_PATH = Path(
    os.environ.get(
        "SCORING_SERVICE_GENE_ROLE_PATH",
        _REPO_ROOT / "preprocessing" / "resources" / "gene_role_reference.csv",
    )
)
UNRESOLVED_SYMBOLS_PATH = Path(
    os.environ.get(
        "SCORING_SERVICE_UNRESOLVED_SYMBOLS_PATH",
        _REPO_ROOT / "preprocessing" / "resources" / "unresolved_gene_symbols.json",
    )
)

DEFAULT_TOP_K = 10

CORS_ORIGINS = os.environ.get(
    "SCORING_SERVICE_CORS_ORIGINS", "http://localhost:5173,http://127.0.0.1:5173"
).split(",")
