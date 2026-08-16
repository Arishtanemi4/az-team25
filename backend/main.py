"""backend/main.py -- single FastAPI app combining gene_service, scoring_service, rag_service,
data_service, and research_service into one process on port 8000. Each service's own
gene_app/scoring_app/rag_app/data_app package still works standalone for its own test suite; this
file only imports their routers, schemas, and services and combines their startup lifespans -- no
route gets a new path prefix, so every URL is identical to when the four ran as separate processes.
data_service has no per-process state to load (every read happens at request time), so it needs no
lifespan entry. research_service is flat (not a routers/schemas/services `<x>_app` package like the
other four) -- it was moved in from extensions/research/, always on now, no env var required."""

import os
import sys
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

# Each service's inner package (gene_app/scoring_app/rag_app) uses bare imports rooted at its
# own package name, not a proper installed package, so its own directory must be on sys.path --
# same convention scoring_service/rag_service already use to make scoring/ and rag/ importable.
_BACKEND_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _BACKEND_DIR.parent
for _service_dir in ("gene_service", "scoring_service", "rag_service", "data_service", "research_service"):
    sys.path.insert(0, str(_BACKEND_DIR / _service_dir))
# research_service/service.py and research_export.py also import `backend.research_service.routes`
# by its full absolute path (routes/ has its own internal cross-references, and "routes" alone is
# not a safe bare name -- backend/research_service/tests/routes/ collides with it during pytest
# collection). That absolute import only resolves if the repo root itself is on sys.path, which
# happens automatically when this process is launched as `uvicorn backend.main:app` from the repo
# root (cwd), but not when launched as `cd backend && uvicorn main:app --reload` (cwd is backend/
# itself) -- computed via __file__, so it is correct in either case.
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from gene_app.config import CELL_LINES_CSV, GENE_REFERENCE_CSV
from gene_app.routers import filters as gene_filters
from gene_app.routers import genes as gene_genes
from gene_app.services.gene_reference_service import GeneReferenceService

from scoring_app.routers import ranking as scoring_ranking
from scoring_app.services.ranking_service import RankingService

from rag_app.routers import graph as rag_graph
from rag_app.routers import literature as rag_literature
from rag_app.routers import methodology as rag_methodology
from rag_app.routers import narration as rag_narration
from rag_app.routers import query_expansion as rag_query_expansion
from rag_app.services.rag_service import RagService

from data_app.routers import eda as data_eda
from data_app.routers import relationships as data_relationships
from data_app.routers import tables as data_tables
from data_app.routers import validation as data_validation

from api import router as research_router
from service import ResearchService

CORS_ORIGINS = os.environ.get(
    "BACKEND_CORS_ORIGINS", "http://localhost:5173,http://127.0.0.1:5173"
).split(",")


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Loaded once at process startup, mirroring what each service's own lifespan did standalone.
    app.state.gene_reference_service = GeneReferenceService(GENE_REFERENCE_CSV, CELL_LINES_CSV)
    app.state.research_service = ResearchService()
    app.state.ranking_service = RankingService(
        result_callback=app.state.research_service.capture_native_result
    )
    app.state.rag_service = RagService()
    yield


app = FastAPI(title="CellLineSelector backend", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)

app.include_router(gene_genes.router)
app.include_router(gene_filters.router)
app.include_router(scoring_ranking.router)
app.include_router(rag_narration.router)
app.include_router(rag_methodology.router)
app.include_router(rag_query_expansion.router)
app.include_router(rag_literature.router)
app.include_router(rag_graph.router)
app.include_router(data_tables.router)
app.include_router(data_relationships.router)
app.include_router(data_eda.router)
app.include_router(data_validation.router)
app.include_router(research_router)


@app.get("/health")
def health():
    return {"status": "ok"}
