"""gene_service -- reference/lookup data only (gene search, lineage/disease filters).
Runs standalone on its own port; scoring_service is the separate process that does the
actual ranking math. See docker-compose.yml for how the two are wired together."""

from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from gene_app.config import CELL_LINES_CSV, CORS_ORIGINS, GENE_REFERENCE_CSV
from gene_app.routers import filters, genes
from gene_app.services.gene_reference_service import GeneReferenceService


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Loaded once at process startup, not per-request -- these files change only when
    # preprocessing/ is re-run, not on every gene search.
    app.state.gene_reference_service = GeneReferenceService(GENE_REFERENCE_CSV, CELL_LINES_CSV)
    yield


app = FastAPI(title="CellLineSelector gene_service", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_methods=["GET"],
    allow_headers=["*"],
)

app.include_router(genes.router)
app.include_router(filters.router)


@app.get("/health")
def health():
    return {"status": "ok"}
