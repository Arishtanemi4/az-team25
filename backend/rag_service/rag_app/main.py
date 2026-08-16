"""rag_service -- wraps rag/'s five components (result narrator, methodology Q&A agent,
query-expansion agent, literature agent, grounding verifier) behind four POST endpoints. Runs
standalone on its own port; gene_service and scoring_service are the separate processes that
handle gene lookups and the ranking math. See docker-compose.yml for how the three are wired
together."""

from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from rag_app.config import CORS_ORIGINS
from rag_app.routers import graph, literature, methodology, narration, query_expansion
from rag_app.services.rag_service import RagService


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Loaded once at process startup -- citations.build_registry() is cheap but there is no
    # reason to rebuild it on every request (app/services/rag_service.py).
    app.state.rag_service = RagService()
    yield


app = FastAPI(title="CellLineSelector rag_service", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_methods=["GET", "POST"],
    allow_headers=["*"],
)

app.include_router(narration.router)
app.include_router(methodology.router)
app.include_router(query_expansion.router)
app.include_router(literature.router)
app.include_router(graph.router)


@app.get("/health")
def health():
    return {"status": "ok"}
