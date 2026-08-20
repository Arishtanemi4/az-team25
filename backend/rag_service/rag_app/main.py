"""rag_service -- wraps rag/'s five components (result narrator, methodology Q&A agent,
query-expansion agent, literature agent, grounding verifier) behind four POST endpoints. Runs
standalone on its own port; gene_service and scoring_service are the separate processes that
handle gene lookups and the ranking math. See docker-compose.yml for how the three are wired
together."""

from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

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


@app.exception_handler(Exception)
async def unhandled_exception_handler(request: Request, exc: Exception):
    # Safety net: FastAPI/Starlette's own default for an uncaught exception is a plain-text 500
    # body, which the frontend's fetch clients (postJson in ragServiceClient.ts) can't parse for
    # a `detail` message -- they'd fall back to a generic "Request failed (500)". Every router
    # here already catches the specific exceptions it expects (RuntimeError, etc.); this only
    # fires for whatever gap remains, and still gives the UI something readable.
    return JSONResponse(status_code=500, content={"detail": f"{type(exc).__name__}: {exc}"})
