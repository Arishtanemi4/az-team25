"""scoring_service -- wraps scoring/score.py::score_panel behind POST /rank. Runs standalone on
its own port; gene_service is the separate process that handles gene search and filter lookups.
See docker-compose.yml for how the two are wired together."""

from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from scoring_app.config import CORS_ORIGINS
from scoring_app.routers import ranking
from scoring_app.services.ranking_service import RankingService


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Loaded once at process startup -- score_panel's own large per-gene evidence tables are
    # still streamed fresh per query inside score.py, only the small reference tables are cached
    # here.
    app.state.ranking_service = RankingService()
    yield


app = FastAPI(title="CellLineSelector scoring_service", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_methods=["POST"],
    allow_headers=["*"],
)

app.include_router(ranking.router)


@app.get("/health")
def health():
    return {"status": "ok"}
