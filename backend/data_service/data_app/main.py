"""data_service -- read-only access to the processed data tables (schema, bounded table
preview, relationship/EDA aggregates, validation results). Runs standalone on its own port for
its own test suite; backend/main.py mounts the same routers into the combined app. Never scores,
re-ranks, or imports `scoring/` (PRODUCT_SURFACE.md SS3) -- it only shows and explains data."""

from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from data_app.config import CORS_ORIGINS
from data_app.routers import eda, relationships, tables, validation


@asynccontextmanager
async def lifespan(app: FastAPI):
    yield


app = FastAPI(title="CellLineSelector data_service", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ORIGINS,
    allow_methods=["GET"],
    allow_headers=["*"],
)

app.include_router(tables.router)
app.include_router(relationships.router)
app.include_router(eda.router)
app.include_router(validation.router)


@app.get("/health")
def health():
    return {"status": "ok"}
