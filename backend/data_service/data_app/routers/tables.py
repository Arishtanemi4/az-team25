"""GET /data/schema, GET /data/tables, GET /data/tables/{table}/preview -- schema and bounded
tabular access over the 14 `data/processed/` tables (PRODUCT_SURFACE.md SS3.1)."""

from fastapi import APIRouter, Query

from data_app.schemas.tables import PreviewResponse, SchemaResponse, TablesResponse
from data_app.services import table_schema_service

router = APIRouter(prefix="/data", tags=["data"])


@router.get("/schema", response_model=SchemaResponse)
def get_schema():
    return table_schema_service.get_schema()


@router.get("/tables", response_model=TablesResponse)
def get_tables():
    return table_schema_service.get_tables()


@router.get("/tables/{table}/preview", response_model=PreviewResponse)
def preview_table(
    table: str,
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    return table_schema_service.get_preview(table, limit=limit, offset=offset)
