"""Response models for /data/schema, /data/tables, /data/tables/{table}/preview."""

from pydantic import BaseModel, Field


class ColumnInfo(BaseModel):
    name: str
    role: str  # join_key | value | state | context
    dtype: str


class JoinInfo(BaseModel):
    to: str
    on: str
    cardinality: str


class TableInfo(BaseModel):
    name: str
    grain: str
    row_count: int
    columns: list[ColumnInfo]
    joins: list[JoinInfo]


class SchemaResponse(BaseModel):
    tables: list[TableInfo]


class TableSummary(BaseModel):
    name: str
    row_count: int


class TablesResponse(BaseModel):
    tables: list[TableSummary]


class PreviewResponse(BaseModel):
    table: str
    grain: str
    total_rows: int
    returned_rows: int
    columns: list[str]
    rows: list[dict]
