"""Backs GET /data/schema, GET /data/tables, GET /data/tables/{table}/preview."""

from fastapi import HTTPException

from data_app.services import table_io
from data_app.services.table_registry import ALLOWED_TABLES, TABLE_REGISTRY, column_role


def _row_count(table: str, manifest_counts: dict) -> int:
    return int(manifest_counts.get(f"{table}.csv", 0))


def get_schema() -> dict:
    manifest_counts = table_io.read_manifest_counts()
    tables = []
    for table in ALLOWED_TABLES:
        entry = TABLE_REGISTRY[table]
        dtypes = table_io.read_dtypes(table)
        columns = [
            {"name": name, "role": column_role(table, name), "dtype": dtypes.get(name, "unknown")}
            for name in table_io.read_columns(table)
        ]
        tables.append({
            "name": table,
            "grain": entry["grain"],
            "row_count": _row_count(table, manifest_counts),
            "columns": columns,
            "joins": entry["joins"],
        })
    return {"tables": tables}


def get_tables() -> dict:
    manifest_counts = table_io.read_manifest_counts()
    return {
        "tables": [
            {"name": table, "row_count": _row_count(table, manifest_counts)}
            for table in ALLOWED_TABLES
        ]
    }


def get_preview(table: str, limit: int, offset: int) -> dict:
    if table not in ALLOWED_TABLES:
        raise HTTPException(status_code=404, detail=f"unknown table {table!r}")
    manifest_counts = table_io.read_manifest_counts()
    rows, columns = table_io.preview_table(table, limit=limit, offset=offset)
    return {
        "table": table,
        "grain": TABLE_REGISTRY[table]["grain"],
        "total_rows": _row_count(table, manifest_counts),
        "returned_rows": len(rows),
        "columns": columns,
        "rows": rows,
    }
