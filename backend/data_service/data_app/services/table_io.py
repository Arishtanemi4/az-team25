"""Read-only access to `data/processed/*.parquet` (preferring the mirror, falling back to CSV)
for data_service. Deliberately duplicated from `scoring/io_utils.py` rather than imported --
`data_service` never imports `scoring/` (PRODUCT_SURFACE.md SS3), so the two stay independently
deployable microservices, same convention `gene_service`'s `gene_reference_service.py` already
follows for the same reason.
"""

import json
import os

import pandas as pd
import pyarrow.parquet as pq

from data_app.config import BUILD_MANIFEST_JSON, DATA_DIR
from data_app.services.table_registry import ALLOWED_TABLES

_PREVIEW_BATCH_SIZE = 50_000


def _csv_path(table: str) -> str:
    if table not in ALLOWED_TABLES:
        raise ValueError(f"unknown table {table!r}")
    return str(DATA_DIR / f"{table}.csv")


def _parquet_path(table: str) -> str:
    return _csv_path(table).replace(".csv", ".parquet")


def table_exists(table: str) -> bool:
    return os.path.exists(_parquet_path(table)) or os.path.exists(_csv_path(table))


def read_manifest_counts() -> dict:
    """`build_manifest.json`'s per-table row counts -- the authoritative total row count, so a
    preview or a relationship query never has to scan a 79-million-row table just to report how
    big it really is."""
    if not os.path.exists(BUILD_MANIFEST_JSON):
        return {}
    with open(BUILD_MANIFEST_JSON, encoding="utf-8") as fh:
        manifest = json.load(fh)
    return manifest.get("counts", {})


def read_columns(table: str) -> list[str]:
    """Column names only, no data read -- Parquet footer metadata is enough."""
    parquet_path = _parquet_path(table)
    if os.path.exists(parquet_path):
        return pq.ParquetFile(parquet_path).schema_arrow.names
    return list(pd.read_csv(_csv_path(table), nrows=0).columns)


def read_dtypes(table: str) -> dict[str, str]:
    """Column name -> a short dtype label, read from Parquet's declared schema when the mirror
    exists (accurate, no data read) or inferred from a small CSV sample otherwise."""
    parquet_path = _parquet_path(table)
    if os.path.exists(parquet_path):
        schema = pq.ParquetFile(parquet_path).schema_arrow
        return {name: str(schema.field(name).type) for name in schema.names}
    sample = pd.read_csv(_csv_path(table), nrows=1000)
    return {col: str(dtype) for col, dtype in sample.dtypes.items()}


def _json_safe(df: pd.DataFrame) -> list[dict]:
    """pandas NaN is not valid JSON -- convert to None so FastAPI's encoder does not choke."""
    return df.astype(object).where(pd.notnull(df), None).to_dict(orient="records")


def preview_table(table: str, limit: int, offset: int) -> tuple[list[dict], list[str]]:
    """Returns up to `limit` rows starting at `offset`, reading only the Parquet row-group
    batches that overlap the requested window rather than materialising the whole table --
    the same access-shape reasoning as V6-10 (docs/plan/EXECUTE.md), scaled down to a bounded
    row window instead of a gene filter."""
    parquet_path = _parquet_path(table)
    want_end = offset + limit

    if os.path.exists(parquet_path):
        pf = pq.ParquetFile(parquet_path)
        columns = pf.schema_arrow.names
        frames = []
        seen = 0
        for batch in pf.iter_batches(batch_size=_PREVIEW_BATCH_SIZE):
            batch_start = seen
            batch_end = seen + batch.num_rows
            seen = batch_end
            if batch_end <= offset:
                continue  # entirely before the window -- skip without materialising to pandas
            batch_df = batch.to_pandas()
            local_start = max(0, offset - batch_start)
            local_end = min(batch.num_rows, want_end - batch_start)
            frames.append(batch_df.iloc[local_start:local_end])
            if batch_end >= want_end:
                break
        window = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns)
        return _json_safe(window), columns

    columns = list(pd.read_csv(_csv_path(table), nrows=0).columns)
    frames = []
    seen = 0
    for chunk in pd.read_csv(_csv_path(table), chunksize=_PREVIEW_BATCH_SIZE):
        chunk_start = seen
        chunk_end = seen + len(chunk)
        seen = chunk_end
        if chunk_end <= offset:
            continue
        local_start = max(0, offset - chunk_start)
        local_end = min(len(chunk), want_end - chunk_start)
        frames.append(chunk.iloc[local_start:local_end])
        if chunk_end >= want_end:
            break
    window = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns)
    return _json_safe(window), columns


def read_full_column(table: str, column: str) -> list:
    """Reads one column's distinct values only -- used to size gene-intersection sets before a
    full read, without pulling every other column of a 79-million-row table into memory."""
    parquet_path = _parquet_path(table)
    if os.path.exists(parquet_path):
        arr = pq.read_table(parquet_path, columns=[column])[column]
        return arr.unique().to_pylist()
    return pd.read_csv(_csv_path(table), usecols=[column])[column].unique().tolist()


def read_full(table: str) -> pd.DataFrame:
    """Reads a whole (small) table -- only ever called on cell_lines/coverage/gene_reference,
    each a couple thousand to tens of thousands of rows."""
    parquet_path = _parquet_path(table)
    if os.path.exists(parquet_path):
        return pd.read_parquet(parquet_path, engine="pyarrow")
    return pd.read_csv(_csv_path(table))


def read_gene_filtered(table: str, ensembl_ids) -> pd.DataFrame:
    """Reads one of the large gene-keyed tables kept to just the given genes, via Parquet
    predicate pushdown on ensembl_id -- the same access shape `scoring/io_utils.py` proves (a
    one-gene query touches ~2,000 rows, not 79 million), duplicated here so data_service has no
    import dependency on scoring/."""
    wanted = list(set(ensembl_ids))
    parquet_path = _parquet_path(table)
    if os.path.exists(parquet_path):
        return pd.read_parquet(parquet_path, engine="pyarrow", filters=[("ensembl_id", "in", wanted)])

    matches = []
    for chunk in pd.read_csv(_csv_path(table), chunksize=1_000_000):
        matches.append(chunk[chunk["ensembl_id"].isin(wanted)])
    return pd.concat(matches, ignore_index=True) if matches else pd.DataFrame()
