"""Loads the small reference tables (gene_reference.csv, cell_lines.csv) once at startup and
answers search/filter lookups from memory. This service never touches scoring math -- it only
helps the frontend build valid queries (which gene tokens exist, which lineages/diseases exist).
"""

import os

import pandas as pd


def _read_table_or_parquet(csv_path) -> pd.DataFrame:
    """Prefers the `.parquet` mirror preprocessing/preprocess.py writes alongside every
    data/processed/*.csv table (binary decode is faster than text parsing), falling back to the
    CSV when no mirror exists (a fresh clone; mirrors are gitignored). Deliberately duplicated
    from scoring/io_utils.py::read_table rather than imported -- gene_service has no other
    dependency on scoring/, and the two stay independently deployable microservices (root
    `_.md` SS9)."""
    csv_path = str(csv_path)
    parquet_path = csv_path.replace(".csv", ".parquet")
    if os.path.exists(parquet_path):
        df = pd.read_parquet(parquet_path, engine="pyarrow")
        print(f"[GeneReferenceService] read {parquet_path} (parquet, {len(df)} rows)")
        return df

    df = pd.read_csv(csv_path)
    print(f"[GeneReferenceService] read {csv_path} (csv fallback -- no parquet mirror found, {len(df)} rows)")
    return df


class GeneReferenceService:
    """Wraps gene_reference.csv and cell_lines.csv for search-as-you-type gene lookup and
    dropdown filter options. Loading both tables once in __init__ avoids re-reading these
    files from disk on every request -- they change only when preprocessing/ is re-run."""

    def __init__(self, gene_reference_csv, cell_lines_csv):
        self.gene_reference_df = _read_table_or_parquet(gene_reference_csv)
        self.cell_lines_df = _read_table_or_parquet(cell_lines_csv)

    def search_genes(self, query: str, limit: int = 20) -> list[dict]:
        """Case-insensitive substring match on gene symbol, falling back to ensembl_id, so a
        researcher can type either kind of token. Returns at most `limit` matches, symbol
        matches ranked before ensembl_id matches since researchers usually think in symbols.
        An empty query (e.g. clicking into an untouched field) returns the first `limit` genes
        alphabetically by symbol, so the dropdown never shows nothing on an empty click."""
        query = query.strip()
        if not query:
            return (
                self.gene_reference_df.sort_values("symbol")[["ensembl_id", "symbol"]]
                .head(limit)
                .to_dict(orient="records")
            )

        symbol_mask = self.gene_reference_df["symbol"].str.contains(query, case=False, na=False)
        id_mask = self.gene_reference_df["ensembl_id"].str.contains(query, case=False, na=False)

        symbol_matches = self.gene_reference_df[symbol_mask]
        id_only_matches = self.gene_reference_df[id_mask & ~symbol_mask]
        combined = pd.concat([symbol_matches, id_only_matches]).head(limit)

        return combined[["ensembl_id", "symbol"]].to_dict(orient="records")

    def list_lineages(self) -> list[str]:
        """Every distinct lineage in cell_lines.csv, sorted, for the lineage filter dropdown."""
        return sorted(self.cell_lines_df["lineage"].dropna().unique().tolist())

    def list_primary_diseases(self, lineage: str | None = None) -> list[str]:
        """Distinct primary_disease values, optionally scoped to one lineage so the disease
        dropdown can cascade from whatever lineage the researcher already picked."""
        df = self.cell_lines_df
        if lineage is not None:
            df = df[df["lineage"] == lineage]
        return sorted(df["primary_disease"].dropna().unique().tolist())
