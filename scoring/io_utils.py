import os

import pandas as pd


def read_table(csv_path: str) -> pd.DataFrame:
    parquet_path = csv_path.replace(".csv", ".parquet")
    if os.path.exists(parquet_path):
        df = pd.read_parquet(parquet_path, engine="pyarrow")
        print(f"[io_utils] read_table: {parquet_path} (parquet, {len(df)} rows)")
        return df

    df = pd.read_csv(csv_path)
    print(f"[io_utils] read_table: {csv_path} (csv fallback -- no parquet mirror found, {len(df)} rows)")
    return df


def read_gene_filtered(csv_path: str, ensembl_ids, chunksize: int = 1_000_000) -> pd.DataFrame:
    wanted = set(ensembl_ids)
    parquet_path = csv_path.replace(".csv", ".parquet")
    if os.path.exists(parquet_path):
        df = pd.read_parquet(
            parquet_path, engine="pyarrow", filters=[("ensembl_id", "in", list(wanted))],
        )
        print(f"[io_utils] read_gene_filtered: {parquet_path} (parquet, {len(df)} rows)")
        return df

    matches = []
    for chunk in pd.read_csv(csv_path, chunksize=chunksize):
        matches.append(chunk[chunk["ensembl_id"].isin(wanted)])
    df = pd.concat(matches, ignore_index=True) if matches else pd.DataFrame()
    print(f"[io_utils] read_gene_filtered: {csv_path} (csv fallback -- no parquet mirror found, {len(df)} rows)")
    return df
