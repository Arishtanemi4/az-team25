import shutil
import tempfile
import time
from pathlib import Path

import pandas as pd

import io_utils
import score

DATA_DIR = Path("data/processed")
BENCHMARK_TABLES = ["expression_rna.csv", "dependency.csv", "mutations.csv"]
SAMPLE_SYMBOLS = ["EGFR", "KRAS", "TP53"]


def _time_read(csv_path: Path, ensembl_ids: list[str]) -> float:
    start = time.perf_counter()
    io_utils.read_gene_filtered(str(csv_path), ensembl_ids)
    return time.perf_counter() - start


def benchmark_table(filename: str, ensembl_ids: list[str]) -> None:
    csv_path = DATA_DIR / filename
    parquet_path = csv_path.with_suffix(".parquet")
    if not parquet_path.exists():
        print(f"{filename}: no .parquet mirror found -- run preprocessing/preprocess.py first, skipping")
        return

    with_parquet_seconds = _time_read(csv_path, ensembl_ids)

    with tempfile.TemporaryDirectory() as scratch_dir:
        moved_path = Path(scratch_dir) / parquet_path.name
        shutil.move(str(parquet_path), str(moved_path))
        try:
            csv_fallback_seconds = _time_read(csv_path, ensembl_ids)
        finally:
            shutil.move(str(moved_path), str(parquet_path))

    speedup = csv_fallback_seconds / with_parquet_seconds if with_parquet_seconds > 0 else float("inf")
    print(
        f"{filename:28s} parquet={with_parquet_seconds:6.2f}s  csv={csv_fallback_seconds:6.2f}s  "
        f"speedup={speedup:5.1f}x"
    )


def main() -> None:
    gene_reference = pd.read_csv(DATA_DIR / "gene_reference.csv")
    resolved, ambiguous, unresolved = score.resolve_genes(SAMPLE_SYMBOLS, gene_reference)
    assert not ambiguous and not unresolved, (ambiguous, unresolved)
    ensembl_ids = list(resolved.values())

    print(f"Benchmarking with sample query genes: {SAMPLE_SYMBOLS} -> {ensembl_ids}\n")
    for filename in BENCHMARK_TABLES:
        benchmark_table(filename, ensembl_ids)


if __name__ == "__main__":
    main()
