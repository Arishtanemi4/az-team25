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
