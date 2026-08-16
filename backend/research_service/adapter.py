"""Read-only access to data/processed/ for the extension layer (extensions/plans/CONTRACTS.md C1).

Mirrors scoring/score.py::_load_filtered's Parquet-predicate-pushdown-first, chunked-CSV-fallback
read pattern (V6-10) rather than importing that private helper directly -- an extension may read
a documented baseline output, but "import a private helper merely to bypass a public contract" is
explicitly the thing extensions/README.md forbids, and `_load_filtered` is underscore-private.
Never writes to data/processed/, never re-runs preprocessing, never caches a whole table.
"""

import os
from pathlib import Path

import pandas as pd

# Every layer this adapter can read, mapped to its data/processed/ filename. One place to add a
# new layer, not a constant scattered across every caller.
LAYER_FILES = {
    "rna": "expression_rna.csv",
    "rna_hpa": "expression_rna_hpa.csv",
    "rna_geo": "expression_rna_geo.csv",
    "protein": "protein.csv",
    "dependency": "dependency.csv",
    "copy_number": "copy_number.csv",
    "mutations": "mutations.csv",
    "fusions": "fusions.csv",
    "genome_signatures": "genome_signatures.csv",
    "mirna": "mirna.csv",
    "metabolomics": "metabolomics.csv",
    "coverage": "coverage.csv",
}

# Layers whose grain is one row per (ModelID, ensembl_id). mutations/fusions are event tables and
# may legitimately carry several rows per (ModelID, ensembl_id) -- C1: "retain event-table
# multiplicity" -- so they are the one deliberate exclusion, never asserted unique.
SINGLE_ROW_LAYERS = frozenset(LAYER_FILES) - {"mutations", "fusions"}


def _resolve_data_dir(data_dir=None):
    """SCORING_SERVICE_DATA_DIR is the existing, already-approved read boundary (C1: "reuse
    verified frozen-compatible processed outputs read-only via existing SCORING_SERVICE_DATA_DIR
    and GENE_SERVICE_DATA_DIR") -- reused here rather than inventing a second env var for the
    same directory."""
    if data_dir is not None:
        return Path(data_dir)
    env_value = os.environ.get("SCORING_SERVICE_DATA_DIR")
    if env_value:
        return Path(env_value)
    return Path(__file__).resolve().parents[2] / "data" / "processed"


def _read_one_layer(path, gene_ids, model_ids):
    """Same Parquet-predicate-pushdown-first, chunked-CSV-fallback strategy
    scoring/score.py::_load_filtered uses (V6-10), reimplemented here rather than imported,
    since that function is private and this package may not depend on a private symbol."""
    gene_ids = set(gene_ids)
    model_ids = set(model_ids)
    parquet_path = str(path).replace(".csv", ".parquet")
    if os.path.exists(parquet_path):
        frame = pd.read_parquet(
            parquet_path, engine="pyarrow", filters=[("ensembl_id", "in", list(gene_ids))]
        )
        return frame[frame["ModelID"].isin(model_ids)].reset_index(drop=True)

    if not os.path.exists(path):
        raise FileNotFoundError(f"No processed table at {path} (or its Parquet mirror)")

    matches = []
    for chunk in pd.read_csv(path, chunksize=1_000_000):
        mask = chunk["ensembl_id"].isin(gene_ids) & chunk["ModelID"].isin(model_ids)
        matches.append(chunk[mask])
    return pd.concat(matches, ignore_index=True) if matches else pd.DataFrame()


def read_measurements(model_ids, gene_ids, layers, data_dir=None):
    """Reads a bounded slice of data/processed/ for exactly the requested models/genes/layers.

    Returns `{layer_name: DataFrame}`. An unknown layer name raises immediately -- C1's "report
    unmatched mappings; no alias guessing" -- rather than silently returning nothing for it. A
    layer whose processed table (and Parquet mirror) is entirely absent raises FileNotFoundError
    naming the layer, so a caller sees a missing input, never an empty frame indistinguishable
    from "measured and found nothing" (the not_assayed/measured_absent distinction this whole
    project is built around, extensions/BIOLOGICAL_CONTEXT_RESEARCH.md's preprocessing
    compatibility rules).

    For every single-row-grain layer, a duplicate (ModelID, ensembl_id) combination read back is
    a data-contract violation (preprocess.py itself asserts this uniqueness when it writes these
    tables) -- raised here rather than silently averaged or the first row kept, so a genuine
    upstream regression is never masked by this adapter.
    """
    unknown = set(layers) - set(LAYER_FILES)
    if unknown:
        raise ValueError(f"Unknown layer(s), not in LAYER_FILES: {sorted(unknown)}")

    base_dir = _resolve_data_dir(data_dir)
    result = {}
    for layer in layers:
        path = base_dir / LAYER_FILES[layer]
        frame = _read_one_layer(path, gene_ids, model_ids)
        if layer in SINGLE_ROW_LAYERS and not frame.empty:
            duplicate_mask = frame.duplicated(subset=["ModelID", "ensembl_id"], keep=False)
            if duplicate_mask.any():
                bad_keys = (
                    frame.loc[duplicate_mask, ["ModelID", "ensembl_id"]]
                    .drop_duplicates()
                    .to_records(index=False)
                    .tolist()
                )
                raise ValueError(
                    f"Layer '{layer}' returned duplicate (ModelID, ensembl_id) rows, which its "
                    f"own grain forbids: {bad_keys}"
                )
        result[layer] = frame
    return result


def read_model_level_availability(model_ids, layer, data_dir=None):
    """Whether `layer` has any row for each of `model_ids` -- for the three tables that are
    NOT gene-keyed (metabolomics: ModelID+metabolite_name+value; mirna: ModelID+miRNA+value;
    genome_signatures: ModelID + model-level scalars, no gene dimension at all per
    docs/PROJECT_ARCHITECTURE.md SS6's admission gate). `read_measurements` cannot be reused
    here -- its `_read_one_layer` filters on an `ensembl_id` column these three tables do not
    have. Returns `{model_id: bool}`. Raises FileNotFoundError, same as read_measurements, if
    the table is entirely absent -- never a silent False standing in for 'never checked'."""
    if layer not in ("metabolomics", "mirna", "genome_signatures"):
        raise ValueError(
            f"read_model_level_availability is only for non-gene-keyed layers "
            f"(metabolomics, mirna, genome_signatures), not {layer!r}"
        )
    base_dir = _resolve_data_dir(data_dir)
    path = base_dir / LAYER_FILES[layer]
    parquet_path = str(path).replace(".csv", ".parquet")
    if os.path.exists(parquet_path):
        frame = pd.read_parquet(parquet_path, engine="pyarrow", columns=["ModelID"])
    elif os.path.exists(path):
        frame = pd.read_csv(path, usecols=["ModelID"])
    else:
        raise FileNotFoundError(f"No processed table at {path} (or its Parquet mirror)")

    present_ids = set(frame["ModelID"])
    return {model_id: model_id in present_ids for model_id in model_ids}


def read_model_metadata(model_ids, data_dir=None):
    """Cell-line identity/context rows for exactly the requested ModelIDs -- the source of a C2
    snapshot's `models` section. Returns (matched_frame, unknown_ids); an unknown ModelID is
    reported to the caller rather than silently dropped."""
    base_dir = _resolve_data_dir(data_dir)
    path = base_dir / "cell_lines.csv"
    if not os.path.exists(path):
        raise FileNotFoundError(f"No cell_lines.csv at {path}")
    frame = pd.read_csv(path)
    requested = set(model_ids)
    matched = frame[frame["ModelID"].isin(requested)].reset_index(drop=True)
    unknown_ids = sorted(requested - set(matched["ModelID"]))
    return matched, unknown_ids
