"""The 14 `data/processed/` tables, described once so every data_service endpoint (schema,
preview, relationships) reads the same facts instead of re-deriving them. This is domain
knowledge -- grain, join keys, which column is *the* measurement -- that cannot be recovered
from the CSV/Parquet header alone, so it is hand-maintained here, not inferred.

`ModelID` and `ensembl_id` are the **only two automatic join keys**
(`docs/PROJECT_ARCHITECTURE.md` SS4); every table below marks them (where present) with
role="join_key" and nothing else gets that role.
"""

# name -> set of columns that are role="value" (the table's real measurement, as opposed to a
# join key, a state flag, or contextual metadata).
_VALUE_COLUMNS = {
    "cell_lines": set(),
    "coverage": set(),
    "gene_reference": set(),
    "expression_rna": {"log2tpm1"},
    "expression_rna_hpa": {"log2ntpm1"},
    "expression_rna_geo": {"log1p_expr"},
    "protein": {"zscore"},
    "mutations": set(),
    "fusions": set(),
    "dependency": {"dependency_score"},
    "copy_number": {"copy_number"},
    "genome_signatures": set(),
    "mirna": {"log1p_value"},
    "metabolomics": {"value"},
}

# name -> set of columns that are role="state" -- a controlled-vocabulary assay-state flag, as
# opposed to a continuous value or free-form context. coverage.csv's *_state/*_available columns
# and protein.csv's `detected` are the only ones in this project (`docs/plan/PARAMETERS.md`).
_STATE_COLUMNS = {
    "coverage": {
        "expression_rna_state", "expression_rna_hpa_state", "expression_rna_geo_state",
        "protein_state", "dependency_state", "copy_number_state", "genome_signatures_state",
        "mutations_state", "fusions_state", "metabolomics_available", "mirna_available",
        "wgs_available",
    },
    "protein": {"detected"},
}

# The coverage.csv columns representing a per-model, per-layer assay-state flag -- used by
# eda_service.lineage_coverage to break the global coverage_state_counts down per lineage.
COVERAGE_LAYER_COLUMNS = sorted(_STATE_COLUMNS["coverage"])

# name -> {grain, joins}. `joins` lists every other table this one joins to via ModelID and/or
# ensembl_id, and the cardinality from this table's own grain looking outward.
TABLE_REGISTRY = {
    "cell_lines": {
        "grain": "one row per ModelID",
        "joins": [],
    },
    "gene_reference": {
        "grain": "one row per ensembl_id",
        "joins": [],
    },
    "coverage": {
        "grain": "one row per ModelID",
        "joins": [{"to": "cell_lines", "on": "ModelID", "cardinality": "one_to_one"}],
    },
    "expression_rna": {
        "grain": "one row per (ModelID, ensembl_id)",
        "joins": [
            {"to": "cell_lines", "on": "ModelID", "cardinality": "many_to_one"},
            {"to": "gene_reference", "on": "ensembl_id", "cardinality": "many_to_one"},
        ],
    },
    "expression_rna_hpa": {
        "grain": "one row per (ModelID, ensembl_id)",
        "joins": [
            {"to": "cell_lines", "on": "ModelID", "cardinality": "many_to_one"},
            {"to": "gene_reference", "on": "ensembl_id", "cardinality": "many_to_one"},
        ],
    },
    "expression_rna_geo": {
        "grain": "one row per (ModelID, ensembl_id)",
        "joins": [
            {"to": "cell_lines", "on": "ModelID", "cardinality": "many_to_one"},
            {"to": "gene_reference", "on": "ensembl_id", "cardinality": "many_to_one"},
        ],
    },
    "protein": {
        "grain": "one row per (ModelID, ensembl_id)",
        "joins": [
            {"to": "cell_lines", "on": "ModelID", "cardinality": "many_to_one"},
            {"to": "gene_reference", "on": "ensembl_id", "cardinality": "many_to_one"},
        ],
    },
    "dependency": {
        "grain": "one row per (ModelID, ensembl_id)",
        "joins": [
            {"to": "cell_lines", "on": "ModelID", "cardinality": "many_to_one"},
            {"to": "gene_reference", "on": "ensembl_id", "cardinality": "many_to_one"},
        ],
    },
    "copy_number": {
        "grain": "one row per (ModelID, ensembl_id) sequencing call",
        "joins": [
            {"to": "cell_lines", "on": "ModelID", "cardinality": "many_to_one"},
            {"to": "gene_reference", "on": "ensembl_id", "cardinality": "many_to_one"},
        ],
    },
    "mutations": {
        "grain": "one row per called variant (ModelID, ensembl_id, position)",
        "joins": [
            {"to": "cell_lines", "on": "ModelID", "cardinality": "many_to_one"},
            {"to": "gene_reference", "on": "ensembl_id", "cardinality": "many_to_one"},
        ],
    },
    "fusions": {
        "grain": "one row per called fusion event (ModelID, gene pair)",
        "joins": [
            {"to": "cell_lines", "on": "ModelID", "cardinality": "many_to_one"},
            {"to": "gene_reference", "on": "ensembl_id", "cardinality": "many_to_one"},
        ],
    },
    "genome_signatures": {
        "grain": "one row per ModelID",
        "joins": [{"to": "cell_lines", "on": "ModelID", "cardinality": "one_to_one"}],
    },
    "mirna": {
        "grain": "one row per (ModelID, mirna_name)",
        "joins": [{"to": "cell_lines", "on": "ModelID", "cardinality": "many_to_one"}],
    },
    "metabolomics": {
        "grain": "one row per (ModelID, metabolite_name)",
        "joins": [{"to": "cell_lines", "on": "ModelID", "cardinality": "many_to_one"}],
    },
}

# The allow-list `{table}` path parameters are validated against -- never used to build a path
# straight from user input.
ALLOWED_TABLES = tuple(sorted(TABLE_REGISTRY))


def column_role(table: str, column: str) -> str:
    """Classifies one column of one table as join_key / value / state / context. join_key is
    the same two names everywhere (ModelID, ensembl_id); value/state are the per-table sets
    above; everything else is context."""
    if column in ("ModelID", "ensembl_id"):
        return "join_key"
    if column in _VALUE_COLUMNS.get(table, set()):
        return "value"
    if column in _STATE_COLUMNS.get(table, set()):
        return "state"
    return "context"
