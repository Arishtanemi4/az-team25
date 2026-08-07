import json

import io_utils

DATA_DIR = "data/processed"
GENE_ROLE_REFERENCE_PATH = "preprocessing/resources/gene_role_reference.csv"
UNRESOLVED_SYMBOLS_PATH = "preprocessing/resources/unresolved_gene_symbols.json"

COVERAGE_STATE_COLUMNS = [
    "expression_rna_state", "expression_rna_hpa_state", "expression_rna_geo_state",
    "protein_state", "dependency_state", "copy_number_state",
    "mutations_state", "fusions_state",
]

LAYER_UNITS = {
    "rna": "log2(TPM+1)",
    "protein": "z-score",
    "dependency": "Chronos gene effect (dependency strength = -score)",
    "copy_number": "relative copy number (diploid ~= 1.0)",
    "mutation": "categorical -- VEP impact / driver / hotspot rule (desirability.MUTATION_D_TABLE)",
    "fusion": "categorical -- event presence, binary carve-out (desirability.FUSION_D_TABLE)",
}


def _load_unresolved_symbols(path):
    try:
        with open(path) as f:
            data = json.load(f)
    except FileNotFoundError:
        return {}
    return {layer: set(ids) for layer, ids in data.items() if not layer.startswith("_")}


def resolve_genes(tokens, gene_reference_df):
    resolved, ambiguous, unresolved = {}, [], []
    for token in tokens:
        by_id = gene_reference_df[gene_reference_df["ensembl_id"] == token]
        if len(by_id) == 1:
            resolved[token] = token
            continue

        by_symbol = gene_reference_df[gene_reference_df["symbol"] == token]
        if len(by_symbol) == 1:
            resolved[token] = by_symbol.iloc[0]["ensembl_id"]
        elif len(by_symbol) > 1:
            ambiguous.append(token)
        else:
            unresolved.append(token)
    return resolved, ambiguous, unresolved


def resolve_genes_or_raise(tokens, gene_reference_df):
    resolved, ambiguous, unresolved = resolve_genes(tokens, gene_reference_df)
    if ambiguous:
        raise ValueError(f"Ambiguous gene symbol(s) -- specify by Ensembl ID instead: {ambiguous}")
    if unresolved:
        raise ValueError(f"Gene(s) not found in gene_reference: {unresolved}")
    return resolved

STAGE_LABELS = {
    "resolve_genes": "Resolving requested genes against the reference table",
    "filter_candidates": "Filtering candidate cell lines by lineage, disease, and QC criteria",
    "load_tables": "Loading expression, dependency, mutation, copy-number, and fusion evidence tables",
    "correlation": "Computing inter-gene correlation weights across the reference expression panel",
    "score_candidates": "Scoring each candidate cell line against the requested gene panel",
    "rank": "Ranking candidates and splitting results by confidence tier",
    "similarity": "Identifying similar backup cell lines for top-ranked results",
}


COVERAGE_FILTER_COLUMNS = ["ModelID", "msi_high", "ploidy", "metabolomics_available", "mirna_available"]


def filter_candidates(cell_lines_df, coverage_df=None, lineage=None, primary_disease=None,
                       exclude_problematic=False, msi_high=None, ploidy_min=None, ploidy_max=None,
                       require_metabolomics=False, require_mirna=False):
    needs_coverage = (
        msi_high is not None or ploidy_min is not None or ploidy_max is not None
        or require_metabolomics or require_mirna
    )
    if needs_coverage and coverage_df is None:
        raise ValueError(
            "filter_candidates: msi_high/ploidy_min/ploidy_max/require_metabolomics/"
            "require_mirna need coverage_df (data/processed/coverage.csv) -- none was given."
        )

    subset = cell_lines_df
    if lineage is not None:
        subset = subset[subset["lineage"] == lineage]
    if primary_disease is not None:
        subset = subset[subset["primary_disease"] == primary_disease]
    if exclude_problematic:
        subset = subset[~subset["is_problematic"].astype(bool)]

    if needs_coverage:
        subset = subset.merge(coverage_df[COVERAGE_FILTER_COLUMNS], on="ModelID", how="left")
        if msi_high is not None:
            subset = subset[subset["msi_high"] == msi_high]
        if ploidy_min is not None:
            subset = subset[subset["ploidy"] >= ploidy_min]
        if ploidy_max is not None:
            subset = subset[subset["ploidy"] <= ploidy_max]
        if require_metabolomics:
            subset = subset[subset["metabolomics_available"].astype(bool)]
        if require_mirna:
            subset = subset[subset["mirna_available"].astype(bool)]

    return subset, len(subset)


def _load_filtered(path, genes, model_ids, chunksize=1_000_000):
    model_ids = set(model_ids)
    df = io_utils.read_gene_filtered(path, genes, chunksize=chunksize)
    return df[df["ModelID"].isin(model_ids)].reset_index(drop=True)
