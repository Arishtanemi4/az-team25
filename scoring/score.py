import json

import pandas as pd

import desirability
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


def load_tables_for_query(query_genes, model_ids, gene_reference_df, data_dir=DATA_DIR,
                           rna_constants_path=desirability.RNA_CONSTANTS_PATH,
                           extended_constants_path=desirability.EXTENDED_CONSTANTS_PATH,
                           essentiality_constants_path=desirability.ESSENTIALITY_CONSTANTS_PATH,
                           gene_role_path=GENE_ROLE_REFERENCE_PATH,
                           unresolved_symbols_path=UNRESOLVED_SYMBOLS_PATH,
                           coverage_df=None, rna_constants=None, extended_constants=None,
                           essentiality_constants=None, gene_role_reference=None,
                           unresolved_symbols=None, preloaded_layer_frames=None):
    if preloaded_layer_frames is not None:
        def _filter_preloaded(name):
            df = preloaded_layer_frames[name]
            if df.empty:
                return df
            mask = df["ensembl_id"].isin(query_genes) & df["ModelID"].isin(model_ids)
            return df[mask]

        expr = _filter_preloaded("expression_rna")
        hpa = _filter_preloaded("expression_rna_hpa")
        geo = _filter_preloaded("expression_rna_geo")
        protein = _filter_preloaded("protein")
        dependency = _filter_preloaded("dependency")
        copy_number = _filter_preloaded("copy_number")
        mutations = _filter_preloaded("mutations")
        fusions = _filter_preloaded("fusions")
    else:
        expr = _load_filtered(f"{data_dir}/expression_rna.csv", query_genes, model_ids)
        hpa = _load_filtered(f"{data_dir}/expression_rna_hpa.csv", query_genes, model_ids)
        geo = _load_filtered(f"{data_dir}/expression_rna_geo.csv", query_genes, model_ids)
        protein = _load_filtered(f"{data_dir}/protein.csv", query_genes, model_ids)
        dependency = _load_filtered(f"{data_dir}/dependency.csv", query_genes, model_ids)
        copy_number = _load_filtered(f"{data_dir}/copy_number.csv", query_genes, model_ids)
        mutations = _load_filtered(f"{data_dir}/mutations.csv", query_genes, model_ids)
        fusions = _load_filtered(f"{data_dir}/fusions.csv", query_genes, model_ids)

    coverage = coverage_df if coverage_df is not None else pd.read_csv(f"{data_dir}/coverage.csv")
    coverage = coverage[coverage["ModelID"].isin(model_ids)].copy()
    coverage["_layer_count"] = (coverage[COVERAGE_STATE_COLUMNS] != "not_assayed").sum(axis=1)

    protein_lookup = {}
    if not protein.empty:
        for row in protein.itertuples(index=False):
            protein_lookup[(row.ModelID, row.ensembl_id)] = (row.zscore, row.detected)

    return {
        "expression_rna": expr.set_index(["ModelID", "ensembl_id"])["log2tpm1"].to_dict() if not expr.empty else {},
        "expression_rna_hpa": hpa.set_index(["ModelID", "ensembl_id"])["log2ntpm1"].to_dict() if not hpa.empty else {},
        "expression_rna_geo": geo.set_index(["ModelID", "ensembl_id"])["log1p_expr"].to_dict() if not geo.empty else {},
        "protein": protein_lookup,
        "dependency": dependency.set_index(["ModelID", "ensembl_id"])["dependency_score"].to_dict() if not dependency.empty else {},
        "copy_number": copy_number.set_index(["ModelID", "ensembl_id"])["copy_number"].to_dict() if not copy_number.empty else {},
        "mutations": {k: v for k, v in mutations.groupby(["ModelID", "ensembl_id"])},
        "fusions": {k: v for k, v in fusions.groupby(["ModelID", "ensembl_id"])},
        "mutations_state": coverage.set_index("ModelID")["mutations_state"].to_dict(),
        "fusions_state": coverage.set_index("ModelID")["fusions_state"].to_dict(),
        "coverage_layer_count": coverage.set_index("ModelID")["_layer_count"].to_dict(),
        "gene_symbols": gene_reference_df.set_index("ensembl_id")["symbol"].to_dict(),
        "rna_constants": rna_constants if rna_constants is not None
            else desirability.load_rna_constants(rna_constants_path),
        "extended_constants": extended_constants if extended_constants is not None
            else desirability.load_extended_constants(extended_constants_path),
        "essentiality_constants": essentiality_constants if essentiality_constants is not None
            else desirability.load_essentiality_constants(essentiality_constants_path),

        "gene_class": (
            gene_role_reference if gene_role_reference is not None
            else pd.read_csv(gene_role_path).set_index("ensembl_id")["gene_class"].to_dict()
        ),

        "unresolved_symbols": (
            unresolved_symbols if unresolved_symbols is not None
            else _load_unresolved_symbols(unresolved_symbols_path)
        ),
    }
