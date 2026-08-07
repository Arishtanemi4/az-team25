import json

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
