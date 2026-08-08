import pandas as pd

import io_utils

EXPRESSION_RNA_PATH = "data/processed/expression_rna.csv"


def load_query_expression(ensembl_ids, expression_rna_path=EXPRESSION_RNA_PATH, chunksize=1_000_000):
    long = io_utils.read_gene_filtered(expression_rna_path, ensembl_ids, chunksize=chunksize)
    return long.pivot(index="ModelID", columns="ensembl_id", values="log2tpm1")


def compute_pairwise_spearman(wide_df):
    return wide_df.corr(method="spearman")


def compute_rho_bar(corr_matrix):
    m = len(corr_matrix)
    rho_bar = {}
    for gene in corr_matrix.columns:
        others = corr_matrix.loc[gene].drop(index=gene)
        mean_corr = others.mean() if m > 1 else None
        rho_bar[gene] = max(mean_corr, 0.0) if mean_corr is not None else None
    return rho_bar


def compute_weights(rho_bar, m):
    weights = {}
    for gene, rho in rho_bar.items():
        if rho is None:
            weights[gene] = 1.0
        else:
            vif = 1 + (m - 1) * rho
            weights[gene] = 1 / vif
    m_eff = sum(weights.values())
    return weights, m_eff


SMALL_QUERY_CAVEAT = (
    "CAMERA's variance-inflation correction was validated on gene sets of ~15-580 genes. "
    "This query has fewer genes than that, so rho_bar and m_eff should be read as an informed "
    "estimate of shared signal, not a fully validated statistical guarantee."
)


def get_correlation_weights(ensembl_ids, expression_rna_path=EXPRESSION_RNA_PATH):
    m = len(ensembl_ids)
    if m == 1:
        gene = ensembl_ids[0]
        return {
            "weights": {gene: 1.0},
            "rho_bar": {gene: None},
            "m_eff": 1.0,
            "note": "m=1: no correlation discount is possible, none was applied.",
        }

    wide_df = load_query_expression(ensembl_ids, expression_rna_path)
    corr_matrix = compute_pairwise_spearman(wide_df)
    rho_bar = compute_rho_bar(corr_matrix)
    weights, m_eff = compute_weights(rho_bar, m)

    for gene in ensembl_ids:
        if gene not in weights:
            weights[gene] = 1.0
            rho_bar[gene] = None
            m_eff += 1.0

    return {"weights": weights, "rho_bar": rho_bar, "m_eff": m_eff, "note": SMALL_QUERY_CAVEAT}
