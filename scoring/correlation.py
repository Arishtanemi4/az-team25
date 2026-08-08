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
