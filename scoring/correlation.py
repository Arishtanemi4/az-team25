import pandas as pd

import io_utils

EXPRESSION_RNA_PATH = "data/processed/expression_rna.csv"


def load_query_expression(ensembl_ids, expression_rna_path=EXPRESSION_RNA_PATH, chunksize=1_000_000):
    long = io_utils.read_gene_filtered(expression_rna_path, ensembl_ids, chunksize=chunksize)
    return long.pivot(index="ModelID", columns="ensembl_id", values="log2tpm1")


def compute_pairwise_spearman(wide_df):
    return wide_df.corr(method="spearman")
