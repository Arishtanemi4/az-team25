import pandas as pd


def _query_gene_matrix(ranked_results, tables, query_genes):
    data = {}
    for result in ranked_results:
        model_id = result["model_id"]
        data[model_id] = {
            gene: tables["expression_rna"].get((model_id, gene)) for gene in query_genes
        }
    return pd.DataFrame(data)
