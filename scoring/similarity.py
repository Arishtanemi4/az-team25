import pandas as pd


def _query_gene_matrix(ranked_results, tables, query_genes):
    data = {}
    for result in ranked_results:
        model_id = result["model_id"]
        data[model_id] = {
            gene: tables["expression_rna"].get((model_id, gene)) for gene in query_genes
        }
    return pd.DataFrame(data)


def attach_similar_lines(ranked_results, tables, query_genes, top_k=3):
    if len(query_genes) < 2 or len(ranked_results) < 2:
        for result in ranked_results:
            result["similar_lines"] = []
        return ranked_results

    matrix = _query_gene_matrix(ranked_results, tables, query_genes)
    corr = matrix.corr(method="spearman")

    for result in ranked_results:
        model_id = result["model_id"]
        others = []
        for other in ranked_results:
            other_id = other["model_id"]
            if other_id == model_id or other["D"] is None or other["D"] <= 0.0:
                continue
            correlation_value = corr.loc[model_id, other_id]
            if pd.isna(correlation_value):
                continue  # not enough shared gene evidence between these two lines
            others.append({
                "model_id": other_id, "D": other["D"], "correlation": float(correlation_value),
            })
        others.sort(key=lambda o: -o["correlation"])
        result["similar_lines"] = others[:top_k]

    return ranked_results
