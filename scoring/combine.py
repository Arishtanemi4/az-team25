import math


def combine_across_genes(d_gene, weights):
    for ensembl_id, d in d_gene.items():
        if d == 0.0 and weights.get(ensembl_id, 0.0) > 0.0:
            return 0.0, {"ensembl_id": ensembl_id}

    total_weight = sum(weights[g] for g in d_gene)
    log_sum = sum(weights[g] * math.log(d_gene[g]) for g in d_gene)
    D = math.exp(log_sum / total_weight)
    return D, None
