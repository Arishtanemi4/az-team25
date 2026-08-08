LAYER_WEIGHTS = {
    "rna": 4 / 15,
    "protein": 2 / 15,
    "fusion": 1 / 15,
    "copy_number": 2 / 15,
    "mutation": 2 / 15,
    "dependency": 4 / 15,
}


def combine_gene_evidence(layer_d_values, weights=None):
    active_weights = weights if weights is not None else LAYER_WEIGHTS
    present = {
        layer: d for layer, d in layer_d_values.items()
        if d is not None and layer in active_weights
    }
    if not present:
        return None

    total_weight = sum(active_weights[layer] for layer in present)
    if total_weight == 0:
        return None
    weighted_sum = sum(active_weights[layer] * d for layer, d in present.items())
    return weighted_sum / total_weight
