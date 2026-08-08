def check_hpa_corroboration(rna_values, hpa_values, rna_constants):
    shared = [g for g in rna_values if g in hpa_values]
    if not shared:
        return None

    for gene in shared:
        entry = rna_constants["global"].get(gene)
        if entry is None:
            continue
        floor = entry["L"]
        rna_on = rna_values[gene] >= floor
        hpa_on = hpa_values[gene] > 0  # HPA has no calibrated floor of its own; "on" = nonzero
        if rna_on != hpa_on:
            return False
    return True


def check_geo_corroboration(rna_values, geo_values, rna_constants):
    shared = [g for g in rna_values if g in geo_values]
    if not shared:
        return None

    for gene in shared:
        entry = rna_constants["global"].get(gene)
        if entry is None:
            continue
        floor = entry["L"]
        rna_on = rna_values[gene] >= floor
        geo_on = geo_values[gene] > 0
        if rna_on != geo_on:
            return False
    return True


def confidence_tier(per_gene_layer_counts, hpa_agreement, inclusion_abstention_fraction,
                     insufficient_fraction=0.5, high_min_layers=3, high_fraction_threshold=0.8,
                     moderate_min_layers=2, moderate_fraction_threshold=0.8):
    if inclusion_abstention_fraction >= insufficient_fraction or not per_gene_layer_counts:
        return "Insufficient"

    n = len(per_gene_layer_counts)
    frac_ge_high = sum(1 for c in per_gene_layer_counts.values() if c >= high_min_layers) / n
    frac_ge_moderate = sum(1 for c in per_gene_layer_counts.values() if c >= moderate_min_layers) / n

    if frac_ge_high >= high_fraction_threshold and hpa_agreement is not False:
        return "High"
    if frac_ge_moderate >= moderate_fraction_threshold:
        return "Moderate"
    return "Low"
