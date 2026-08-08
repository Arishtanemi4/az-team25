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
