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


def _mutation_narrative_clause(mutation_ctx):
    row = mutation_ctx["rows"].iloc[0]
    parts = [f"variant {row.get('protein_change')}" if row.get("protein_change") else "a variant"]
    if row.get("vep_impact"):
        parts.append(f"{row['vep_impact']} impact")
    clin_sig = row.get("VepClinSig")
    if isinstance(clin_sig, str) and clin_sig:
        parts.append(f"ClinVar: {clin_sig}")
    civic = row.get("CivicDescription")
    if isinstance(civic, str) and civic:
        parts.append(f"CIViC: {civic[:120]}")
    return f"Somatic variant evidence: {', '.join(parts)}."


def _fusion_narrative_clause(fusion_ctx):
    rows = fusion_ctx["rows"]
    partners = sorted({p for p in rows.get("partner_ensembl_id", []) if isinstance(p, str)})
    n_events = fusion_ctx["events"]
    n_high_conf = int(rows["confidence_high"].sum())
    n_in_frame = int(rows["in_frame"].sum())
    partner_note = f" with partner(s) {', '.join(partners)}" if partners else ""
    return (
        f"Fusion evidence: {n_events} event(s){partner_note}, "
        f"{n_high_conf} high-confidence, {n_in_frame} in-frame."
    )
