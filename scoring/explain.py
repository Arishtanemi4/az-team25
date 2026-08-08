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


GENE_ROLE_DEPENDENT_LAYERS = ("mutation", "copy_number")


def _gene_role_abstained_layers(gene_result):
    abstained = []
    for layer in GENE_ROLE_DEPENDENT_LAYERS:
        info = gene_result["layers"].get(layer, {})
        if info.get("d") is None and info.get("state") != "not_assayed" and info.get("gene_class") == "unknown":
            abstained.append(layer)
    return abstained


IDENTITY_RESOLUTION_LAYERS = ("protein", "dependency", "copy_number")


def _identity_unresolvable_layers(gene_result, unresolved_symbols):
    if not unresolved_symbols:
        return []
    unresolvable = []
    for layer in IDENTITY_RESOLUTION_LAYERS:
        info = gene_result["layers"].get(layer, {})
        if info.get("state") == "not_assayed" and gene_result["ensembl_id"] in unresolved_symbols.get(layer, ()):
            unresolvable.append(layer)
    return unresolvable


FUSION_STATUS_UNKNOWN = (
    "Fusion status unknown -- this line was never RNA-assayed for fusion calling, so an "
    "absent fusion call here is missing evidence, not a confirmed negative."
)
PAN_ESSENTIAL_NOTE = (
    "Dependency evidence withheld: this gene is broadly essential across measured cell lines "
    "(strongly negative in {fraction_strong:.0%} of profiled lines), so a strong Chronos "
    "dependency here is not selective evidence for this query. Never a score bonus, never "
    "silent."
)
GENE_ROLE_UNKNOWN_NOTE = (
    "{layer} evidence withheld: whether this reading is desirable depends on whether the gene "
    "is an oncogene or a tumour suppressor, and that role is not yet classified for this gene "
    "(most genes in this build now have a classified role; this one does not). Abstained "
    "rather than guessed a direction."
)
PROTEIN_NON_DETECTED_NOTE = (
    "Protein assay ran but the target was not detected -- abstained (excluded from the "
    "weighted mean) rather than scored at the imputed detection floor. Biologically "
    "consistent with low expression (protein missingness tracks low RNA 6.1x more often than "
    "high) -- that context is disclosed here, never used to infer or score a specific value."
)
PROTEIN_MIXED_SCALE_NOTE = (
    "Protein evidence is a z-score, panel-relative by construction (it arrives that way and "
    "cannot be undone) -- its floor/target are in z-units, a mixed-scale limitation on this "
    "gene's evidence combination, disclosed here rather than corrected (edge case 14)."
)


def gene_narrative(gene_result, mutation_ctx=None, fusion_ctx=None, pan_essential_ctx=None):
    symbol = gene_result["symbol"]
    role = gene_result["role"]

    dependency_state = gene_result["layers"]["dependency"]["state"]
    pan_essential_clause = ""
    if dependency_state == "excluded_pan_essential" and pan_essential_ctx is not None:
        fraction = pan_essential_ctx.get("fraction_strong") or 0.0
        pan_essential_clause = f" {PAN_ESSENTIAL_NOTE.format(fraction_strong=fraction)}"

    gene_role_layers = _gene_role_abstained_layers(gene_result)
    gene_role_clause = "".join(
        f" {GENE_ROLE_UNKNOWN_NOTE.format(layer=layer.replace('_', ' ').capitalize())}"
        for layer in gene_role_layers
    )

    protein_non_detected_clause = (
        f" {PROTEIN_NON_DETECTED_NOTE}"
        if gene_result["layers"]["protein"]["state"] == "non_detected"
        else ""
    )

    if gene_result["d_gene"] is None:
        narrative = f"{symbol} ({role}): no evidence in any measured layer -- abstained, not scored."
        if gene_result["layers"]["fusion"]["state"] == "not_assayed":
            narrative += f" {FUSION_STATUS_UNKNOWN}"
        narrative += pan_essential_clause + gene_role_clause + protein_non_detected_clause
        return narrative

    layer_parts = [
        f"{layer} (d={info['d']:.2f})" for layer, info in gene_result["layers"].items()
        if info["d"] is not None
    ]
    narrative = f"{symbol} ({role}): d_gene={gene_result['d_gene']:.2f}, from " + ", ".join(layer_parts)
    if gene_result["missing_layers"]:
        narrative += f". Missing: {', '.join(gene_result['missing_layers'])}."
    if gene_result["layers"]["fusion"]["state"] == "not_assayed":
        narrative += f" {FUSION_STATUS_UNKNOWN}"
    narrative += pan_essential_clause + gene_role_clause + protein_non_detected_clause
    if "protein" in [l for l, i in gene_result["layers"].items() if i["d"] is not None]:
        narrative += f" {PROTEIN_MIXED_SCALE_NOTE}"
    if mutation_ctx is not None:
        narrative += f" {_mutation_narrative_clause(mutation_ctx)}"
    if fusion_ctx is not None:
        narrative += f" {_fusion_narrative_clause(fusion_ctx)}"
    return narrative


def missing_evidence_report(per_gene_results, unresolved_symbols=None):
    report = []
    for gene in per_gene_results:
        for layer in _gene_role_abstained_layers(gene):
            report.append({
                "ensembl_id": gene["ensembl_id"],
                "symbol": gene["symbol"],
                "role": gene["role"],
                "severity": "gene_role_unknown",
                "layer": layer,
            })
        for layer in _identity_unresolvable_layers(gene, unresolved_symbols):
            report.append({
                "ensembl_id": gene["ensembl_id"],
                "symbol": gene["symbol"],
                "role": gene["role"],
                "severity": "measured_unresolvable",
                "layer": layer,
            })
    for gene in per_gene_results:
        if gene["d_gene"] is None:
            severity = "cannot_certify_absence" if gene["role"] == "exclusion" else "untested"
            report.append({
                "ensembl_id": gene["ensembl_id"],
                "symbol": gene["symbol"],
                "role": gene["role"],
                "severity": severity,
            })
    return report


TUMOUR_REPRESENTATIVENESS_GAP = (
    "Not available: Jin et al. (2023) tumour representativeness has no corresponding table "
    "in data/processed/ or data/augmented/ (the file is reserved for external validation)."
)
HALLMARK_TAGS_GAP = (
    "Not available: no gene-to-hallmark mapping resource exists anywhere in this repo."
)
RNA_ONLY_CORRELATION_NOTE = (
    "The correlation discount (rho_bar/m_eff) is estimated from RNA co-expression alone, but "
    "is applied to d_gene, which aggregates all six evidence layers. Two genes can be "
    "RNA-co-regulated yet functionally independent, or RNA-independent yet functionally "
    "redundant -- a disclosed structural mismatch, not corrected here."
)


def build_result(model_id, D, veto_info, per_gene_results, rho_bar, m_eff, m,
                  cell_lines_row, hpa_agreement=None, geo_agreement=None,
                  fusion_context=None, mutation_context=None, pan_essential_context=None,
                  tier_params=None, build_narrative=True, unresolved_symbols=None):
    per_gene_layer_counts = {
        g["ensembl_id"]: sum(1 for info in g["layers"].values() if info["d"] is not None)
        for g in per_gene_results if g["d_gene"] is not None
    }
    inclusion_genes = [g for g in per_gene_results if g["role"] == "inclusion"]
    inclusion_abstained = sum(1 for g in inclusion_genes if g["d_gene"] is None)
    inclusion_abstention_fraction = (
        inclusion_abstained / len(inclusion_genes) if inclusion_genes else 1.0
    )

    tier = confidence_tier(
        per_gene_layer_counts, hpa_agreement, inclusion_abstention_fraction,
        **(tier_params or {}),
    )

    veto = None
    if veto_info is not None:
        vetoed_gene = next(g for g in per_gene_results if g["ensembl_id"] == veto_info["ensembl_id"])
        reason = "exclusion_expressed" if vetoed_gene["role"] == "exclusion" else "inclusion_below_floor"
        veto = {"ensembl_id": vetoed_gene["ensembl_id"], "symbol": vetoed_gene["symbol"], "reason": reason}

    mutation_by_gene = {c["ensembl_id"]: c for c in (mutation_context or [])}
    fusion_by_gene = {c["ensembl_id"]: c for c in (fusion_context or [])}
    pan_essential_by_gene = {c["ensembl_id"]: c for c in (pan_essential_context or [])}

    warnings = [
        {
            "type": "pan_essential",
            "ensembl_id": c["ensembl_id"],
            "symbol": c["symbol"],
            "message": PAN_ESSENTIAL_NOTE.format(fraction_strong=c.get("fraction_strong") or 0.0),
        }
        for c in (pan_essential_context or [])
    ]
    if m > 1:
        warnings.append({"type": "correlation_layer_scope", "message": RNA_ONLY_CORRELATION_NOTE})

    return {
        "model_id": model_id,
        "D": D,
        "confidence_tier": tier,
        "hpa_agreement": hpa_agreement,
        "geo_agreement": geo_agreement,
        "veto": veto,
        "rho_bar": rho_bar,
        "m_eff": m_eff,
        "m": m,
        "warnings": warnings,
        "per_gene": [
            {
                **g,
                "narrative": gene_narrative(
                    g, mutation_by_gene.get(g["ensembl_id"]), fusion_by_gene.get(g["ensembl_id"]),
                    pan_essential_by_gene.get(g["ensembl_id"]),
                ) if build_narrative else None,
            }
            for g in per_gene_results
        ],
        "missing_evidence": missing_evidence_report(per_gene_results, unresolved_symbols),
        "is_problematic": bool(cell_lines_row.get("is_problematic", False)),
        "tumour_representativeness": None,
        "tumour_representativeness_note": TUMOUR_REPRESENTATIVENESS_GAP,
        "hallmark_tags": None,
        "hallmark_tags_note": HALLMARK_TAGS_GAP,
        "fusion_context": [
            {"ensembl_id": c["ensembl_id"], "symbol": c["symbol"], "events": c["events"]}
            for c in (fusion_context or [])
        ],
        "similar_lines": [],
    }
