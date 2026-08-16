"""Read-only side-by-side research comparison for one immutable snapshot (C6).

This module assembles already-scored native evidence and any already-attached research packets. It
does not calculate a score, infer missing evidence, compare different units numerically, or decide
a winner.  A missing packet is labelled ``not_computed`` rather than built during comparison.
"""

import copy

COMPARISON_SCHEMA_VERSION = "research-comparison-v1"
PARTITION_BUCKETS = (
    "ranked_cell_lines",
    "ranked_beyond_top_n",
    "low_confidence_lines",
    "insufficient_evidence_lines",
    "disqualified_lines",
)
LAYER_ORDER = ("rna", "protein", "dependency", "copy_number", "mutation", "fusion")


def _validate_selection(model_ids):
    """C6 permits exactly two or three distinct snapshot model IDs, in user-selected order."""
    if not isinstance(model_ids, (list, tuple)) or len(model_ids) not in {2, 3}:
        raise ValueError("Select exactly two or three model IDs for comparison")
    if len(set(model_ids)) != len(model_ids):
        raise ValueError("Comparison model IDs must be distinct")


def _native_by_id(native_result):
    """Map every retained native partition record to its bucket without dropping Low/vetoed lines."""
    indexed = {}
    for bucket in PARTITION_BUCKETS:
        for line in native_result.get(bucket, []) or []:
            model_id = line.get("model_id")
            if model_id in indexed:
                raise ValueError(f"model_id {model_id!r} appears in more than one native partition")
            indexed[model_id] = {"bucket": bucket, "line": line}
    return indexed


def _metadata_by_id(snapshot):
    """Snapshot model metadata is optional context, never an excuse to omit native evidence."""
    rows = snapshot.get("models") or []
    if not isinstance(rows, list):
        return {}
    return {row.get("ModelID"): row for row in rows if isinstance(row, dict) and row.get("ModelID")}


def _context_parts(snapshot):
    """Read the extension-owned attached-packet convention without computing any missing packet."""
    context = snapshot.get("context")
    if not isinstance(context, dict):
        return {}, None
    measured = context.get("measured_by_model") or context.get("measured_contexts") or {}
    graph = context.get("graph_context") or context.get("graph")
    return measured if isinstance(measured, dict) else {}, graph if isinstance(graph, dict) else None


def _matching_measured_context(measured, snapshot, model_id):
    """Reject an attached measured packet from another snapshot rather than mixing query evidence."""
    packet = measured.get(model_id)
    if packet is None:
        return None
    if packet.get("query_id") != snapshot.get("query_id") or packet.get("model_id") != model_id:
        raise ValueError("Attached measured context has a mixed query ID or model ID")
    return packet


def _matching_alternatives(snapshot):
    """Return attached alternatives only when they belong to this snapshot; otherwise reject mixing."""
    alternatives = snapshot.get("alternatives")
    if alternatives is None:
        return None
    if not isinstance(alternatives, dict) or alternatives.get("query_id") != snapshot.get("query_id"):
        raise ValueError("Attached alternatives have a mixed query ID")
    return alternatives


def _row(row_id, section, label, model_ids, values, **metadata):
    """Build one stable comparison row with values kept in the caller's selection order."""
    return {
        "row_id": row_id,
        "section": section,
        "label": label,
        "values": {model_id: copy.deepcopy(values.get(model_id)) for model_id in model_ids},
        **metadata,
    }


def _gene_order_and_roles(snapshot):
    """Keep inclusion then exclusion request order, de-duplicating only repeated canonical IDs."""
    query = snapshot.get("query") or snapshot.get("native_result", {}).get("query") or {}
    ordered, roles = [], {}
    for role, field in (("inclusion", "inclusion_genes"), ("exclusion", "exclusion_genes")):
        for gene_id in query.get(field, []) or []:
            if gene_id not in roles:
                ordered.append(gene_id)
                roles[gene_id] = role
    return ordered, roles


def _per_gene(line, gene_id):
    for gene in line.get("per_gene", []) or []:
        if gene.get("ensembl_id") == gene_id:
            return gene
    return None


def _packet_gene(packet, gene_id):
    if packet is None:
        return None
    for gene in packet.get("genes", []) or []:
        if gene.get("ensembl_id") == gene_id:
            return gene
    return None


def _event_value(packet_gene, layer_name):
    """S04 event rows are additive context. No packet means unknown/not-computed, not no events."""
    if packet_gene is None:
        return {"status": "not_computed", "reason": "Measured context was not attached to this snapshot."}
    layer = packet_gene.get("layers", {}).get(layer_name)
    if layer is None:
        return {"status": "not_available", "reason": "Layer is absent from the attached measured context."}
    return copy.deepcopy(layer.get("events", {"status": "not_computed", "reason": "No event packet was attached."}))


def _pathway_value(graph, gene_id):
    """Graph membership is generic reference context and is never made cell-line-specific here."""
    if graph is None:
        return {"status": "not_computed", "reason": "Graph context was not attached to this snapshot."}
    gene = graph.get("genes", {}).get(gene_id)
    if gene is None:
        return {"status": "unavailable", "reason": "Graph context has no entry for this gene."}
    return {
        "status": "ok",
        "graph_source": copy.deepcopy(graph.get("graph_source")),
        "membership": copy.deepcopy(gene.get("pathway_membership")),
        "claim_type": gene.get("claim_type"),
        "reason": gene.get("reason"),
    }


def compare_models(snapshot, model_ids):
    """Create an ordered C6 comparison from a single snapshot without rescoring.

    Values remain nested by selected ModelID so the JSON comparison preserves types exactly. The
    companion CSV export flattens these rows but does not perform cross-unit subtraction, ranking,
    colouring, or any other winner inference.
    """
    _validate_selection(model_ids)
    native = _native_by_id(snapshot.get("native_result", {}))
    unknown = [model_id for model_id in model_ids if model_id not in native]
    if unknown:
        raise ValueError(f"Unknown or mixed-snapshot comparison model IDs: {unknown}")
    metadata = _metadata_by_id(snapshot)
    measured_by_model, graph = _context_parts(snapshot)
    measured = {model_id: _matching_measured_context(measured_by_model, snapshot, model_id) for model_id in model_ids}
    alternatives = _matching_alternatives(snapshot)
    genes, roles = _gene_order_and_roles(snapshot)
    lines = {model_id: native[model_id]["line"] for model_id in model_ids}
    packets = {model_id: measured[model_id] for model_id in model_ids}

    rows = []
    rows.append(_row("model_id", "model", "Model ID", model_ids, {model_id: model_id for model_id in model_ids}, provenance="snapshot.native_result"))
    for field, label in (
        ("cell_line_name", "Cell-line name"), ("lineage", "Lineage"),
        ("primary_disease", "Primary disease"), ("is_problematic", "Problematic flag"),
    ):
        rows.append(_row(field, "model", label, model_ids, {model_id: metadata.get(model_id, {}).get(field) for model_id in model_ids}, provenance="snapshot.models"))
    rows.append(_row("partition_bucket", "native", "Native partition", model_ids, {model_id: native[model_id]["bucket"] for model_id in model_ids}, provenance="snapshot.native_result"))
    for field, label in (("D", "Native desirability (D)"), ("confidence_tier", "Native confidence tier"), ("veto", "Veto"), ("warnings", "Warnings"), ("hpa_agreement", "HPA agreement"), ("geo_agreement", "GEO agreement")):
        rows.append(_row(field, "native", label, model_ids, {model_id: lines[model_id].get(field) for model_id in model_ids}, provenance="snapshot.native_result"))
    rows.append(_row("measured_context_status", "availability", "Measured context", model_ids, {model_id: "computed" if packets[model_id] else "not_computed" for model_id in model_ids}, provenance="extensions.context"))
    rows.append(_row("graph_context_status", "availability", "Graph context", model_ids, {model_id: "computed" if graph is not None else "not_computed" for model_id in model_ids}, provenance="extensions.graph_context"))
    rows.append(_row("alternatives_status", "availability", "Alternatives", model_ids, {model_id: alternatives.get("status") if alternatives and alternatives.get("anchor_id") == model_id else "not_computed" for model_id in model_ids}, provenance="extensions.alternatives"))

    for gene_id in genes:
        role = roles[gene_id]
        rows.append(_row(f"{gene_id}:role", "gene", "Gene role", model_ids, {model_id: role for model_id in model_ids}, gene_id=gene_id, gene_role=role, provenance="snapshot.query"))
        rows.append(_row(f"{gene_id}:symbol", "gene", "Gene symbol", model_ids, {model_id: (_per_gene(lines[model_id], gene_id) or {}).get("symbol") for model_id in model_ids}, gene_id=gene_id, gene_role=role, provenance="snapshot.native_result"))
        rows.append(_row(f"{gene_id}:d_gene", "gene", "Per-gene desirability", model_ids, {model_id: (_per_gene(lines[model_id], gene_id) or {}).get("d_gene") for model_id in model_ids}, gene_id=gene_id, gene_role=role, provenance="snapshot.native_result"))
        rows.append(_row(f"{gene_id}:missing_evidence", "gene", "Missing evidence", model_ids, {model_id: [copy.deepcopy(item) for item in lines[model_id].get("missing_evidence", []) if item.get("ensembl_id") == gene_id] for model_id in model_ids}, gene_id=gene_id, gene_role=role, provenance="snapshot.native_result"))
        for layer_name in LAYER_ORDER:
            layer_values = {model_id: ((_per_gene(lines[model_id], gene_id) or {}).get("layers", {}).get(layer_name)) for model_id in model_ids}
            rows.append(_row(f"{gene_id}:{layer_name}:value", "layer", f"{layer_name} value", model_ids, {model_id: (layer_values[model_id] or {}).get("value") for model_id in model_ids}, gene_id=gene_id, gene_role=role, layer=layer_name, provenance="snapshot.native_result"))
            rows.append(_row(f"{gene_id}:{layer_name}:unit", "layer", f"{layer_name} unit", model_ids, {model_id: (layer_values[model_id] or {}).get("unit") for model_id in model_ids}, gene_id=gene_id, gene_role=role, layer=layer_name, provenance="snapshot.native_result"))
            rows.append(_row(f"{gene_id}:{layer_name}:d", "layer", f"{layer_name} desirability contribution", model_ids, {model_id: (layer_values[model_id] or {}).get("d") for model_id in model_ids}, gene_id=gene_id, gene_role=role, layer=layer_name, provenance="snapshot.native_result"))
            rows.append(_row(f"{gene_id}:{layer_name}:state", "layer", f"{layer_name} state", model_ids, {model_id: (layer_values[model_id] or {}).get("state") for model_id in model_ids}, gene_id=gene_id, gene_role=role, layer=layer_name, provenance="snapshot.native_result"))
        for layer_name in ("mutation", "fusion"):
            rows.append(_row(f"{gene_id}:{layer_name}:events", "events", f"Observed {layer_name} events", model_ids, {model_id: _event_value(_packet_gene(packets[model_id], gene_id), layer_name) for model_id in model_ids}, gene_id=gene_id, gene_role=role, layer=layer_name, provenance="extensions.context"))
        rows.append(_row(f"{gene_id}:pathways", "graph", "Source-labelled pathway membership", model_ids, {model_id: _pathway_value(graph, gene_id) for model_id in model_ids}, gene_id=gene_id, gene_role=role, provenance="extensions.graph_context"))

    limitations = []
    if graph is None:
        limitations.append("Graph context not_computed; no pathway membership was fabricated.")
    if any(packet is None for packet in packets.values()):
        limitations.append("Measured context not_computed for at least one selected model; event rows are not inferred from absence.")
    if alternatives is None:
        limitations.append("Alternatives not_computed; native legacy similar_lines remain in the snapshot.")
    return {
        "schema_version": COMPARISON_SCHEMA_VERSION,
        "query_id": snapshot.get("query_id"),
        "selected_model_ids": list(model_ids),
        "rows": rows,
        "limitations": limitations,
        "boundary": "Rows preserve source units and states. No cross-unit subtraction, colour ranking, winner badge, or causal claim is computed.",
    }
