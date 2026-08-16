"""Complete research-record export for S08 (C6), using only an immutable snapshot and comparison."""

import copy
import csv
import io
import json

import compare
import serialization
from backend.research_service.routes import serialization as route_serialization

EXPORT_SCHEMA_VERSION = "research-export-v1"
PARTITION_BUCKETS = (
    "ranked_cell_lines",
    "ranked_beyond_top_n",
    "low_confidence_lines",
    "insufficient_evidence_lines",
    "disqualified_lines",
)


def _attached_or_not_computed(value, label):
    if value is None:
        return {"status": "not_computed", "reason": f"{label} was not attached to this snapshot."}
    return copy.deepcopy(value)


def _method_policy(snapshot):
    """Export an attached policy verbatim, never reconstruct a policy from a partial result."""
    if isinstance(snapshot.get("method_policy"), dict):
        return copy.deepcopy(snapshot["method_policy"])
    alternatives = snapshot.get("alternatives")
    if isinstance(alternatives, dict) and isinstance(alternatives.get("policy"), dict):
        return copy.deepcopy(alternatives["policy"])
    return {"status": "not_computed", "reason": "Similarity policy was not attached to this snapshot."}


def _json_value(value):
    """Use strict JSON text for nested CSV cells while leaving JSON export values untouched."""
    if isinstance(value, (dict, list)):
        return json.dumps(serialization.to_json_safe(value), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return value


def _csv_safe(value):
    """Prefix spreadsheet-formula-looking text only in CSV, never in authoritative JSON."""
    value = _json_value(value)
    if value is None:
        return ""
    text = str(value)
    return f"'{text}" if text.startswith(("=", "+", "-", "@")) else text


def _comparison_csv(comparison):
    """Flatten ordered C6 rows without changing their values or adding a comparison score."""
    output = io.StringIO(newline="")
    model_ids = comparison["selected_model_ids"]
    fieldnames = ["row_id", "section", "label", "gene_id", "gene_role", "layer", "unit", "provenance", *model_ids]
    writer = csv.DictWriter(output, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    for row in comparison["rows"]:
        record = {
            "row_id": row["row_id"], "section": row["section"], "label": row["label"],
            "gene_id": row.get("gene_id"), "gene_role": row.get("gene_role"),
            "layer": row.get("layer"), "unit": row.get("unit"), "provenance": row.get("provenance"),
        }
        for model_id in model_ids:
            record[model_id] = _csv_safe(row["values"][model_id])
        writer.writerow(record)
    return output.getvalue()


_ROUTE_CSV_FIELDS = [
    "record_type", "research_query_id", "relationship_id", "route_policy_version",
    "route_family", "route_status", "route_qualified", "source_gene", "target_gene",
    "context_type", "ModelID", "ensembl_id", "source_id", "source_version", "data_id",
    "raw_value", "unit", "transformation", "denominator", "evidence_state", "qc_state",
    "limitation", "source_record_id", "evidence_id",
    "automatic_task_id", "automatic_task_status", "automatic_task_reason",
    "automatic_pair", "automatic_policy_identity", "automatic_counts",
    "automatic_omitted_pairs", "automatic_result", "automatic_conflicts",
]


def _automatic_or_not_computed(snapshot):
    """Expose persisted AUTO03 state exactly as saved; exports never start background work."""
    value = snapshot.get("automatic_routes")
    if value is None:
        return {
            "status": "not_computed",
            "reason": "Automatic route execution has not been started for this saved query.",
        }
    return copy.deepcopy(value)


def _route_csv(snapshot):
    """Flatten stored long ledgers only; an empty attachment writes an explicit non-fabricated row."""
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=_ROUTE_CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    attachments = route_serialization.attachments(snapshot.get("routes"))
    rows = []
    for relationship_id, entry in attachments["by_relationship_id"].items():
        route = entry.get("result", {})
        relationship = entry.get("request_orientation") or route.get("relationship") or {}
        source = relationship.get("source_gene") or {}
        target = relationship.get("target_gene") or {}
        base = {
            "record_type": "route_evidence", "research_query_id": snapshot.get("query_id"),
            "relationship_id": relationship_id, "route_policy_version": entry.get("route_policy_version"),
            "route_family": route.get("route_family"), "route_status": route.get("route_status"),
            "route_qualified": route.get("route_qualified"), "source_gene": source.get("symbol"),
            "target_gene": target.get("symbol"), "context_type": relationship.get("context_type"),
        }
        for ledger in route.get("evidence_ledger", []) or []:
            row = dict(base)
            row.update({field: ledger.get(field) for field in _ROUTE_CSV_FIELDS if field in ledger})
            rows.append(row)
    if not rows:
        rows.append({
            "record_type": "route_state", "research_query_id": snapshot.get("query_id"),
            "evidence_state": attachments["status"],
            "limitation": "No relationship-specific route was declared; no evidence rows were fabricated.",
        })
    automatic = _automatic_or_not_computed(snapshot)
    if automatic.get("status") != "not_computed":
        plan = automatic.get("plan", {})
        state = plan.get("state", {})
        rows.append({
            "record_type": "automatic_route_plan", "research_query_id": snapshot.get("query_id"),
            "route_policy_version": plan.get("policy_version"),
            "automatic_policy_identity": plan.get("policy_identity"),
            "automatic_counts": state,
            "automatic_omitted_pairs": plan.get("omitted_pairs", []),
            "automatic_conflicts": automatic.get("conflicts", []),
            "limitation": "; ".join(automatic.get("limitations", [])),
        })
        results = automatic.get("results_by_task_id", {})
        for task in plan.get("tasks", []):
            rows.append({
                "record_type": "automatic_route_task", "research_query_id": snapshot.get("query_id"),
                "route_policy_version": plan.get("policy_version"),
                "route_family": task.get("route_family"), "context_type": task.get("context_type"),
                "automatic_task_id": task.get("task_id"), "automatic_task_status": task.get("status"),
                "automatic_task_reason": task.get("reason"), "automatic_pair": task.get("canonical_gene_pair"),
                "automatic_policy_identity": plan.get("policy_identity"),
                "automatic_result": results.get(task.get("task_id")),
            })
    for row in sorted(rows, key=lambda item: (
        str(item.get("record_type") or ""), str(item.get("automatic_task_id") or ""),
        str(item.get("relationship_id") or ""), str(item.get("ModelID") or ""),
        str(item.get("ensembl_id") or ""), str(item.get("source_id") or ""), str(item.get("evidence_id") or ""),
    )):
        writer.writerow({field: _csv_safe(row.get(field)) for field in _ROUTE_CSV_FIELDS})
    return output.getvalue()


def export_routes(snapshot, format, route_registry=None, route_controls=None):
    """Authoritative BR07 route record from saved attachments, with no adapter or engine call."""
    attachments = route_serialization.attachments(snapshot.get("routes"))
    if format == "json":
        payload = {
            "schema_version": "biological-routes-export-v1",
            "snapshot_identity": {
                "query_id": snapshot.get("query_id"),
                "snapshot_schema_version": snapshot.get("schema_version"),
                "manifest": copy.deepcopy(snapshot.get("manifest")),
                "original_research_snapshot_link": f"/research/queries/{snapshot.get('query_id')}",
            },
            "original_research_query": copy.deepcopy(snapshot.get("query")),
            "routes": copy.deepcopy(attachments),
            "automatic_routes": _automatic_or_not_computed(snapshot),
            "route_registry": copy.deepcopy(route_registry) if route_registry is not None else None,
            "route_controls": copy.deepcopy(route_controls) if route_controls is not None else None,
            "native_molecular_context": {
                "native_result": copy.deepcopy(snapshot.get("native_result")),
                "context": copy.deepcopy(snapshot.get("context")) if snapshot.get("context") is not None else {
                    "status": "not_computed", "reason": "No graph/pathway context was attached to this snapshot."
                },
                "statement": "Native molecular context is copied from the immutable snapshot; no molecular evidence was recomputed.",
            },
            "limitations": copy.deepcopy(snapshot.get("limitations", [])),
        }
        return {"format": "json", "media_type": "application/json", "content": json.dumps(
            serialization.to_json_safe(payload), indent=2, sort_keys=True, allow_nan=False
        )}
    if format == "csv":
        return {"format": "csv", "media_type": "text/csv", "content": _route_csv(snapshot)}
    raise ValueError("format must be either 'json' or 'csv'")


def export_research(snapshot, comparison_ids, format):
    """Return a JSON or flat-CSV research export without invoking a scorer or data reader.

    JSON is authoritative and retains the complete native partition, snapshot, attached research
    payloads, manifest, and limitations. CSV intentionally contains only the selected comparison
    rows and their unit/state/provenance columns.
    """
    comparison = compare.compare_models(snapshot, comparison_ids)
    native_result = snapshot.get("native_result", {})
    if format == "json":
        payload = {
            "schema_version": EXPORT_SCHEMA_VERSION,
            "selected_comparison_ids": list(comparison_ids),
            "snapshot": copy.deepcopy(snapshot),
            "native_partitions": {bucket: copy.deepcopy(native_result.get(bucket, [])) for bucket in PARTITION_BUCKETS},
            "comparison": comparison,
            "research_payloads": {
                "context": _attached_or_not_computed(snapshot.get("context"), "Context"),
                "alternatives": _attached_or_not_computed(snapshot.get("alternatives"), "Alternatives"),
                "routes": _attached_or_not_computed(snapshot.get("routes"), "Biological routes"),
                "automatic_routes": _automatic_or_not_computed(snapshot),
            },
            "method_policy": _method_policy(snapshot),
            "build_identity": {"manifest": copy.deepcopy(snapshot.get("manifest")), "query_id": snapshot.get("query_id")},
            "limitations": copy.deepcopy(snapshot.get("limitations", [])) + copy.deepcopy(comparison["limitations"]),
        }
        return {"format": "json", "media_type": "application/json", "content": json.dumps(serialization.to_json_safe(payload), indent=2, allow_nan=False)}
    if format == "csv":
        return {"format": "csv", "media_type": "text/csv", "content": _comparison_csv(comparison)}
    raise ValueError("format must be either 'json' or 'csv'")
