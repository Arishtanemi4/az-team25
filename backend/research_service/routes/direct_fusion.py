"""BR03 direct-fusion evaluation over the immutable BR02 evidence projection.

The evaluator deliberately does not read processed tables, change the adapter ledger, or calculate
a technical score.  It only asks whether an exact unordered resolved pair is retained for each
model, while returning the original source fields as separately labelled evidence.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping

from .contracts import EvidenceState, ModelSupportState, Relationship, ResolvedGene, RouteFamily, RouteStatus

_RESEARCH_SERVICE_DIR = str(Path(__file__).resolve().parents[1])
if _RESEARCH_SERVICE_DIR not in sys.path:
    sys.path.insert(0, _RESEARCH_SERVICE_DIR)

from serialization import canonical_bytes, to_json_safe  # noqa: E402 -- sys.path must be set up first


DIRECT_FUSION_SCHEMA_VERSION = "br03-direct-fusion-v1"
_ASSAYED_STATES = {EvidenceState.OBSERVED.value, EvidenceState.MEASURED_ABSENT.value}
def _as_pair(value: Any) -> tuple[str, str] | None:
    """Return a canonical two-gene pair only when a source row names two non-empty IDs."""
    if not isinstance(value, Mapping):
        return None
    first = value.get("ensembl_id")
    second = value.get("partner_ensembl_id")
    if not isinstance(first, str) or not first or not isinstance(second, str) or not second:
        return None
    return tuple(sorted((first, second)))


def _technical_evidence(raw_value: Mapping[str, Any]) -> dict[str, Any]:
    """Keep technical call attributes distinct; no confidence/frame/read composite is made."""
    breakpoints = {
        key: raw_value.get(key)
        for key in sorted(raw_value)
        if "breakpoint" in key.lower()
    }
    domains = {
        key: raw_value.get(key)
        for key in sorted(raw_value)
        if "domain" in key.lower()
    }
    return {
        "TotalReadsSupportingFusion": raw_value.get("TotalReadsSupportingFusion"),
        "FFPM": raw_value.get("FFPM"),
        "confidence": raw_value.get("confidence"),
        "confidence_high": raw_value.get("confidence_high"),
        "reading_frame": raw_value.get("reading_frame"),
        "in_frame": raw_value.get("in_frame"),
        "breakpoints": breakpoints,
        "retained_domains": domains,
    }


def _relationship(adapter_result: Mapping[str, Any]) -> tuple[dict[str, Any], tuple[str, str]]:
    """Validate the BR02 direct-fusion relationship without guessing identifiers or orientation."""
    relationship = adapter_result.get("relationship")
    if not isinstance(relationship, Mapping):
        raise ValueError("resolved BR02 direct-fusion evidence requires a relationship")
    if relationship.get("route_family") != RouteFamily.DIRECT_FUSION.value:
        raise ValueError("direct fusion evaluator requires a direct_fusion adapter result")
    source = relationship.get("source_gene")
    target = relationship.get("target_gene")
    if not isinstance(source, Mapping) or not isinstance(target, Mapping):
        raise ValueError("relationship requires resolved source_gene and target_gene")
    source_id = source.get("ensembl_id")
    target_id = target.get("ensembl_id")
    if not isinstance(source_id, str) or not source_id or not isinstance(target_id, str) or not target_id:
        raise ValueError("relationship genes require non-empty canonical Ensembl IDs")
    if source_id == target_id:
        raise ValueError("direct fusion requires two distinct canonical Ensembl IDs")
    return dict(relationship), tuple(sorted((source_id, target_id)))


def _unresolvable_result(adapter_result: Mapping[str, Any]) -> dict[str, Any]:
    """Return an explicit unavailable result when BR02 could not resolve an input identity."""
    query_id = adapter_result.get("research_query_id")
    if not isinstance(query_id, str) or not query_id:
        raise ValueError("adapter result requires a non-empty research_query_id")
    return {
        "schema_version": DIRECT_FUSION_SCHEMA_VERSION,
        "research_query_id": query_id,
        "route_family": RouteFamily.DIRECT_FUSION.value,
        "route_policy_version": adapter_result.get("route_policy_version"),
        "relationship": None,
        "relationship_id": None,
        "route_status": RouteStatus.UNAVAILABLE.value,
        "route_qualified": False,
        "qualification_statistic": {
            "name": "exact_retained_pair_call_count",
            "value": None,
            "direction": "At least one exact unordered canonical pair call is required.",
            "denominator": {"snapshot_models": 0, "fusion_assayed_models": 0},
        },
        "qualification_direction": "At least one exact unordered canonical pair call is required.",
        "denominator": {"snapshot_models": 0, "fusion_assayed_models": 0},
        "model_packets": [],
        "evidence_ledger": [],
        "source": {"route_input_sources": list(adapter_result.get("route_input_sources", []))},
        "qc": {"technical_fields_are_uncombined": True},
        "missingness": {"evidence_state": EvidenceState.UNRESOLVABLE.value},
        "limitations": ["An unresolved or ambiguous frozen gene-reference entry is not guessed."],
    }


def evaluate_direct_fusion(adapter_result: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate exact direct-fusion calls from one immutable BR02 adapter result.

    ``adapter_result`` is the complete, read-only projection returned by
    :func:`adapt_route_evidence`.  Exact pair matching is orientation-independent; the supplied
    relationship remains unchanged in the response so the submitted source/target order is clear.
    """
    if not isinstance(adapter_result, Mapping):
        raise ValueError("adapter_result must be a mapping")
    original_bytes = canonical_bytes(adapter_result)
    if adapter_result.get("route_family") != RouteFamily.DIRECT_FUSION.value:
        raise ValueError("direct fusion evaluator requires a direct_fusion adapter result")
    if adapter_result.get("evidence_state") == EvidenceState.UNRESOLVABLE.value:
        result = _unresolvable_result(adapter_result)
        json.dumps(to_json_safe(result), sort_keys=True, separators=(",", ":"), allow_nan=False)
        return result

    query_id = adapter_result.get("research_query_id")
    policy_version = adapter_result.get("route_policy_version")
    if not isinstance(query_id, str) or not query_id:
        raise ValueError("adapter result requires a non-empty research_query_id")
    if not isinstance(policy_version, str) or not policy_version:
        raise ValueError("adapter result requires a non-empty route_policy_version")
    relationship, expected_pair = _relationship(adapter_result)
    relationship_id = adapter_result.get("relationship_id")
    if not isinstance(relationship_id, str) or not relationship_id:
        # BR02 retains the relationship object rather than an attachment key.  The immutable
        # contracts own the deterministic ID calculation, so use it without changing orientation.
        relationship_id = Relationship(
            RouteFamily.DIRECT_FUSION,
            ResolvedGene(str(relationship["source_gene"]["symbol"]), str(relationship["source_gene"]["ensembl_id"])),
            ResolvedGene(str(relationship["target_gene"]["symbol"]), str(relationship["target_gene"]["ensembl_id"])),
        ).relationship_id(policy_version)

    model_spine = adapter_result.get("model_spine")
    ledger = adapter_result.get("ledger")
    if not isinstance(model_spine, list) or not isinstance(ledger, list):
        raise ValueError("resolved BR02 direct-fusion evidence requires model_spine and ledger lists")
    models: dict[str, dict[str, Any]] = {}
    for item in model_spine:
        if not isinstance(item, Mapping) or not isinstance(item.get("ModelID"), str) or not item["ModelID"]:
            raise ValueError("model_spine requires non-empty ModelID values")
        models[item["ModelID"]] = dict(item)
    if len(models) != len(model_spine):
        raise ValueError("model_spine ModelID values must be unique")

    rows_by_model: dict[str, list[dict[str, Any]]] = defaultdict(list)
    exact_by_model: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in ledger:
        if not isinstance(row, Mapping):
            raise ValueError("ledger must contain mappings")
        model_id = row.get("ModelID")
        if not isinstance(model_id, str) or model_id not in models:
            raise ValueError("ledger ModelID must belong to the immutable model_spine")
        source_id = row.get("source_id")
        if source_id != "processed_fusions":
            raise ValueError("direct-fusion adapter ledger may only contain processed_fusions evidence")
        copied = dict(row)
        rows_by_model[model_id].append(copied)
        raw_value = copied.get("raw_value")
        if copied.get("evidence_state") == EvidenceState.OBSERVED.value and _as_pair(raw_value) == expected_pair:
            exact_by_model[model_id].append(copied)

    packets: list[dict[str, Any]] = []
    assayed_models = 0
    exact_event_count = 0
    for model_id in sorted(models):
        model_rows = rows_by_model[model_id]
        exact_rows = sorted(exact_by_model[model_id], key=lambda row: (str(row.get("evidence_id", "")), canonical_bytes(row)))
        is_assayed = any(row.get("evidence_state") in _ASSAYED_STATES for row in model_rows)
        if is_assayed:
            assayed_models += 1
        exact_event_count += len(exact_rows)
        if exact_rows:
            support_state = ModelSupportState.SUPPORT.value
            evidence_state = EvidenceState.OBSERVED.value
        elif is_assayed:
            support_state = ModelSupportState.CONFLICT.value
            evidence_state = EvidenceState.MEASURED_ABSENT.value
        else:
            support_state = ModelSupportState.UNKNOWN.value
            evidence_state = EvidenceState.NOT_ASSAYED.value
        matching_events = []
        for row in exact_rows:
            raw_value = row.get("raw_value")
            if not isinstance(raw_value, Mapping):
                raise ValueError("observed fusion ledger rows require a complete source-record mapping")
            matching_events.append({
                "evidence_id": row.get("evidence_id"),
                "source_record_id": row.get("source_record_id"),
                "technical_evidence": _technical_evidence(raw_value),
                "source_record": dict(raw_value),
            })
        packets.append({
            "ModelID": model_id,
            "lineage": models[model_id].get("lineage"),
            "support_state": support_state,
            "evidence_state": evidence_state,
            "fusion_assay_state": "assayed" if is_assayed else "not_assayed",
            "route_native_value": {
                "exact_retained_pair_event_count": len(matching_events),
                "matching_events": matching_events,
            },
            "limitation": (
                "All matching retained source events are shown separately; their technical fields are not combined."
                if matching_events else
                "No exact retained pair row is available for this model; not-assayed is distinct from an assayed no-call."
            ),
        })

    exact_model_count = sum(bool(events) for events in exact_by_model.values())
    if exact_event_count:
        route_status = RouteStatus.SUPPORTED.value
        route_qualified = True
    elif assayed_models:
        route_status = RouteStatus.CONTRADICTED.value
        route_qualified = False
    else:
        route_status = RouteStatus.UNAVAILABLE.value
        route_qualified = False
    result = {
        "schema_version": DIRECT_FUSION_SCHEMA_VERSION,
        "research_query_id": query_id,
        "route_family": RouteFamily.DIRECT_FUSION.value,
        "route_policy_version": policy_version,
        "relationship": relationship,
        "relationship_id": relationship_id,
        "route_status": route_status,
        "route_qualified": route_qualified,
        "qualification_statistic": {
            "name": "exact_retained_pair_call_count",
            "value": exact_event_count,
            "direction": "At least one exact unordered canonical pair call supports this direct-fusion route.",
            "denominator": {
                "snapshot_models": len(models),
                "fusion_assayed_models": assayed_models,
                "exact_pair_models": exact_model_count,
                "exact_pair_events": exact_event_count,
            },
        },
        "qualification_direction": "At least one exact unordered canonical pair call supports this direct-fusion route.",
        "denominator": {
            "snapshot_models": len(models),
            "fusion_assayed_models": assayed_models,
            "exact_pair_models": exact_model_count,
            "exact_pair_events": exact_event_count,
        },
        "model_packets": packets,
        "evidence_ledger": sorted((dict(row) for row in ledger), key=lambda row: (str(row.get("ModelID", "")), str(row.get("evidence_id", "")))),
        "source": {
            "route_input_sources": list(adapter_result.get("route_input_sources", [])),
            "adapter_schema_version": adapter_result.get("schema_version"),
            "data_id": "processed_fusions retained source rows via BR02 evidence ledger",
        },
        "qc": {
            "technical_fields_are_uncombined": True,
            "exact_pair_rule": "unordered canonical Ensembl pair equality",
            "low_confidence_and_out_of_frame_retained": True,
        },
        "missingness": {
            "fusion_assayed_models": assayed_models,
            "not_assayed_models": len(models) - assayed_models,
            "measured_no_exact_pair_models": sum(packet["support_state"] == ModelSupportState.CONFLICT.value for packet in packets),
        },
        "limitations": [
            "A retained processed fusion call is observed processed evidence, not independent biological validation.",
            "Technical fusion fields remain separate and are not a fusion probability, vote, score, causal claim, protein claim, or ranking signal.",
            "The route does not change native D, rank, tier, veto, candidate partitions, filters, or preprocessing.",
        ],
    }
    safe_result = to_json_safe(result)
    json.dumps(safe_result, sort_keys=True, separators=(",", ":"), allow_nan=False)
    if canonical_bytes(adapter_result) != original_bytes:
        raise RuntimeError("BR03 direct-fusion evaluation attempted to mutate immutable adapter evidence")
    return safe_result


# A descriptive alias keeps the route-engine API obvious to future BR07 routing code.
evaluate_direct_fusion_route = evaluate_direct_fusion
