"""BR04/BR05 conditional-dependency evaluation over immutable BR02 evidence.

This module deliberately reads no processed tables.  It classifies the rows already retained by
the BR02 adapter, then reports a lineage-adjusted descriptive comparison of target Chronos values.
It is not a score, a ranking input, or a claim of mechanism or model suitability.
"""

from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .contracts import ContextType, EvidenceState, ModelSupportState, Relationship, ResolvedGene, RouteFamily, RouteStatus

_RESEARCH_SERVICE_DIR = str(Path(__file__).resolve().parents[1])
if _RESEARCH_SERVICE_DIR not in sys.path:
    sys.path.insert(0, _RESEARCH_SERVICE_DIR)

from serialization import canonical_bytes, to_json_safe  # noqa: E402 -- sys.path must be set up first


CONDITIONAL_DEPENDENCY_SCHEMA_VERSION = "br05-conditional-dependency-v2"
_ASSAYED_STATES = {EvidenceState.OBSERVED.value, EvidenceState.MEASURED_ABSENT.value}
_MUTATION_TRUE_FIELDS = (
    "LikelyLoF",
    "is_driver",
    "is_hotspot",
    "OncogeneHighImpact",
    "TumorSuppressorHighImpact",
    "TranscriptLikelyLof",
)
_TRUE_WORDS = {"true", "t", "yes", "y", "1"}
_FALSE_WORDS = {"false", "f", "no", "n", "0", "", "none", "null", "nan", "na", "n/a"}
_REFERENCE_LABEL = "No qualifying event recorded in a mutation-assayed model; this is not definitive wild type without locus-level callability."


def _explicit_true(value: Any) -> bool:
    """Accept only explicit true encodings; absent or unfamiliar annotations never become true."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        return int(value) == 1
    if isinstance(value, (float, np.floating)):
        return math.isfinite(float(value)) and float(value) == 1.0
    if isinstance(value, str):
        return value.strip().lower() in _TRUE_WORDS
    return False


def _qualifying_mutation(value: Any) -> bool:
    """Apply the frozen damaging/driver-like annotation rule to one retained raw event row."""
    if not isinstance(value, Mapping):
        return False
    if any(_explicit_true(value.get(field)) for field in _MUTATION_TRUE_FIELDS):
        return True
    impact = value.get("vep_impact")
    return isinstance(impact, str) and impact.strip().upper() == "HIGH"


def _finite_number(value: Any) -> float | None:
    """Return an observed finite Chronos value only; booleans and NaN are not measurements."""
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if math.isfinite(numeric) else None


def _known_problematic(value: Any) -> bool | None:
    """Parse the model QC flag without treating an absent/ambiguous flag as a pass."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool) and int(value) in {0, 1}:
        return bool(int(value))
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in _TRUE_WORDS:
            return True
        if normalized in _FALSE_WORDS - {"", "none", "null", "nan", "na", "n/a"}:
            return False
    return None


def _relationship(adapter_result: Mapping[str, Any], context_type: ContextType) -> dict[str, Any]:
    """Validate one resolved directional context without re-resolving frozen genes."""
    relationship = adapter_result.get("relationship")
    if not isinstance(relationship, Mapping):
        raise ValueError("resolved BR02 conditional evidence requires a relationship")
    if relationship.get("route_family") != RouteFamily.CONDITIONAL_DEPENDENCY.value:
        raise ValueError("conditional evaluator requires a conditional_dependency adapter result")
    if relationship.get("context_type") != context_type.value:
        raise ValueError(f"conditional evaluator requires {context_type.value} context")
    source, target = relationship.get("source_gene"), relationship.get("target_gene")
    if not isinstance(source, Mapping) or not isinstance(target, Mapping):
        raise ValueError("relationship requires resolved source_gene and target_gene")
    for gene in (source, target):
        if not isinstance(gene.get("symbol"), str) or not gene["symbol"] or not isinstance(gene.get("ensembl_id"), str) or not gene["ensembl_id"]:
            raise ValueError("relationship genes require non-empty canonical IDs")
    if source["ensembl_id"] == target["ensembl_id"]:
        raise ValueError("conditional dependency requires distinct canonical genes")
    return dict(relationship)


def _unresolvable_result(adapter_result: Mapping[str, Any]) -> dict[str, Any]:
    """Keep BR02 identity failure explicit rather than inventing an evaluable mutation cohort."""
    query_id = adapter_result.get("research_query_id")
    if not isinstance(query_id, str) or not query_id:
        raise ValueError("adapter result requires a non-empty research_query_id")
    return {
        "schema_version": CONDITIONAL_DEPENDENCY_SCHEMA_VERSION,
        "research_query_id": query_id,
        "route_family": RouteFamily.CONDITIONAL_DEPENDENCY.value,
        "route_policy_version": adapter_result.get("route_policy_version"),
        "relationship": None,
        "relationship_id": None,
        "route_status": RouteStatus.UNAVAILABLE.value,
        "route_qualified": False,
        "qualification_statistic": {"name": "lineage_adjusted_mutation_context_chronos_shift", "value": None,
                                      "unit": "Chronos dependency score", "direction": "Positive means context more target-dependent."},
        "denominator": {"snapshot_models": 0, "context_models": 0, "reference_models": 0, "evaluable_models": 0, "eligible_lineages": 0},
        "model_packets": [], "evidence_ledger": [],
        "source": {"route_input_sources": list(adapter_result.get("route_input_sources", []))},
        "qc": {"screen_specific_dependency_qc": "unavailable in frozen BR02 output; no QC pass is inferred."},
        "limitations": ["An unresolved or ambiguous frozen gene-reference entry is not guessed."],
    }


def _stratified_bootstrap(rows: list[dict[str, Any]], replicates: int, seed: int) -> list[float]:
    """Resample adjusted Chronos values within fixed lineage-by-group strata."""
    strata: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        strata[(row["lineage"], row["group"])].append(row["adjusted_value"])
    ordered = [(key, np.asarray(sorted(values), dtype=float)) for key, values in sorted(strata.items())]
    generator = np.random.default_rng(seed)
    shifts: list[float] = []
    for _ in range(replicates):
        sampled: dict[str, list[float]] = {"context": [], "reference": []}
        for (_, group), values in ordered:
            sampled[group].extend(generator.choice(values, size=len(values), replace=True).tolist())
        shifts.append(float(np.median(sampled["reference"]) - np.median(sampled["context"])))
    return shifts


def _statistic(rows: list[dict[str, Any]]) -> float | None:
    """Calculate the locked adjusted median difference only when both groups are present."""
    context = [row["adjusted_value"] for row in rows if row["group"] == "context"]
    reference = [row["adjusted_value"] for row in rows if row["group"] == "reference"]
    if not context or not reference:
        return None
    return float(np.median(reference) - np.median(context))


def _evaluate_mutation(adapter_result: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate BR04 mutation context using only a resolved, immutable BR02 adapter payload.

    The statistic is ``median(reference adjusted Chronos) - median(context adjusted Chronos)``.
    A positive shift is therefore more negative target dependency in the mutation-context group.
    """
    if not isinstance(adapter_result, Mapping):
        raise ValueError("adapter_result must be a mapping")
    original_bytes = canonical_bytes(adapter_result)
    if adapter_result.get("route_family") != RouteFamily.CONDITIONAL_DEPENDENCY.value:
        raise ValueError("conditional evaluator requires a conditional_dependency adapter result")
    if adapter_result.get("evidence_state") == EvidenceState.UNRESOLVABLE.value:
        return _unresolvable_result(adapter_result)

    query_id, policy_version = adapter_result.get("research_query_id"), adapter_result.get("route_policy_version")
    if not isinstance(query_id, str) or not query_id or not isinstance(policy_version, str) or not policy_version:
        raise ValueError("resolved BR02 result requires query and policy versions")
    relationship = _relationship(adapter_result, ContextType.MUTATION)
    source_id, target_id = relationship["source_gene"]["ensembl_id"], relationship["target_gene"]["ensembl_id"]
    relationship_id = adapter_result.get("relationship_id")
    if not isinstance(relationship_id, str) or not relationship_id:
        relationship_id = Relationship(
            RouteFamily.CONDITIONAL_DEPENDENCY,
            ResolvedGene(relationship["source_gene"]["symbol"], source_id),
            ResolvedGene(relationship["target_gene"]["symbol"], target_id), ContextType.MUTATION,
        ).relationship_id(policy_version)

    spine, ledger = adapter_result.get("model_spine"), adapter_result.get("ledger")
    if not isinstance(spine, list) or not isinstance(ledger, list):
        raise ValueError("resolved BR02 conditional evidence requires model_spine and ledger lists")
    models: dict[str, dict[str, Any]] = {}
    for row in spine:
        if not isinstance(row, Mapping) or not isinstance(row.get("ModelID"), str) or not row["ModelID"]:
            raise ValueError("model_spine requires non-empty unique ModelID values")
        models[row["ModelID"]] = dict(row)
    if len(models) != len(spine):
        raise ValueError("model_spine ModelID values must be unique")

    mutation_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    dependency_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in ledger:
        if not isinstance(row, Mapping) or row.get("ModelID") not in models:
            raise ValueError("ledger ModelID must belong to immutable model_spine")
        copied = dict(row)
        if copied.get("source_id") == "processed_mutations" and copied.get("ensembl_id") == source_id:
            mutation_rows[copied["ModelID"]].append(copied)
        elif copied.get("source_id") == "processed_dependency" and copied.get("ensembl_id") == target_id:
            dependency_rows[copied["ModelID"]].append(copied)

    preliminary: list[dict[str, Any]] = []
    packets_by_id: dict[str, dict[str, Any]] = {}
    for model_id in sorted(models):
        model = models[model_id]
        mutations = sorted(mutation_rows[model_id], key=lambda row: (str(row.get("evidence_id", "")), canonical_bytes(row)))
        dependencies = sorted(dependency_rows[model_id], key=lambda row: (str(row.get("evidence_id", "")), canonical_bytes(row)))
        mutation_assayed = any(row.get("evidence_state") in _ASSAYED_STATES for row in mutations)
        positive = any(row.get("evidence_state") == EvidenceState.OBSERVED.value and _qualifying_mutation(row.get("raw_value")) for row in mutations)
        group = "context" if positive else "reference" if mutation_assayed else None
        observed_dependency = [row for row in dependencies if row.get("evidence_state") == EvidenceState.OBSERVED.value and _finite_number(row.get("raw_value")) is not None]
        value = _finite_number(observed_dependency[0]["raw_value"]) if len(observed_dependency) == 1 else None
        problematic = _known_problematic(model.get("is_problematic"))
        lineage = model.get("lineage") if isinstance(model.get("lineage"), str) and model.get("lineage") else None
        evaluable = group is not None and value is not None and problematic is False and lineage is not None
        preliminary.append({"ModelID": model_id, "lineage": lineage, "group": group, "raw_value": value, "evaluable": evaluable})
        packets_by_id[model_id] = {
            "ModelID": model_id, "lineage": lineage, "support_state": ModelSupportState.UNKNOWN.value,
            "evidence_state": EvidenceState.NOT_ASSAYED.value if not mutation_assayed else EvidenceState.OBSERVED.value,
            "mutation_context": "qualifying_event_recorded" if positive else "no_qualifying_event_recorded" if mutation_assayed else "unknown_not_assayed",
            "mutation_reference_label": _REFERENCE_LABEL if group == "reference" else None,
            "mutation_assay_state": "assayed" if mutation_assayed else "not_assayed",
            "dependency_assay_state": "measured" if value is not None else "unresolved_or_not_measured",
            "problematic_model_state": "problematic" if problematic is True else "not_problematic" if problematic is False else "unknown",
            "qc_state": "screen_specific_dependency_qc_unavailable_in_frozen_output",
            "route_native_value": {"raw_dependency_score": value, "raw_unit": "Chronos dependency score", "lineage_adjusted_dependency_score": None,
                                   "adjusted_unit": "Chronos dependency score", "is_probability": False},
            "evidence_ids": {"mutation": [row.get("evidence_id") for row in mutations], "dependency": [row.get("evidence_id") for row in dependencies]},
            "limitation": "Screen-specific dependency QC is unavailable in the frozen adapter output; no QC pass is inferred.",
        }

    by_lineage: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in preliminary:
        if row["evaluable"]:
            by_lineage[row["lineage"]].append(row)
    eligible_lineages = sorted(lineage for lineage, rows in by_lineage.items() if len(rows) >= 15)
    eligible_set = set(eligible_lineages)
    analysis_rows: list[dict[str, Any]] = []
    for lineage in eligible_lineages:
        rows = by_lineage[lineage]
        median = float(np.median([row["raw_value"] for row in rows]))
        for row in rows:
            row = dict(row)
            row["adjusted_value"] = float(row["raw_value"] - median)
            analysis_rows.append(row)
            packets_by_id[row["ModelID"]]["route_native_value"]["lineage_adjusted_dependency_score"] = row["adjusted_value"]

    context_rows = [row for row in analysis_rows if row["group"] == "context"]
    reference_rows = [row for row in analysis_rows if row["group"] == "reference"]
    raw_context = [row["raw_value"] for row in context_rows]
    raw_reference = [row["raw_value"] for row in reference_rows]
    adjusted_context = [row["adjusted_value"] for row in context_rows]
    adjusted_reference = [row["adjusted_value"] for row in reference_rows]
    shift = _statistic(analysis_rows)
    gates = {
        "context_minimum_n": len(context_rows) >= 10,
        "reference_minimum_n": len(reference_rows) >= 10,
        "total_minimum_n": len(analysis_rows) >= 30,
        "eligible_lineages_minimum": len(eligible_lineages) >= 2,
    }
    estimable = all(gates.values()) and shift is not None
    bootstrap = _stratified_bootstrap(analysis_rows, 500, 25) if estimable else []
    ci = [float(np.quantile(bootstrap, 0.025)), float(np.quantile(bootstrap, 0.975))] if bootstrap else [None, None]
    loo: list[dict[str, Any]] = []
    full_sign = 1 if shift is not None and shift > 0 else -1 if shift is not None and shift < 0 else 0
    for held_out in eligible_lineages:
        retained = [row for row in analysis_rows if row["lineage"] != held_out]
        retained_lineages = sorted({row["lineage"] for row in retained})
        retained_context = sum(row["group"] == "context" for row in retained)
        retained_reference = sum(row["group"] == "reference" for row in retained)
        holdout_value = _statistic(retained)
        holdout_estimable = (retained_context >= 10 and retained_reference >= 10 and len(retained) >= 30 and len(retained_lineages) >= 2 and holdout_value is not None)
        holdout_sign = 1 if holdout_value is not None and holdout_value > 0 else -1 if holdout_value is not None and holdout_value < 0 else 0
        loo.append({"held_out_lineage": held_out, "shift": holdout_value if holdout_estimable else None,
                    "estimable": holdout_estimable, "matches_full_sign": bool(holdout_estimable and full_sign != 0 and holdout_sign == full_sign),
                    "denominator": {"context_models": retained_context, "reference_models": retained_reference, "evaluable_models": len(retained), "eligible_lineages": len(retained_lineages)}})
    estimable_loo = [row for row in loo if row["estimable"]]
    consistency = (sum(row["matches_full_sign"] for row in estimable_loo) / len(estimable_loo)) if estimable_loo else None

    essential = (adapter_result.get("common_essential") or {}).get(target_id)
    common_essential = essential.get("is_common_essential") if isinstance(essential, Mapping) and isinstance(essential.get("is_common_essential"), bool) else None
    if not estimable:
        status = RouteStatus.UNAVAILABLE.value
    elif common_essential is True or common_essential is None:
        status = RouteStatus.EXPLORATORY.value
    elif shift > 0 and ci[0] is not None and ci[0] > 0 and consistency is not None and consistency >= 0.80:
        status = RouteStatus.SUPPORTED.value
    elif shift < 0 and ci[1] is not None and ci[1] < 0 and consistency is not None and consistency >= 0.80:
        status = RouteStatus.CONTRADICTED.value
    else:
        status = RouteStatus.EXPLORATORY.value

    for row in preliminary:
        packet = packets_by_id[row["ModelID"]]
        if packet["problematic_model_state"] != "not_problematic" or row["group"] is None or row["raw_value"] is None:
            packet["support_state"] = ModelSupportState.UNKNOWN.value
        elif row["group"] == "reference":
            packet["support_state"] = ModelSupportState.NOT_APPLICABLE.value
            packet["evidence_state"] = EvidenceState.NOT_APPLICABLE.value
        elif row["lineage"] not in eligible_set:
            packet["support_state"] = ModelSupportState.UNKNOWN.value
        elif row["raw_value"] <= -0.5:
            packet["support_state"] = ModelSupportState.SUPPORT.value
        else:
            packet["support_state"] = ModelSupportState.CONFLICT.value

    denominator = {"snapshot_models": len(models), "mutation_assayed_models": sum(row["group"] is not None for row in preliminary),
                   "dependency_measured_models": sum(row["raw_value"] is not None for row in preliminary), "evaluable_models": sum(row["evaluable"] for row in preliminary),
                   "context_models": len(context_rows), "reference_models": len(reference_rows), "total_models": len(analysis_rows), "eligible_lineages": len(eligible_lineages)}
    limitations = [
        "Observed processed mutation and Chronos evidence is internal processed-data evaluation, not independent biological validation.",
        _REFERENCE_LABEL,
        "Screen-specific dependency QC is unavailable in the frozen BR02 output; measured target Chronos is analysed without inventing a QC pass.",
        "The route does not change native D, rank, tier, veto, candidate partitions, filters, preprocessing, scores, or ranks.",
    ]
    if common_essential is True:
        limitations.append("The target is flagged common-essential; universal essentiality is not selective conditional evidence and this route cannot qualify.")
    elif common_essential is None:
        limitations.append("The frozen adapter did not provide a resolvable target common-essential flag, so this route cannot qualify.")
    result = {
        "schema_version": CONDITIONAL_DEPENDENCY_SCHEMA_VERSION, "research_query_id": query_id,
        "route_family": RouteFamily.CONDITIONAL_DEPENDENCY.value, "route_policy_version": policy_version,
        "relationship": relationship, "relationship_id": relationship_id, "route_status": status,
        "route_qualified": status == RouteStatus.SUPPORTED.value,
        "qualification_statistic": {"name": "lineage_adjusted_mutation_context_chronos_shift", "value": shift,
                                      "unit": "Chronos dependency score", "direction": "median(reference adjusted Chronos) - median(context adjusted Chronos); positive means context more target-dependent.",
                                      "raw_group_medians": {"context": float(np.median(raw_context)) if raw_context else None, "reference": float(np.median(raw_reference)) if raw_reference else None},
                                      "adjusted_group_medians": {"context": float(np.median(adjusted_context)) if adjusted_context else None, "reference": float(np.median(adjusted_reference)) if adjusted_reference else None},
                                      "gates": gates,
                                      "bootstrap": {"replicates": 500, "seed": 25, "stratification": "lineage x mutation context/reference group", "confidence_interval_95": ci},
                                      "leave_one_lineage_out": {"values": loo, "sign_consistency": consistency, "minimum": 0.80, "unresolved_when_no_estimable_holdout": not estimable_loo}},
        "qualification_direction": "Positive adjusted shift means context more target-dependent.", "denominator": denominator,
        "model_packets": [packets_by_id[model_id] for model_id in sorted(packets_by_id)],
        "evidence_ledger": sorted((dict(row) for row in ledger), key=lambda row: (str(row.get("ModelID", "")), str(row.get("evidence_id", "")))),
        "source": {"route_input_sources": list(adapter_result.get("route_input_sources", [])), "adapter_schema_version": adapter_result.get("schema_version"),
                   "data_id": "processed_mutations and processed_dependency retained source rows via BR02 evidence ledger", "units": {"mutation": "processed mutation event", "target_dependency": "Chronos dependency score"}},
        "qc": {"problematic_models_excluded_from_qualification_and_model_support": True, "screen_specific_dependency_qc": "screen-specific QC unavailable in frozen BR02 output; no QC pass is inferred."},
        "common_essential": {"target": dict(essential) if isinstance(essential, Mapping) else None, "is_common_essential": common_essential,
                             "warning": "Common-essential targets never qualify this route."},
        "missingness": {"mutation_not_assayed_models": sum(row["group"] is None for row in preliminary), "dependency_unresolved_models": sum(row["raw_value"] is None for row in preliminary),
                        "problematic_models": sum(packet["problematic_model_state"] == "problematic" for packet in packets_by_id.values()), "problematic_state_unknown_models": sum(packet["problematic_model_state"] == "unknown" for packet in packets_by_id.values()),
                        "lineage_ineligible_models": sum(row["evaluable"] and row["lineage"] not in eligible_set for row in preliminary)},
        "limitations": limitations,
    }
    safe = to_json_safe(result)
    json.dumps(safe, sort_keys=True, separators=(",", ":"), allow_nan=False)
    if canonical_bytes(adapter_result) != original_bytes:
        raise RuntimeError("BR04 conditional-dependency evaluation attempted to mutate immutable adapter evidence")
    return safe


def _rna_context_rows(
    models: Mapping[str, Mapping[str, Any]], ledger: list[Any], source_id: str, target_id: str,
) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]], dict[str, list[dict[str, Any]]]]:
    """Read only BR02 source-RNA and target-dependency rows, retaining their raw finite values."""
    rna_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    dependency_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in ledger:
        if not isinstance(row, Mapping) or row.get("ModelID") not in models:
            raise ValueError("ledger ModelID must belong to immutable model_spine")
        copied = dict(row)
        if copied.get("source_id") == "processed_expression_rna" and copied.get("ensembl_id") == source_id:
            rna_rows[copied["ModelID"]].append(copied)
        elif copied.get("source_id") == "processed_dependency" and copied.get("ensembl_id") == target_id:
            dependency_rows[copied["ModelID"]].append(copied)

    preliminary: list[dict[str, Any]] = []
    for model_id in sorted(models):
        rna = sorted(rna_rows[model_id], key=lambda row: (str(row.get("evidence_id", "")), canonical_bytes(row)))
        dependency = sorted(dependency_rows[model_id], key=lambda row: (str(row.get("evidence_id", "")), canonical_bytes(row)))
        observed_rna = [row for row in rna if row.get("evidence_state") == EvidenceState.OBSERVED.value and _finite_number(row.get("raw_value")) is not None]
        observed_dependency = [row for row in dependency if row.get("evidence_state") == EvidenceState.OBSERVED.value and _finite_number(row.get("raw_value")) is not None]
        rna_value = _finite_number(observed_rna[0]["raw_value"]) if len(observed_rna) == 1 else None
        dependency_value = _finite_number(observed_dependency[0]["raw_value"]) if len(observed_dependency) == 1 else None
        model = models[model_id]
        lineage = model.get("lineage") if isinstance(model.get("lineage"), str) and model.get("lineage") else None
        problematic = _known_problematic(model.get("is_problematic"))
        preliminary.append({
            "ModelID": model_id, "lineage": lineage, "rna_value": rna_value,
            "raw_value": dependency_value, "problematic": problematic,
            "evaluable": rna_value is not None and dependency_value is not None and problematic is False and lineage is not None,
            "rna_evidence": rna, "dependency_evidence": dependency,
        })
    return preliminary, rna_rows, dependency_rows


def _rna_lineage_groups(preliminary: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Assign Q80/Q20 groups independently in each lineage without pooling any threshold."""
    by_lineage: dict[str, list[dict[str, Any]]] = defaultdict(list)
    all_lineages = sorted({row["lineage"] for row in preliminary if row["lineage"] is not None})
    for row in preliminary:
        if row["evaluable"]:
            by_lineage[row["lineage"]].append(row)

    thresholds: list[dict[str, Any]] = []
    for lineage in all_lineages:
        rows = by_lineage[lineage]
        eligible_n = len(rows)
        threshold = {
            "lineage": lineage, "eligible_n": eligible_n, "q20_log2tpm1": None,
            "q80_log2tpm1": None, "quantile_method": "NumPy linear", "non_discriminating_reason": None,
        }
        if eligible_n < 15:
            threshold["non_discriminating_reason"] = "eligible_n_below_15"
            thresholds.append(threshold)
            continue
        values = np.asarray(sorted(row["rna_value"] for row in rows), dtype=float)
        q20 = float(np.quantile(values, 0.20, method="linear"))
        q80 = float(np.quantile(values, 0.80, method="linear"))
        threshold["q20_log2tpm1"], threshold["q80_log2tpm1"] = q20, q80
        lower = {row["ModelID"] for row in rows if row["rna_value"] <= q20}
        upper = {row["ModelID"] for row in rows if row["rna_value"] >= q80}
        if q80 <= q20:
            threshold["non_discriminating_reason"] = "q80_lte_q20"
        elif lower & upper:
            threshold["non_discriminating_reason"] = "q20_q80_membership_overlap"
        else:
            for row in rows:
                row["group"] = "context" if row["ModelID"] in upper else "reference" if row["ModelID"] in lower else "middle"
        thresholds.append(threshold)
    return thresholds, [row for row in preliminary if row.get("group") in {"context", "reference", "middle"}]


def _evaluate_rna_high(adapter_result: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate locked within-lineage RNA Q80/Q20 context using only BR02 RNA and Chronos rows."""
    original_bytes = canonical_bytes(adapter_result)
    query_id, policy_version = adapter_result.get("research_query_id"), adapter_result.get("route_policy_version")
    if not isinstance(query_id, str) or not query_id or not isinstance(policy_version, str) or not policy_version:
        raise ValueError("resolved BR02 result requires query and policy versions")
    relationship = _relationship(adapter_result, ContextType.RNA_HIGH)
    source_id, target_id = relationship["source_gene"]["ensembl_id"], relationship["target_gene"]["ensembl_id"]
    relationship_id = adapter_result.get("relationship_id")
    if not isinstance(relationship_id, str) or not relationship_id:
        relationship_id = Relationship(
            RouteFamily.CONDITIONAL_DEPENDENCY,
            ResolvedGene(relationship["source_gene"]["symbol"], source_id),
            ResolvedGene(relationship["target_gene"]["symbol"], target_id), ContextType.RNA_HIGH,
        ).relationship_id(policy_version)
    spine, ledger = adapter_result.get("model_spine"), adapter_result.get("ledger")
    if not isinstance(spine, list) or not isinstance(ledger, list):
        raise ValueError("resolved BR02 conditional evidence requires model_spine and ledger lists")
    models: dict[str, dict[str, Any]] = {}
    for row in spine:
        if not isinstance(row, Mapping) or not isinstance(row.get("ModelID"), str) or not row["ModelID"]:
            raise ValueError("model_spine requires non-empty unique ModelID values")
        models[row["ModelID"]] = dict(row)
    if len(models) != len(spine):
        raise ValueError("model_spine ModelID values must be unique")

    preliminary, _, _ = _rna_context_rows(models, ledger, source_id, target_id)
    thresholds, threshold_evaluable = _rna_lineage_groups(preliminary)
    threshold_by_lineage = {row["lineage"]: row for row in thresholds}
    eligible_lineages = sorted(row["lineage"] for row in thresholds if row["non_discriminating_reason"] is None)
    eligible_set = set(eligible_lineages)
    analysis_rows: list[dict[str, Any]] = []
    adjusted_by_model: dict[str, float] = {}
    for lineage in eligible_lineages:
        rows = [row for row in threshold_evaluable if row["lineage"] == lineage]
        median = float(np.median([row["raw_value"] for row in rows]))
        for row in rows:
            row["adjusted_value"] = float(row["raw_value"] - median)
            adjusted_by_model[row["ModelID"]] = row["adjusted_value"]
            if row["group"] in {"context", "reference"}:
                analysis_rows.append(row)

    context_rows = [row for row in analysis_rows if row["group"] == "context"]
    reference_rows = [row for row in analysis_rows if row["group"] == "reference"]
    raw_context = [row["raw_value"] for row in context_rows]
    raw_reference = [row["raw_value"] for row in reference_rows]
    adjusted_context = [row["adjusted_value"] for row in context_rows]
    adjusted_reference = [row["adjusted_value"] for row in reference_rows]
    shift = _statistic(analysis_rows)
    gates = {
        "context_minimum_n": len(context_rows) >= 10,
        "reference_minimum_n": len(reference_rows) >= 10,
        "total_minimum_n": len(analysis_rows) >= 30,
        "eligible_lineages_minimum": len(eligible_lineages) >= 2,
    }
    estimable = all(gates.values()) and shift is not None
    bootstrap = _stratified_bootstrap(analysis_rows, 500, 25) if estimable else []
    ci = [float(np.quantile(bootstrap, 0.025)), float(np.quantile(bootstrap, 0.975))] if bootstrap else [None, None]
    loo: list[dict[str, Any]] = []
    full_sign = 1 if shift is not None and shift > 0 else -1 if shift is not None and shift < 0 else 0
    for held_out in eligible_lineages:
        retained = [row for row in analysis_rows if row["lineage"] != held_out]
        retained_lineages = sorted({row["lineage"] for row in retained})
        retained_context = sum(row["group"] == "context" for row in retained)
        retained_reference = sum(row["group"] == "reference" for row in retained)
        holdout_value = _statistic(retained)
        holdout_estimable = (retained_context >= 10 and retained_reference >= 10 and len(retained) >= 30 and len(retained_lineages) >= 2 and holdout_value is not None)
        holdout_sign = 1 if holdout_value is not None and holdout_value > 0 else -1 if holdout_value is not None and holdout_value < 0 else 0
        loo.append({"held_out_lineage": held_out, "shift": holdout_value if holdout_estimable else None,
                    "estimable": holdout_estimable, "matches_full_sign": bool(holdout_estimable and full_sign != 0 and holdout_sign == full_sign),
                    "denominator": {"context_models": retained_context, "reference_models": retained_reference, "evaluable_models": len(retained), "eligible_lineages": len(retained_lineages)}})
    estimable_loo = [row for row in loo if row["estimable"]]
    consistency = (sum(row["matches_full_sign"] for row in estimable_loo) / len(estimable_loo)) if estimable_loo else None

    essential = (adapter_result.get("common_essential") or {}).get(target_id)
    common_essential = essential.get("is_common_essential") if isinstance(essential, Mapping) and isinstance(essential.get("is_common_essential"), bool) else None
    if not estimable:
        status = RouteStatus.UNAVAILABLE.value
    elif common_essential is True or common_essential is None:
        status = RouteStatus.EXPLORATORY.value
    elif shift > 0 and ci[0] is not None and ci[0] > 0 and consistency is not None and consistency >= 0.80:
        status = RouteStatus.SUPPORTED.value
    elif shift < 0 and ci[1] is not None and ci[1] < 0 and consistency is not None and consistency >= 0.80:
        status = RouteStatus.CONTRADICTED.value
    else:
        status = RouteStatus.EXPLORATORY.value

    packets: list[dict[str, Any]] = []
    for row in preliminary:
        threshold = threshold_by_lineage.get(row["lineage"]) if row["lineage"] is not None else None
        group = row.get("group")
        if row["problematic"] is not False or row["rna_value"] is None or row["raw_value"] is None or threshold is None or threshold["non_discriminating_reason"] is not None:
            support, evidence_state = ModelSupportState.UNKNOWN.value, EvidenceState.NOT_ASSAYED.value
        elif group in {"reference", "middle"}:
            support, evidence_state = ModelSupportState.NOT_APPLICABLE.value, EvidenceState.NOT_APPLICABLE.value
        elif row["raw_value"] <= -0.5:
            support, evidence_state = ModelSupportState.SUPPORT.value, EvidenceState.OBSERVED.value
        else:
            support, evidence_state = ModelSupportState.CONFLICT.value, EvidenceState.OBSERVED.value
        packets.append({
            "ModelID": row["ModelID"], "lineage": row["lineage"], "support_state": support, "evidence_state": evidence_state,
            "rna_context": "rna_high" if group == "context" else "rna_low" if group == "reference" else "middle" if group == "middle" else "unresolved",
            "rna_assay_state": "measured" if row["rna_value"] is not None else "unresolved_or_not_measured",
            "dependency_assay_state": "measured" if row["raw_value"] is not None else "unresolved_or_not_measured",
            "problematic_model_state": "problematic" if row["problematic"] is True else "not_problematic" if row["problematic"] is False else "unknown",
            "lineage_threshold": dict(threshold) if threshold is not None else None,
            "qc_state": "screen_specific_dependency_qc_unavailable_in_frozen_output",
            "route_native_value": {"raw_rna_log2tpm1": row["rna_value"], "rna_unit": "log2(TPM+1)", "raw_dependency_score": row["raw_value"], "raw_unit": "Chronos dependency score", "lineage_adjusted_dependency_score": adjusted_by_model.get(row["ModelID"]), "adjusted_unit": "Chronos dependency score", "is_probability": False},
            "evidence_ids": {"source_rna": [item.get("evidence_id") for item in row["rna_evidence"]], "dependency": [item.get("evidence_id") for item in row["dependency_evidence"]]},
            "limitation": "RNA context uses only frozen processed source log2(TPM+1) and target Chronos rows; screen-specific dependency QC is unavailable and no QC pass is inferred.",
        })
    relevant_ledger = [dict(row) for row in ledger if (row.get("source_id") == "processed_expression_rna" and row.get("ensembl_id") == source_id) or (row.get("source_id") == "processed_dependency" and row.get("ensembl_id") == target_id)]
    denominator = {"snapshot_models": len(models), "rna_measured_models": sum(row["rna_value"] is not None for row in preliminary), "dependency_measured_models": sum(row["raw_value"] is not None for row in preliminary), "evaluable_models": sum(row["evaluable"] for row in preliminary), "context_models": len(context_rows), "reference_models": len(reference_rows), "middle_models": sum(row.get("group") == "middle" for row in threshold_evaluable), "total_models": len(analysis_rows), "eligible_lineages": len(eligible_lineages)}
    limitations = [
        "Observed processed RNA and Chronos evidence is internal processed-data evaluation, not independent biological validation.",
        "RNA thresholds are standard linear Q20/Q80 calculated independently within each exact lineage; HPA, GEO, protein, graph and validation data are not route inputs.",
        "Screen-specific dependency QC is unavailable in the frozen BR02 output; measured target Chronos is analysed without inventing a QC pass.",
        "The route does not change native D, rank, tier, veto, candidate partitions, filters, preprocessing, scores, or ranks.",
    ]
    if common_essential is True:
        limitations.append("The target is flagged common-essential; universal essentiality is not selective conditional evidence and this route cannot qualify.")
    elif common_essential is None:
        limitations.append("The frozen adapter did not provide a resolvable target common-essential flag, so this route cannot qualify.")
    result = {
        "schema_version": CONDITIONAL_DEPENDENCY_SCHEMA_VERSION, "research_query_id": query_id, "route_family": RouteFamily.CONDITIONAL_DEPENDENCY.value,
        "route_policy_version": policy_version, "relationship": relationship, "relationship_id": relationship_id, "route_status": status, "route_qualified": status == RouteStatus.SUPPORTED.value,
        "qualification_statistic": {"name": "lineage_adjusted_rna_high_context_chronos_shift", "value": shift, "unit": "Chronos dependency score", "direction": "median(reference adjusted Chronos) - median(context adjusted Chronos); positive means RNA-high context more target-dependent.", "raw_group_medians": {"context": float(np.median(raw_context)) if raw_context else None, "reference": float(np.median(raw_reference)) if raw_reference else None}, "adjusted_group_medians": {"context": float(np.median(adjusted_context)) if adjusted_context else None, "reference": float(np.median(adjusted_reference)) if adjusted_reference else None}, "gates": gates, "bootstrap": {"replicates": 500, "seed": 25, "stratification": "lineage x RNA Q80/Q20 context/reference group", "confidence_interval_95": ci}, "leave_one_lineage_out": {"values": loo, "sign_consistency": consistency, "minimum": 0.80, "unresolved_when_no_estimable_holdout": not estimable_loo}},
        "qualification_direction": "Positive adjusted shift means RNA-high context more target-dependent.", "denominator": denominator, "lineage_thresholds": thresholds,
        "model_packets": sorted(packets, key=lambda packet: packet["ModelID"]), "evidence_ledger": sorted(relevant_ledger, key=lambda row: (str(row.get("ModelID", "")), str(row.get("evidence_id", "")))),
        "source": {"route_input_sources": ["processed_expression_rna", "processed_dependency"], "adapter_schema_version": adapter_result.get("schema_version"), "data_id": "processed_expression_rna and processed_dependency retained source rows via BR02 evidence ledger", "units": {"source_rna": "log2(TPM+1)", "target_dependency": "Chronos dependency score"}},
        "qc": {"problematic_models_excluded_from_qualification_and_model_support": True, "screen_specific_dependency_qc": "screen-specific QC unavailable in frozen BR02 output; no QC pass is inferred."},
        "common_essential": {"target": dict(essential) if isinstance(essential, Mapping) else None, "is_common_essential": common_essential, "warning": "Common-essential targets never qualify this route."},
        "missingness": {"rna_unresolved_models": sum(row["rna_value"] is None for row in preliminary), "dependency_unresolved_models": sum(row["raw_value"] is None for row in preliminary), "problematic_models": sum(row["problematic"] is True for row in preliminary), "problematic_state_unknown_models": sum(row["problematic"] is None for row in preliminary), "lineage_ineligible_models": sum(row["evaluable"] and (threshold_by_lineage.get(row["lineage"], {}).get("non_discriminating_reason") is not None) for row in preliminary)},
        "limitations": limitations,
    }
    safe = to_json_safe(result)
    json.dumps(safe, sort_keys=True, separators=(",", ":"), allow_nan=False)
    if canonical_bytes(adapter_result) != original_bytes:
        raise RuntimeError("BR05 conditional-dependency evaluation attempted to mutate immutable adapter evidence")
    return safe


def evaluate_conditional_dependency(adapter_result: Mapping[str, Any]) -> dict[str, Any]:
    """Dispatch the frozen directional mutation or RNA-high context without changing either policy."""
    if not isinstance(adapter_result, Mapping):
        raise ValueError("adapter_result must be a mapping")
    if adapter_result.get("route_family") != RouteFamily.CONDITIONAL_DEPENDENCY.value:
        raise ValueError("conditional evaluator requires a conditional_dependency adapter result")
    if adapter_result.get("evidence_state") == EvidenceState.UNRESOLVABLE.value:
        return _unresolvable_result(adapter_result)
    relationship = adapter_result.get("relationship")
    context = relationship.get("context_type") if isinstance(relationship, Mapping) else None
    if context == ContextType.MUTATION.value:
        return _evaluate_mutation(adapter_result)
    if context == ContextType.RNA_HIGH.value:
        return _evaluate_rna_high(adapter_result)
    raise ValueError("conditional evaluator requires mutation or rna_high context")


evaluate_conditional_dependency_route = evaluate_conditional_dependency
