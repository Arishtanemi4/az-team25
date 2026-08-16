"""BR06 CRISPR codependency over one immutable BR02 evidence projection.

This route reports a descriptive, signed association of two frozen processed Chronos vectors.  It
does not reread processed data, construct a joint score, infer a mechanism, or alter native ranking.
"""

from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy.stats import spearmanr

from .contracts import EvidenceState, ModelSupportState, Relationship, ResolvedGene, RouteFamily, RouteStatus

_RESEARCH_SERVICE_DIR = str(Path(__file__).resolve().parents[1])
if _RESEARCH_SERVICE_DIR not in sys.path:
    sys.path.insert(0, _RESEARCH_SERVICE_DIR)

from serialization import canonical_bytes, to_json_safe  # noqa: E402 -- sys.path must be set up first


CODEPENDENCY_SCHEMA_VERSION = "br06-crispr-codependency-v1"
_DEPENDENCY_SOURCE = "processed_dependency"
_SEED = 25
_REPLICATES = 500
_LINEAGE_MINIMUM = 15
_SHARED_MINIMUM = 30
_LINEAGE_COUNT_MINIMUM = 2


def _finite_number(value: Any) -> float | None:
    """Return a finite non-boolean numeric measurement, never a coercion of missing evidence."""
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if math.isfinite(numeric) else None


def _known_problematic(value: Any) -> bool | None:
    """Require an explicit frozen model QC state before treating a model as non-problematic."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool) and int(value) in {0, 1}:
        return bool(int(value))
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "t", "yes", "y", "1"}:
            return True
        if normalized in {"false", "f", "no", "n", "0"}:
            return False
    return None


def _rho(first: list[float] | np.ndarray, second: list[float] | np.ndarray) -> float | None:
    """Calculate signed Spearman rho only for finite, non-constant vectors of equal size."""
    first_values = np.asarray(first, dtype=float)
    second_values = np.asarray(second, dtype=float)
    if len(first_values) < 2 or len(first_values) != len(second_values):
        return None
    if not np.isfinite(first_values).all() or not np.isfinite(second_values).all():
        return None
    if np.ptp(first_values) == 0 or np.ptp(second_values) == 0:
        return None
    value = float(spearmanr(first_values, second_values).statistic)
    return value if math.isfinite(value) and -1.0 <= value <= 1.0 else None


def _relationship(adapter_result: Mapping[str, Any]) -> tuple[dict[str, Any], tuple[str, str]]:
    """Validate the resolved symmetric pair while retaining its request orientation in output."""
    relationship = adapter_result.get("relationship")
    if not isinstance(relationship, Mapping):
        raise ValueError("resolved BR02 codependency evidence requires a relationship")
    if relationship.get("route_family") != RouteFamily.CRISPR_CODEPENDENCY.value:
        raise ValueError("codependency evaluator requires a crispr_codependency adapter result")
    if relationship.get("context_type") is not None:
        raise ValueError("crispr_codependency does not accept a context_type")
    source, target = relationship.get("source_gene"), relationship.get("target_gene")
    if not isinstance(source, Mapping) or not isinstance(target, Mapping):
        raise ValueError("relationship requires resolved source_gene and target_gene")
    source_id, target_id = source.get("ensembl_id"), target.get("ensembl_id")
    if not isinstance(source_id, str) or not source_id or not isinstance(target_id, str) or not target_id:
        raise ValueError("relationship genes require non-empty canonical Ensembl IDs")
    if source_id == target_id:
        raise ValueError("crispr codependency requires two distinct canonical Ensembl IDs")
    return dict(relationship), tuple(sorted((source_id, target_id)))


def _unresolvable_result(adapter_result: Mapping[str, Any]) -> dict[str, Any]:
    """Return unavailable identity output instead of guessing an unresolved gene's vector."""
    query_id = adapter_result.get("research_query_id")
    if not isinstance(query_id, str) or not query_id:
        raise ValueError("adapter result requires a non-empty research_query_id")
    return {
        "schema_version": CODEPENDENCY_SCHEMA_VERSION,
        "research_query_id": query_id,
        "route_family": RouteFamily.CRISPR_CODEPENDENCY.value,
        "route_policy_version": adapter_result.get("route_policy_version"),
        "relationship": None,
        "relationship_id": None,
        "route_status": RouteStatus.UNAVAILABLE.value,
        "route_qualified": False,
        "evidence_class": "exploratory functional association",
        "qualification_statistic": {
            "name": "lineage_median_adjusted_spearman_rho",
            "value": None,
            "raw_spearman_rho": None,
            "direction": "Signed association only; it is not a causal, mechanistic, physical-interaction, or synthetic-lethality claim.",
        },
        "denominator": {"snapshot_models": 0, "shared_evaluable_models": 0, "eligible_lineages": 0},
        "model_packets": [],
        "evidence_ledger": [],
        "source": {"route_input_sources": list(adapter_result.get("route_input_sources", []))},
        "qc": {"screen_specific_dependency_qc": "unavailable in frozen BR02 output; no QC pass is inferred."},
        "common_essential": {"source": None, "target": None, "source_is_common_essential": None, "target_is_common_essential": None,
                              "warning": "Either common-essential gene prevents route qualification."},
        "limitations": ["An unresolved or ambiguous frozen gene-reference entry is not guessed."],
    }


def _bootstrap_adjusted_rho(rows: list[dict[str, Any]]) -> tuple[list[float], int]:
    """Bootstrap adjusted rho 500 times by resampling rows independently within each lineage."""
    by_lineage: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_lineage[row["lineage"]].append(row)
    ordered = [(lineage, sorted(values, key=lambda item: item["ModelID"])) for lineage, values in sorted(by_lineage.items())]
    generator = np.random.default_rng(_SEED)
    values: list[float] = []
    for _ in range(_REPLICATES):
        first: list[float] = []
        second: list[float] = []
        for _, lineage_rows in ordered:
            indices = generator.integers(0, len(lineage_rows), size=len(lineage_rows))
            first.extend(lineage_rows[int(index)]["adjusted_effects"][0] for index in indices)
            second.extend(lineage_rows[int(index)]["adjusted_effects"][1] for index in indices)
        statistic = _rho(first, second)
        if statistic is not None:
            values.append(statistic)
    return values, _REPLICATES


def _holdouts(rows: list[dict[str, Any]], lineages: list[str], full_rho: float | None) -> tuple[list[dict[str, Any]], float | None]:
    """Expose every lineage holdout and apply the accepted conservative requalification gates."""
    full_sign = 1 if full_rho is not None and full_rho > 0 else -1 if full_rho is not None and full_rho < 0 else 0
    output: list[dict[str, Any]] = []
    for held_out in lineages:
        retained = [row for row in rows if row["lineage"] != held_out]
        retained_lineages = sorted({row["lineage"] for row in retained})
        statistic = _rho([row["adjusted_effects"][0] for row in retained], [row["adjusted_effects"][1] for row in retained])
        eligible = len(retained) >= _SHARED_MINIMUM and len(retained_lineages) >= _LINEAGE_COUNT_MINIMUM and statistic is not None
        sign = 1 if statistic is not None and statistic > 0 else -1 if statistic is not None and statistic < 0 else 0
        output.append({
            "held_out_lineage": held_out,
            "adjusted_spearman_rho": statistic if eligible else None,
            "estimable": eligible,
            "matches_full_sign": bool(eligible and full_sign != 0 and sign == full_sign),
            "denominator": {"shared_evaluable_models": len(retained), "eligible_lineages": len(retained_lineages)},
        })
    estimable = [row for row in output if row["estimable"]]
    consistency = sum(row["matches_full_sign"] for row in estimable) / len(estimable) if estimable else None
    return output, consistency


def evaluate_codependency(adapter_result: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate descriptive CRISPR codependency from a resolved immutable BR02 payload only."""
    if not isinstance(adapter_result, Mapping):
        raise ValueError("adapter_result must be a mapping")
    original_bytes = canonical_bytes(adapter_result)
    if adapter_result.get("route_family") != RouteFamily.CRISPR_CODEPENDENCY.value:
        raise ValueError("codependency evaluator requires a crispr_codependency adapter result")
    if adapter_result.get("evidence_state") == EvidenceState.UNRESOLVABLE.value:
        result = _unresolvable_result(adapter_result)
        json.dumps(to_json_safe(result), sort_keys=True, separators=(",", ":"), allow_nan=False)
        return result

    query_id, policy_version = adapter_result.get("research_query_id"), adapter_result.get("route_policy_version")
    if not isinstance(query_id, str) or not query_id or not isinstance(policy_version, str) or not policy_version:
        raise ValueError("resolved BR02 result requires query and policy versions")
    relationship, pair = _relationship(adapter_result)
    relationship_id = adapter_result.get("relationship_id")
    if not isinstance(relationship_id, str) or not relationship_id:
        relationship_id = Relationship(
            RouteFamily.CRISPR_CODEPENDENCY,
            ResolvedGene(str(relationship["source_gene"]["symbol"]), str(relationship["source_gene"]["ensembl_id"])),
            ResolvedGene(str(relationship["target_gene"]["symbol"]), str(relationship["target_gene"]["ensembl_id"])),
        ).relationship_id(policy_version)

    spine, ledger = adapter_result.get("model_spine"), adapter_result.get("ledger")
    if not isinstance(spine, list) or not isinstance(ledger, list):
        raise ValueError("resolved BR02 codependency evidence requires model_spine and ledger lists")
    models: dict[str, dict[str, Any]] = {}
    for row in spine:
        if not isinstance(row, Mapping) or not isinstance(row.get("ModelID"), str) or not row["ModelID"]:
            raise ValueError("model_spine requires non-empty unique ModelID values")
        models[row["ModelID"]] = dict(row)
    if len(models) != len(spine):
        raise ValueError("model_spine ModelID values must be unique")

    dependencies: dict[str, dict[str, list[dict[str, Any]]]] = defaultdict(lambda: defaultdict(list))
    for row in ledger:
        if not isinstance(row, Mapping) or row.get("ModelID") not in models:
            raise ValueError("ledger ModelID must belong to immutable model_spine")
        if row.get("source_id") == _DEPENDENCY_SOURCE and row.get("ensembl_id") in pair:
            dependencies[row["ModelID"]][row["ensembl_id"]].append(dict(row))

    preliminary: list[dict[str, Any]] = []
    packets_by_id: dict[str, dict[str, Any]] = {}
    for model_id in sorted(models):
        model = models[model_id]
        by_gene = {gene_id: sorted(dependencies[model_id][gene_id], key=lambda row: (str(row.get("evidence_id", "")), canonical_bytes(row))) for gene_id in pair}
        raw_effects: dict[str, float | None] = {}
        for gene_id, rows in by_gene.items():
            observed = [row for row in rows if row.get("evidence_state") == EvidenceState.OBSERVED.value and _finite_number(row.get("raw_value")) is not None]
            raw_effects[gene_id] = _finite_number(observed[0].get("raw_value")) if len(observed) == 1 else None
        problematic = _known_problematic(model.get("is_problematic"))
        lineage = model.get("lineage") if isinstance(model.get("lineage"), str) and model.get("lineage") else None
        evaluable = all(raw_effects[gene_id] is not None for gene_id in pair) and problematic is False and lineage is not None
        preliminary.append({"ModelID": model_id, "lineage": lineage, "problematic": problematic, "raw_effects": raw_effects, "evaluable": evaluable})
        packets_by_id[model_id] = {
            "ModelID": model_id, "lineage": lineage, "support_state": ModelSupportState.UNKNOWN.value,
            "evidence_state": EvidenceState.OBSERVED.value if all(raw_effects[gene_id] is not None for gene_id in pair) else EvidenceState.NOT_ASSAYED.value,
            "dependency_assay_state": {gene_id: "measured" if raw_effects[gene_id] is not None else "unresolved_or_not_measured" for gene_id in pair},
            "problematic_model_state": "problematic" if problematic is True else "not_problematic" if problematic is False else "unknown",
            "qc_state": "screen_specific_dependency_qc_unavailable_in_frozen_output",
            "route_native_value": {"raw_effects_by_gene": dict(raw_effects), "raw_unit": "Chronos dependency score",
                                   "lineage_adjusted_effects_by_gene": {gene_id: None for gene_id in pair}, "adjusted_unit": "Chronos dependency score", "is_probability": False},
            "evidence_ids": {gene_id: [row.get("evidence_id") for row in by_gene[gene_id]] for gene_id in pair},
            "limitation": "Screen-specific dependency QC is unavailable in the frozen adapter output; no QC pass is inferred.",
        }

    by_lineage: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in preliminary:
        if row["evaluable"]:
            by_lineage[row["lineage"]].append(row)
    eligible_lineages = sorted(lineage for lineage, rows in by_lineage.items() if len(rows) >= _LINEAGE_MINIMUM)
    eligible_set = set(eligible_lineages)
    analysis_rows: list[dict[str, Any]] = []
    for lineage in eligible_lineages:
        rows = by_lineage[lineage]
        medians = {gene_id: float(np.median([row["raw_effects"][gene_id] for row in rows])) for gene_id in pair}
        for row in rows:
            adjusted = tuple(float(row["raw_effects"][gene_id] - medians[gene_id]) for gene_id in pair)
            analysis = dict(row)
            analysis["adjusted_effects"] = adjusted
            analysis_rows.append(analysis)
            packets_by_id[row["ModelID"]]["route_native_value"]["lineage_adjusted_effects_by_gene"] = {pair[0]: adjusted[0], pair[1]: adjusted[1]}

    raw_rho = _rho([row["raw_effects"][pair[0]] for row in analysis_rows], [row["raw_effects"][pair[1]] for row in analysis_rows])
    adjusted_rho = _rho([row["adjusted_effects"][0] for row in analysis_rows], [row["adjusted_effects"][1] for row in analysis_rows])
    gates = {"shared_minimum_n": len(analysis_rows) >= _SHARED_MINIMUM, "eligible_lineages_minimum": len(eligible_lineages) >= _LINEAGE_COUNT_MINIMUM,
             "raw_vector_variance": raw_rho is not None, "adjusted_vector_variance": adjusted_rho is not None}
    estimable = all(gates.values())
    bootstrap, requested_replicates = _bootstrap_adjusted_rho(analysis_rows) if estimable else ([], _REPLICATES)
    ci = [float(np.quantile(bootstrap, 0.025)), float(np.quantile(bootstrap, 0.975))] if len(bootstrap) == requested_replicates else [None, None]
    loo, consistency = _holdouts(analysis_rows, eligible_lineages, adjusted_rho) if estimable else ([], None)

    essential_values = adapter_result.get("common_essential") if isinstance(adapter_result.get("common_essential"), Mapping) else {}
    source_essential, target_essential = essential_values.get(relationship["source_gene"]["ensembl_id"]), essential_values.get(relationship["target_gene"]["ensembl_id"])
    source_common = source_essential.get("is_common_essential") if isinstance(source_essential, Mapping) and isinstance(source_essential.get("is_common_essential"), bool) else None
    target_common = target_essential.get("is_common_essential") if isinstance(target_essential, Mapping) and isinstance(target_essential.get("is_common_essential"), bool) else None
    common_essential = source_common is True or target_common is True
    common_unknown = source_common is None or target_common is None
    interval_excludes_zero = ci[0] is not None and (ci[0] > 0 or ci[1] < 0)
    if not estimable:
        status = RouteStatus.UNAVAILABLE.value
    elif common_essential or common_unknown:
        status = RouteStatus.EXPLORATORY.value
    elif interval_excludes_zero and consistency is not None and consistency >= 0.80:
        status = RouteStatus.SUPPORTED.value
    else:
        status = RouteStatus.EXPLORATORY.value

    for row in preliminary:
        packet = packets_by_id[row["ModelID"]]
        first, second = row["raw_effects"][pair[0]], row["raw_effects"][pair[1]]
        if row["problematic"] is not False or first is None or second is None or row["lineage"] not in eligible_set:
            packet["support_state"] = ModelSupportState.UNKNOWN.value
        elif first <= -0.5 and second <= -0.5:
            packet["support_state"] = ModelSupportState.SUPPORT.value
        elif first <= -0.5 or second <= -0.5:
            packet["support_state"] = ModelSupportState.NEUTRAL.value
        else:
            packet["support_state"] = ModelSupportState.CONFLICT.value

    relevant_ledger = [dict(row) for row in ledger if row.get("source_id") == _DEPENDENCY_SOURCE and row.get("ensembl_id") in pair]
    limitations = [
        "CRISPR codependency is exploratory functional association evidence even when this internal route_status is supported; correlation is not causal, mechanistic, physical-interaction, or synthetic-lethality evidence.",
        "Observed processed Chronos evidence is internal processed-data evaluation, not independent biological validation.",
        "Screen-specific dependency QC is unavailable in the frozen BR02 output; measured values are analysed without claiming a QC pass.",
        "The route does not change native D, rank, tier, veto, candidate partitions, filters, preprocessing, scores, or ranks.",
    ]
    if common_essential:
        limitations.append("At least one gene is flagged common-essential; this route cannot qualify.")
    elif common_unknown:
        limitations.append("At least one frozen common-essential flag is unresolved; this route cannot qualify.")
    result = {
        "schema_version": CODEPENDENCY_SCHEMA_VERSION, "research_query_id": query_id, "route_family": RouteFamily.CRISPR_CODEPENDENCY.value,
        "route_policy_version": policy_version, "relationship": relationship, "relationship_id": relationship_id, "route_status": status,
        "route_qualified": status == RouteStatus.SUPPORTED.value, "evidence_class": "exploratory functional association",
        "qualification_statistic": {
            "name": "lineage_median_adjusted_spearman_rho", "value": adjusted_rho, "raw_spearman_rho": raw_rho,
            "direction": "Signed association only; positive and negative rho are both descriptive and neither is a causal, mechanistic, physical-interaction, or synthetic-lethality claim.",
            "gates": gates, "bootstrap": {"replicates": _REPLICATES, "estimable_replicates": len(bootstrap), "seed": _SEED, "stratification": "lineage", "confidence_interval_95": ci},
            "leave_one_lineage_out": {"values": loo, "sign_consistency": consistency, "minimum": 0.80, "requalification_rule": "retained shared N >= 30 and >= 2 eligible lineages with non-constant adjusted vectors", "unresolved_when_no_estimable_holdout": not any(row["estimable"] for row in loo)},
        },
        "qualification_direction": "The adjusted Spearman sign is retained; no directional biological hypothesis is imposed.",
        "denominator": {"snapshot_models": len(models), "both_dependency_measured_models": sum(all(row["raw_effects"][gene_id] is not None for gene_id in pair) for row in preliminary),
                        "shared_evaluable_models": len(analysis_rows), "eligible_lineages": len(eligible_lineages), "lineage_minimum_n": _LINEAGE_MINIMUM},
        "model_packets": [packets_by_id[model_id] for model_id in sorted(packets_by_id)],
        "evidence_ledger": sorted(relevant_ledger, key=lambda row: (str(row.get("ModelID", "")), str(row.get("ensembl_id", "")), str(row.get("evidence_id", "")))),
        "source": {"route_input_sources": ["processed_dependency"], "adapter_schema_version": adapter_result.get("schema_version"), "data_id": "processed_dependency retained source rows via BR02 evidence ledger", "units": {"dependency": "Chronos dependency score"}},
        "qc": {"problematic_models_excluded_from_qualification_and_model_support": True, "screen_specific_dependency_qc": "screen-specific QC unavailable in frozen BR02 output; no QC pass is inferred."},
        "common_essential": {"source": dict(source_essential) if isinstance(source_essential, Mapping) else None, "target": dict(target_essential) if isinstance(target_essential, Mapping) else None,
                             "source_is_common_essential": source_common, "target_is_common_essential": target_common, "warning": "Either common-essential gene prevents route qualification."},
        "missingness": {"source_dependency_unresolved_models": sum(row["raw_effects"][pair[0]] is None for row in preliminary), "target_dependency_unresolved_models": sum(row["raw_effects"][pair[1]] is None for row in preliminary),
                        "problematic_models": sum(row["problematic"] is True for row in preliminary), "problematic_state_unknown_models": sum(row["problematic"] is None for row in preliminary),
                        "lineage_ineligible_models": sum(row["evaluable"] and row["lineage"] not in eligible_set for row in preliminary)},
        "limitations": limitations,
    }
    safe = to_json_safe(result)
    json.dumps(safe, sort_keys=True, separators=(",", ":"), allow_nan=False)
    if canonical_bytes(adapter_result) != original_bytes:
        raise RuntimeError("BR06 codependency evaluation attempted to mutate immutable adapter evidence")
    return safe


evaluate_codependency_route = evaluate_codependency
