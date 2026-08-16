"""Outcome-blind planning for AUTO01's fixed automatic biological-route battery.

The planner reads only the canonical query retained in a trusted saved snapshot.  It does not
read native result partitions, desirability values, candidate lists, route attachments, source
data, or route outcomes; AUTO03 will later persist execution against this deterministic state.
"""

from __future__ import annotations

import hashlib
import json
import re
from itertools import combinations
from typing import Any, Mapping

from .registry import automatic_policy_bytes, load_automatic_route_policy


PLAN_SCHEMA_VERSION = "automatic-biological-route-plan-v1"
PLAN_STATUS = "planned"
PENDING_STATUS = "pending"
NOT_APPLICABLE_STATUS = "not_applicable"
SINGLE_GENE_REASON = "not_applicable_requires_two_inclusion_genes"
PAIR_LIMIT_REASON = "automatic_pair_limit"
_QUERY_ID_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_TASK_STATUS = frozenset({"pending", "completed", "failed"})


def _canonical_bytes(value: Any) -> bytes:
    """Encode an ID payload in one strict JSON form so hashes are reproducible."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _stable_id(prefix: str, payload: Mapping[str, Any]) -> str:
    """Return a readable, deterministic identifier without introducing a clock or counter."""
    return f"{prefix}_{hashlib.sha256(_canonical_bytes(payload)).hexdigest()}"


def _canonical_gene_ids(value: Any, field_name: str) -> list[str]:
    """Accept only canonical Ensembl-shaped identifiers already stored in the snapshot.

    Symbol resolution belongs to native query capture, before the snapshot exists.  This narrow
    check rejects an unresolved symbol while retaining synthetic Ensembl-shaped IDs used in
    extension fixtures (for example ``ENSG-SYN``).
    """
    if not isinstance(value, list):
        raise ValueError(f"Snapshot query {field_name} must be a list of canonical Ensembl IDs")
    ids: set[str] = set()
    for gene_id in value:
        if not isinstance(gene_id, str) or gene_id != gene_id.strip() or not gene_id.startswith("ENSG"):
            raise ValueError(f"Snapshot query {field_name} contains a non-canonical Ensembl ID")
        ids.add(gene_id)
    return sorted(ids)


def _canonical_snapshot_query(snapshot: Mapping[str, Any]) -> tuple[str, list[str], list[str]]:
    """Validate the saved query boundary without falling back to native output fields.

    Requiring the content-addressed saved query makes a symbol-only browser request insufficient
    for planning and keeps a malformed snapshot from silently creating a different route plan.
    """
    if not isinstance(snapshot, Mapping):
        raise ValueError("Automatic planning requires a saved snapshot mapping")
    query_id = snapshot.get("query_id")
    if not isinstance(query_id, str) or not _QUERY_ID_PATTERN.fullmatch(query_id):
        raise ValueError("Automatic planning requires a saved snapshot query_id")
    query = snapshot.get("query")
    if not isinstance(query, Mapping):
        raise ValueError("Automatic planning requires the canonical saved snapshot query")
    inclusion_ids = _canonical_gene_ids(query.get("inclusion_genes"), "inclusion_genes")
    exclusion_ids = _canonical_gene_ids(query.get("exclusion_genes"), "exclusion_genes")
    if not inclusion_ids:
        raise ValueError("Automatic planning requires at least one canonical inclusion gene")
    overlap = sorted(set(inclusion_ids).intersection(exclusion_ids))
    if overlap:
        raise ValueError("Canonical inclusion and exclusion genes cannot overlap")
    return query_id, inclusion_ids, exclusion_ids


def _pair_record(plan_id: str, first_id: str, second_id: str) -> dict[str, str]:
    """Build the serializable identity shared by a retained pair and all six of its tasks."""
    pair = [first_id, second_id]
    return {
        "pair_id": _stable_id("automatic_pair", {"plan_id": plan_id, "canonical_gene_pair": pair}),
        "first_ensembl_id": first_id,
        "second_ensembl_id": second_id,
    }


def _task_record(plan_id: str, pair: Mapping[str, str], declaration: Mapping[str, Any]) -> dict[str, Any]:
    """Translate one locked policy declaration into one pending, direction-explicit task."""
    direction = declaration["pair_direction"]
    first_id, second_id = pair["first_ensembl_id"], pair["second_ensembl_id"]
    source_id, target_id = (second_id, first_id) if direction == "second_to_first" else (first_id, second_id)
    task_identity = {
        "plan_id": plan_id,
        "pair_id": pair["pair_id"],
        "declared_task_id": declaration["task_id"],
        "source_ensembl_id": source_id,
        "target_ensembl_id": target_id,
    }
    return {
        "task_id": _stable_id("automatic_task", task_identity),
        "pair_id": pair["pair_id"],
        "canonical_gene_pair": [first_id, second_id],
        "declared_task_id": declaration["task_id"],
        "route_family": declaration["route_family"],
        "context_type": declaration["context_type"],
        "pair_direction": direction,
        "source_ensembl_id": source_id,
        "target_ensembl_id": target_id,
        "status": PENDING_STATUS,
        "reason": "automatic_fixed_battery",
    }


def _plan_identity(plan: Mapping[str, Any]) -> dict[str, Any]:
    """Keep mutable execution state out of the deterministic plan identity."""
    return {
        "schema_version": plan.get("schema_version"),
        "plan_id": plan.get("plan_id"),
        "query_id": plan.get("query_id"),
        "policy_version": plan.get("policy_version"),
        "policy_identity": plan.get("policy_identity"),
        "invocation": plan.get("invocation"),
        "selection_timing": plan.get("selection_timing"),
        "gene_context": plan.get("gene_context"),
        "pair_routes": plan.get("pair_routes"),
        "pairs": plan.get("pairs"),
        "omitted_pairs": plan.get("omitted_pairs"),
        "native_ranking_boundary": plan.get("native_ranking_boundary"),
        "tasks": [{key: value for key, value in task.items() if key != "status"} for task in plan.get("tasks", [])],
    }


def plan_identity(plan: Mapping[str, Any]) -> dict[str, Any]:
    """Expose the immutable persisted-plan fields for AUTO03 fail-closed comparisons."""
    if not isinstance(plan, Mapping):
        raise ValueError("automatic execution requires an automatic plan object")
    return _plan_identity(plan)


def validate_persisted_plan(plan: Mapping[str, Any]) -> None:
    """Reject a malformed saved plan without consulting route outcomes or creating a new plan."""
    identity = plan_identity(plan)
    if identity["schema_version"] != PLAN_SCHEMA_VERSION:
        raise ValueError("automatic execution has an unsupported plan schema")
    if not isinstance(identity["plan_id"], str) or not identity["plan_id"].startswith("automatic_plan_"):
        raise ValueError("automatic execution plan has an invalid plan_id")
    if not isinstance(identity["policy_identity"], str) or not re.fullmatch(r"[0-9a-f]{64}", identity["policy_identity"]):
        raise ValueError("automatic execution plan has an invalid policy identity")
    if not isinstance(plan.get("tasks"), list) or not isinstance(plan.get("state"), Mapping):
        raise ValueError("automatic execution plan requires tasks and state")
    task_ids = []
    pair_ids = {pair.get("pair_id") for pair in plan.get("pairs", []) if isinstance(pair, Mapping)}
    for task in plan["tasks"]:
        if not isinstance(task, Mapping) or task.get("pair_id") not in pair_ids:
            raise ValueError("automatic execution plan task has an unknown pair")
        expected_id = _stable_id("automatic_task", {
            "plan_id": plan["plan_id"], "pair_id": task.get("pair_id"),
            "declared_task_id": task.get("declared_task_id"),
            "source_ensembl_id": task.get("source_ensembl_id"),
            "target_ensembl_id": task.get("target_ensembl_id"),
        })
        if task.get("task_id") != expected_id:
            raise ValueError("automatic execution plan task ID does not match its declaration")
        if task.get("status") not in _TASK_STATUS:
            raise ValueError("automatic execution plan task has an unsupported status")
        task_ids.append(task["task_id"])
    if len(task_ids) != len(set(task_ids)):
        raise ValueError("automatic execution plan has duplicate task IDs")
    counts = {status: sum(task["status"] == status for task in plan["tasks"]) for status in _TASK_STATUS}
    state = plan["state"]
    if (state.get("planned_task_count"), state.get("pending_task_count"),
            state.get("completed_task_count"), state.get("failed_task_count")) != (
        len(plan["tasks"]), counts["pending"], counts["completed"], counts["failed"]
    ):
        raise ValueError("automatic execution plan state counts do not match task state")
def create_automatic_plan(snapshot: Mapping[str, Any], policy_path=None) -> dict[str, Any]:
    """Create the complete fixed-battery plan from a canonical saved query and AUTO01 policy.

    The output is strict JSON data and includes mutable-by-AUTO03 task state, but this function
    itself creates no route evaluation, persistence operation, biological score, or native-result
    change.  Pair selection is solely sorted canonical inclusion IDs plus the declared cap.
    """
    policy = load_automatic_route_policy(policy_path)
    query_id, inclusion_ids, exclusion_ids = _canonical_snapshot_query(snapshot)
    plan_identity = {
        "schema_version": PLAN_SCHEMA_VERSION,
        "query_id": query_id,
        "policy_version": policy["policy_version"],
        "policy_identity": hashlib.sha256(automatic_policy_bytes(policy_path)).hexdigest(),
        "inclusion_genes": inclusion_ids,
        "exclusion_genes": exclusion_ids,
    }
    plan_id = _stable_id("automatic_plan", plan_identity)
    all_pairs = list(combinations(inclusion_ids, 2))
    retained_pairs = all_pairs[:policy["maximum_unordered_pairs"]]
    omitted_pairs = all_pairs[policy["maximum_unordered_pairs"]:]
    pairs = [_pair_record(plan_id, first_id, second_id) for first_id, second_id in retained_pairs]
    tasks = [
        _task_record(plan_id, pair, declaration)
        for pair in pairs
        for declaration in policy["fixed_task_battery"]
    ]

    if len(inclusion_ids) < 2:
        pair_routes = {
            "status": NOT_APPLICABLE_STATUS,
            "reason": SINGLE_GENE_REASON,
            "message": "Pair-specific automatic routes require two inclusion genes; gene context remains planned.",
        }
    else:
        pair_routes = {
            "status": PLAN_STATUS,
            "reason": "automatic_fixed_battery",
            "message": "Pair-specific routes were selected from canonical inclusion genes before route outcomes.",
        }

    return {
        "schema_version": PLAN_SCHEMA_VERSION,
        "plan_id": plan_id,
        "query_id": query_id,
        "policy_version": policy["policy_version"],
        "policy_identity": hashlib.sha256(automatic_policy_bytes(policy_path)).hexdigest(),
        "invocation": policy["invocation"],
        "selection_timing": policy["selection_timing"],
        "gene_context": {
            "status": PLAN_STATUS,
            "reason": "canonical_saved_query_genes",
            "inclusion_genes": inclusion_ids,
            "exclusion_genes": exclusion_ids,
        },
        "pair_routes": pair_routes,
        "pairs": pairs,
        "omitted_pairs": [
            {
                "canonical_gene_pair": [first_id, second_id],
                "status": "omitted",
                "reason": PAIR_LIMIT_REASON,
            }
            for first_id, second_id in omitted_pairs
        ],
        "tasks": tasks,
        "state": {
            "status": PLAN_STATUS,
            "planned_pair_count": len(pairs),
            "planned_task_count": len(tasks),
            "omitted_pair_count": len(omitted_pairs),
            "pending_task_count": len(tasks),
            "completed_task_count": 0,
            "failed_task_count": 0,
        },
        "native_ranking_boundary": policy["native_ranking_boundary"],
    }
