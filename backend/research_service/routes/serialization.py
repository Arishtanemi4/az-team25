"""Deterministic storage helpers for BR07 route attachments.

Route evaluators own scientific results.  This module only gives those results a stable,
extension-owned snapshot envelope; it deliberately adds no clock value or derived statistic.
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path
from typing import Any, Mapping

_RESEARCH_SERVICE_DIR = str(Path(__file__).resolve().parents[1])
if _RESEARCH_SERVICE_DIR not in sys.path:
    sys.path.insert(0, _RESEARCH_SERVICE_DIR)

from serialization import to_json_safe  # noqa: E402 -- sys.path must be set up first


ATTACHMENT_SCHEMA_VERSION = "br07-route-attachments-v1"


def empty_attachments() -> dict[str, Any]:
    """Return the one deterministic representation for a snapshot with no evaluated routes."""
    return {
        "schema_version": ATTACHMENT_SCHEMA_VERSION,
        "status": "not_computed",
        "by_relationship_id": {},
    }


def attachments(value: Any) -> dict[str, Any]:
    """Copy and validate the extension-owned attachment envelope without repairing data."""
    if value is None:
        return empty_attachments()
    if not isinstance(value, Mapping) or value.get("schema_version") != ATTACHMENT_SCHEMA_VERSION:
        raise ValueError("snapshot routes attachment has an unsupported schema")
    entries = value.get("by_relationship_id")
    if not isinstance(entries, Mapping):
        raise ValueError("snapshot routes attachment requires by_relationship_id")
    copied = copy.deepcopy(dict(value))
    copied["by_relationship_id"] = {
        key: copied["by_relationship_id"][key] for key in sorted(copied["by_relationship_id"])
    }
    return copied


def make_attachment(result: Mapping[str, Any]) -> dict[str, Any]:
    """Store a route result exactly once, retaining the evaluator's first request orientation."""
    relationship_id = result.get("relationship_id")
    policy_version = result.get("route_policy_version")
    relationship = result.get("relationship")
    if not isinstance(relationship_id, str) or not relationship_id:
        raise ValueError("route result requires a deterministic relationship_id")
    if not isinstance(policy_version, str) or not policy_version:
        raise ValueError("route result requires a route_policy_version")
    if not isinstance(relationship, Mapping):
        raise ValueError("route result requires its resolved relationship")
    return {
        "relationship_id": relationship_id,
        "route_policy_version": policy_version,
        "request_orientation": copy.deepcopy(dict(relationship)),
        "result": copy.deepcopy(to_json_safe(dict(result))),
    }


def put_attachment(value: Any, result: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any], bool]:
    """Persist by canonical ID or return the first immutable result for an idempotent POST."""
    envelope = attachments(value)
    entry = make_attachment(result)
    key = entry["relationship_id"]
    existing = envelope["by_relationship_id"].get(key)
    if existing is not None:
        if existing.get("route_policy_version") != entry["route_policy_version"]:
            raise ValueError("relationship ID is already attached under a different policy version")
        return envelope, copy.deepcopy(existing), True
    envelope["by_relationship_id"][key] = entry
    envelope["by_relationship_id"] = {
        item_key: envelope["by_relationship_id"][item_key]
        for item_key in sorted(envelope["by_relationship_id"])
    }
    envelope["status"] = "computed"
    return envelope, copy.deepcopy(entry), False


def summaries(value: Any) -> list[dict[str, Any]]:
    """Return bounded collection cards; complete ledgers remain in the explicit export."""
    result = []
    for relationship_id, entry in attachments(value)["by_relationship_id"].items():
        route = entry.get("result", {})
        result.append({
            "relationship_id": relationship_id,
            "route_policy_version": entry.get("route_policy_version"),
            "relationship": copy.deepcopy(entry.get("request_orientation")),
            "route_status": route.get("route_status"),
            "route_qualified": route.get("route_qualified"),
            "qualification_statistic": copy.deepcopy(route.get("qualification_statistic")),
            "denominator": copy.deepcopy(route.get("denominator")),
            "limitations": copy.deepcopy(route.get("limitations", [])),
        })
    return result
