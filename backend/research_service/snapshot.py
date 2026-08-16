"""Immutable query snapshots (extensions/plans/CONTRACTS.md C2). A snapshot is the one evidence
object every later extension stage (S04-S09) reads -- never a fresh scoring call, never a
client-supplied score JSON (C2: "Backend creates snapshots from its own results"). The native
result inside it is deep-copied once at creation and never mutated again; derived fields
(context, alternatives, comparison -- later stages) sit beside it, never folded into it.
"""

import copy
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

# Explicit sys.path insertion + bare-name import, matching this codebase's convention
# throughout scoring/, backend/, and backend/research_service/tests/ (none of it relies on relative package
# imports) -- rather than `from . import serialization`, which would require this module to
# always be loaded as part of the `backend.research_service` package.
_RESEARCH_DIR = str(Path(__file__).resolve().parent)
if _RESEARCH_DIR not in sys.path:
    sys.path.insert(0, _RESEARCH_DIR)

import serialization  # noqa: E402 -- sys.path must be set up first

SCHEMA_VERSION = "research-query-v1"

_VALID_QUERY_ID = re.compile(r"^[0-9a-f]{64}$")

# Keys stripped from `manifest` before hashing, so a caller's own wall-clock stamp on the
# manifest can never make an otherwise-identical query/data/method produce a different ID.
_TIMESTAMP_KEYS = frozenset({"created_at", "generated_at", "timestamp", "built_at"})

# backend/research_service/.gitignore excludes this whole tree -- generated, per-machine, rebuildable from the
# same query against the same baseline.
_DEFAULT_SNAPSHOT_DIR = Path(__file__).resolve().parent / "runtime" / "snapshots"


def _snapshot_dir(snapshot_dir=None):
    return Path(snapshot_dir) if snapshot_dir is not None else _DEFAULT_SNAPSHOT_DIR


def compute_query_id(native_result, manifest):
    """SHA-256 over the canonical (query+filters, input/method identity, native output) --
    excluding any timestamp key `manifest` may carry, and excluding `created_at` implicitly
    since `native_result` (scoring/export.py's own return value) never carries one. The query
    itself is read from `native_result["query"]` rather than passed separately -- S02's export
    contract already carries `inclusion_genes`/`exclusion_genes`/`filters` there (C2's
    ten-item output contract item 1), so there is one place a query's shape is defined."""
    manifest_without_timestamps = {
        key: value for key, value in manifest.items() if key not in _TIMESTAMP_KEYS
    }
    identity_payload = {"native_result": native_result, "manifest": manifest_without_timestamps}
    return serialization.content_hash(identity_payload)


def create_snapshot(native_result, manifest, model_metadata):
    """Assembles one immutable snapshot.

    `native_result` is deep-copied here, once, so no caller-held reference to the original dict
    can mutate a snapshot after creation (C2: "preserve native_result; derived payloads are
    separate and cannot mutate it"). `context`/`alternatives` are left `None` -- S03's own scope
    is the snapshot object and the read-only adapter, not populating them; S04-S08 fill them in
    beside this immutable core, never inside `native_result` itself, so a graph/context failure
    can never erase or alter `D` (C2: "graph failure must not erase D")."""
    native_copy = copy.deepcopy(native_result)
    query_id = compute_query_id(native_copy, manifest)
    return {
        "schema_version": SCHEMA_VERSION,
        "query_id": query_id,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "manifest": manifest,
        "query": native_copy.get("query"),
        "native_result": native_copy,
        "models": model_metadata,
        "context": None,
        "alternatives": None,
        "limitations": [],
    }


def _validate_query_id(query_id):
    """Rejects anything that is not exactly a 64-character lowercase hex string before it is
    ever used to build a filesystem path -- a `../`, an absolute path, or a null byte cannot
    survive this regex, closing the traversal case C2/S03 both name explicitly."""
    if not isinstance(query_id, str) or not _VALID_QUERY_ID.match(query_id):
        raise ValueError(
            f"Invalid snapshot query_id (must be exactly 64 lowercase hex characters): "
            f"{query_id!r}"
        )


def save_snapshot(snapshot, snapshot_dir=None):
    """Writes one snapshot to <snapshot_dir>/<query_id>.json as strict JSON (`allow_nan=False`).
    The path is content-addressed, so saving the same snapshot twice writes the same bytes to
    the same path -- an overwrite is a no-op, never a race between two versions of one query."""
    query_id = snapshot["query_id"]
    _validate_query_id(query_id)
    directory = _snapshot_dir(snapshot_dir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{query_id}.json"
    safe = serialization.to_json_safe(snapshot)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(safe, handle, indent=2, allow_nan=False)
    return path


class SnapshotNotFoundError(LookupError):
    """A well-formed but unknown query_id. C2: "Unknown ID -> 404 with rerun instruction" -- the
    API layer (a later stage, S09) is expected to catch this and translate it to an HTTP 404;
    this package only distinguishes "malformed ID" (ValueError, checked first) from "valid ID,
    nothing saved under it" (this), since it has no HTTP layer of its own to answer with."""


def load_snapshot(query_id, snapshot_dir=None):
    """Reads back a previously saved snapshot by its content-addressed ID."""
    _validate_query_id(query_id)
    directory = _snapshot_dir(snapshot_dir)
    path = directory / f"{query_id}.json"
    if not path.exists():
        raise SnapshotNotFoundError(f"No snapshot saved for query_id {query_id}")
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)
