"""Strict, deterministic JSON for the extension layer (extensions/plans/CONTRACTS.md C1/C2).

Two distinct needs, kept separate: `to_json_safe` makes an object writable as JSON at all
(numpy scalars, +-Infinity/NaN -> null) by reusing scoring/export.py's own **public**
`sanitize_for_json` -- the public export boundary an extension must use, per
extensions/README.md's "expose/use existing public export boundary... do not call private
sanitizers from new extension code" (scoring/export.py:S02 renamed `_sanitize` to this public
name for exactly this reuse). `canonical_bytes` additionally makes that JSON *reproducible*:
sorted keys and fixed separators, so the same logical object always hashes to the same bytes
regardless of dict insertion order -- what C2's content-addressed snapshot ID needs.
"""

import hashlib
import json
import sys
from pathlib import Path

_SCORING_DIR = str(Path(__file__).resolve().parents[2] / "scoring")
if _SCORING_DIR not in sys.path:
    sys.path.insert(0, _SCORING_DIR)

import export as scoring_export  # noqa: E402 -- sys.path must be set up first


def to_json_safe(obj):
    """Re-exports scoring/export.py's public sanitize_for_json under this module's own name, so
    extension code has one import path (`from backend.research_service import serialization`) rather
    than a reach into scoring/ scattered across every caller."""
    return scoring_export.sanitize_for_json(obj)


def canonical_bytes(obj):
    """UTF-8 bytes of `obj`'s JSON-safe form, sorted keys, fixed separators, `allow_nan=False`
    -- two calls with the same logical content always produce byte-identical output, which is
    what content_hash below needs to hash."""
    safe = to_json_safe(obj)
    return json.dumps(safe, sort_keys=True, separators=(",", ":"), allow_nan=False).encode(
        "utf-8"
    )


def content_hash(obj):
    """SHA-256 hex digest of `obj`'s canonical bytes."""
    return hashlib.sha256(canonical_bytes(obj)).hexdigest()
