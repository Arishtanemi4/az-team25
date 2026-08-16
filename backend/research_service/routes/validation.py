"""BR09 deterministic, validation-only route evaluation.

This module is deliberately downstream of the BR02--BR06 route engines.  It freezes the
registered controls and analytic identity before any result is calculated, writes only beneath
``backend/research_service/runtime/routes/evaluation``, and never supplies validation-only data to an engine.
It is an internal processed-data audit, not a calibration or a biological recommendation.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
import platform
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
from scipy.stats import spearmanr

from .codependency import evaluate_codependency
from .conditional_dependency import evaluate_conditional_dependency
from .direct_fusion import evaluate_direct_fusion
from .evidence_adapter import adapt_route_evidence
from .registry import (
    CONTROL_QUERIES_PATH,
    ROUTE_REGISTRY_PATH,
    VALIDATION_SOURCES_PATH,
    load_controls,
    load_route_registry,
)

_RESEARCH_SERVICE_DIR = str(Path(__file__).resolve().parents[1])
if _RESEARCH_SERVICE_DIR not in sys.path:
    sys.path.insert(0, _RESEARCH_SERVICE_DIR)

from serialization import canonical_bytes, to_json_safe  # noqa: E402 -- sys.path must be set up first


VALIDATION_SCHEMA_VERSION = "br09-validation-v1"
EVALUATION_DIRECTORY = Path(__file__).resolve().parents[1] / "runtime" / "routes" / "evaluation"
DEFAULT_DATA_DIRECTORY = Path(__file__).resolve().parents[3] / "data" / "processed"
DEFAULT_SNAPSHOT_DIRECTORY = Path(__file__).resolve().parents[1] / "runtime" / "snapshots"
DEFAULT_SNAPSHOT_ID = "30904f54abe2208d012df985723a4b465d6cab0cdf78fb3ce681fb16d2644f48"
_FIXTURE_KINDS = {"unrelated_fixture", "small_cohort_fixture", "common_essential_fixture", "missing_assay_fixture"}
_PREREGISTERED_MANIFEST = "preregistered_manifest.json"
_FINAL_REPAIR_CATEGORIES = [
    "post_freeze_output_serialization_code_hash_repair",
    "strict_final_manifest_drift_verification",
    "sensitivity_reconstruction_alignment_with_locked_route_statistics",
    "controlled_verdict_enum_normalization",
]


def _sha256_file(path: Path) -> str:
    """Hash a governed input byte-for-byte without changing it."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _lfs_pointer(path: Path) -> bool:
    """Recognise an unmaterialised Git LFS source before treating a stub as data."""
    if not path.is_file() or path.stat().st_size > 1024:
        return False
    return path.read_text(encoding="utf-8", errors="replace").startswith("version https://git-lfs.github.com/spec/v1")


def _file_identity(path: Path) -> dict[str, Any]:
    """Return deterministic input identity, including an explicit missing/LFS state."""
    if not path.exists():
        return {"path": path.as_posix(), "state": "absent", "sha256": None, "size_bytes": None}
    if _lfs_pointer(path):
        return {"path": path.as_posix(), "state": "lfs_pointer_unusable", "sha256": _sha256_file(path), "size_bytes": path.stat().st_size}
    return {"path": path.as_posix(), "state": "available", "sha256": _sha256_file(path), "size_bytes": path.stat().st_size}


def fixed_sensitivity_grid() -> dict[str, Any]:
    """Return the predeclared bounded grid; this is never selected after outcomes are seen."""
    return {
        "bootstrap_replicates": [250, 500, 1000],
        "dependency_support_boundary_model_only": [-0.6, -0.5, -0.4],
        "conditional_context_reference_minimum_n": [8, 10, 12],
        "conditional_total_minimum_n": [24, 30, 36],
        "codependency_shared_minimum_n": [24, 30, 36],
        "lineage_minimum_n": [12, 15, 18],
        "eligible_lineages_minimum": [2],
        "rna_quantile_pairs": [[0.15, 0.85], [0.20, 0.80], [0.25, 0.75]],
        "leave_one_lineage_out_stability_minimum": [0.70, 0.80, 0.90],
        "policy_status_rule": "one-factor-at-a-time around the locked base policy; model support boundary is model-only and never changes route status",
    }


def _route_code_paths() -> list[Path]:
    directory = Path(__file__).resolve().parent
    return sorted(path for path in directory.glob("*.py") if path.name != "validation.py")


def _governed_input_paths(root: Path, data_dir: Path, snapshot_id: str) -> list[Path]:
    """Return the complete fixed identity surface consumed by BR09 or its route engines."""
    return [
        data_dir / "cell_lines.parquet", data_dir / "gene_reference.parquet", data_dir / "coverage.parquet",
        data_dir / "fusions.parquet", data_dir / "mutations.parquet", data_dir / "expression_rna.parquet",
        data_dir / "dependency.parquet", root / "scoring" / "resources" / "common_essential_genes.json",
        root / "backend" / "research_service" / "runtime" / "snapshots" / f"{snapshot_id}.json",
    ]


def _manifest_surface(*, root: Path, snapshot_id: str, data_dir: Path) -> dict[str, Any]:
    """Calculate the current governed identity without writing or evaluating an outcome."""
    controls = load_controls()
    registry = load_route_registry()
    snapshot = root / "backend" / "research_service" / "runtime" / "snapshots" / f"{snapshot_id}.json"
    code_paths = [Path(__file__), *_route_code_paths()]
    return {
        "analytic_identity": {
            "controls_sha256": _sha256_file(CONTROL_QUERIES_PATH),
            "registry_sha256": _sha256_file(ROUTE_REGISTRY_PATH),
            "validation_sources_sha256": _sha256_file(VALIDATION_SOURCES_PATH),
            "policy_version": registry["policy_version"],
            "random_seed": 25,
            "snapshot_id": snapshot_id,
            "snapshot_sha256": _sha256_file(snapshot),
            "control_ids": [row["control_id"] for row in controls["controls"]],
            "sensitivity_grid": fixed_sensitivity_grid(),
        },
        "code_hashes": {path.relative_to(root).as_posix(): _sha256_file(path) for path in code_paths},
        "input_hashes": {item["path"]: item for item in (_file_identity(path) for path in _governed_input_paths(root, data_dir, snapshot_id))},
        "source_availability": discover_validation_sources(root),
    }


def assert_validation_import_boundary() -> None:
    """Prove engines do not depend on this validation-only module."""
    for path in _route_code_paths():
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imports.append(node.module)
        if any(name == "validation" or name.endswith(".validation") for name in imports):
            raise RuntimeError(f"route engine imports validation-only module: {path.name}")


def discover_validation_sources(root: Path) -> dict[str, Any]:
    """Inventory only governed Sanger/PRISM paths; no guessed treatment mapping is permitted."""
    prism_raw = root / "data" / "external" / "prism_repurposing" / "raw"
    prism_derived = root / "data" / "external" / "prism_repurposing" / "derived"
    prism_files = [_file_identity(path) for path in sorted([*prism_raw.glob("*"), *prism_derived.glob("*")]) if path.is_file()]
    usable_prism = [item for item in prism_files if item["state"] == "available"]
    if not prism_files:
        prism_state, prism_reason = "absent", "No governed PRISM source files were found."
    elif not usable_prism:
        prism_state, prism_reason = "unavailable", "All discovered PRISM files are Git LFS pointers and therefore have no usable metadata, response matrix, or gene-to-compound mapping."
    else:
        prism_state, prism_reason = "not_interpretable", "PRISM files exist, but this audit does not infer a gene-to-compound mapping; a governed explicit mapping is required before response can be interpreted."
    names = [path.as_posix() for path in root.joinpath("data").rglob("*") if path.is_file() and ("sanger" in path.name.lower() or "project_score" in path.name.lower())]
    sanger_state = "available" if names else "absent"
    return {
        "schema_version": VALIDATION_SCHEMA_VERSION,
        "source_role": "validation_only",
        "sanger_project_score": {"state": sanger_state, "paths": sorted(names), "reason": "No governed local Sanger Project Score file was found." if not names else "Present but not used as a route input."},
        "prism_repurposing": {"state": prism_state, "files": prism_files, "reason": prism_reason, "gene_to_compound_mapping": "unavailable"},
        "direct_fusion_independent_corroboration": {"state": "unavailable", "reason": "No governed external fusion source independent of the processed DepMap-derived fusion table was found locally. The processed source is not an independent corroboration."},
    }


def freeze_execution_manifest(
    output_dir: Path = EVALUATION_DIRECTORY,
    *,
    root: Path | None = None,
    snapshot_id: str = DEFAULT_SNAPSHOT_ID,
    data_dir: Path = DEFAULT_DATA_DIRECTORY,
) -> dict[str, Any]:
    """Freeze a chained final identity before a clean BR09 rerun, then write canonical JSON."""
    root = root or Path(__file__).resolve().parents[3]
    output_dir.mkdir(parents=True, exist_ok=True)
    preregistered = output_dir / _PREREGISTERED_MANIFEST
    if not preregistered.exists():
        raise RuntimeError("final BR09 manifest requires the preserved preregistered manifest")
    surface = _manifest_surface(root=root, snapshot_id=snapshot_id, data_dir=data_dir)
    manifest = {
        "schema_version": VALIDATION_SCHEMA_VERSION,
        **surface,
        "manifest_chain": {
            "preregistered_manifest_path": _PREREGISTERED_MANIFEST,
            "preregistered_manifest_sha256": _sha256_file(preregistered),
            "preregistration_scope": "controls, registry, validation sources, snapshot, seed and sensitivity grid were frozen before first outcome inspection",
            "final_manifest_scope": "frozen before the final clean control rerun, not before the first-ever outcome",
            "non_policy_repair_categories": _FINAL_REPAIR_CATEGORIES,
        },
        "runtime": {"python": sys.version.split()[0], "implementation": platform.python_implementation(), "platform": platform.platform(), "numpy": np.__version__},
    }
    payload = to_json_safe(manifest)
    (output_dir / "execution_manifest.json").write_bytes(canonical_bytes(payload) + b"\n")
    return payload


def _load_frozen_manifest(output_dir: Path) -> dict[str, Any]:
    path = output_dir / "execution_manifest.json"
    if not path.exists():
        raise RuntimeError("BR09 execution manifest must be frozen before evaluating controls")
    return json.loads(path.read_text(encoding="utf-8"))


def verify_execution_manifest(
    output_dir: Path = EVALUATION_DIRECTORY,
    *,
    snapshot_id: str = DEFAULT_SNAPSHOT_ID,
    data_dir: Path = DEFAULT_DATA_DIRECTORY,
) -> dict[str, Any]:
    """Fail closed if the final manifest or any governed identity drifted after freezing."""
    manifest = _load_frozen_manifest(output_dir)
    root = Path(__file__).resolve().parents[3]
    chain = manifest.get("manifest_chain")
    if not isinstance(chain, Mapping) or chain.get("non_policy_repair_categories") != _FINAL_REPAIR_CATEGORIES:
        raise RuntimeError("BR09 final manifest chain is missing or has been altered")
    preregistered = output_dir / str(chain.get("preregistered_manifest_path", ""))
    if not preregistered.is_file() or _sha256_file(preregistered) != chain.get("preregistered_manifest_sha256"):
        raise RuntimeError("BR09 preregistered manifest is absent or has changed")
    expected = _manifest_surface(root=root, snapshot_id=snapshot_id, data_dir=data_dir)
    for section in ("analytic_identity", "code_hashes", "input_hashes", "source_availability"):
        if manifest.get(section) != expected[section]:
            raise RuntimeError(f"BR09 final manifest {section} differs from current governed identity")
    return manifest


def _packet_state_counts(result: Mapping[str, Any], boundary: float = -0.5) -> dict[str, int]:
    """Count retained packet states; optional boundary is explicitly model-only."""
    counts: Counter[str] = Counter()
    family = result.get("route_family")
    for packet in result.get("model_packets", []):
        state = packet.get("support_state", "unknown")
        if boundary != -0.5 and family in {"conditional_dependency", "crispr_codependency"}:
            values = packet.get("route_native_value", {})
            if family == "conditional_dependency" and packet.get("rna_context", packet.get("mutation_context")) in {"rna_high", "qualifying_event_recorded"}:
                raw = values.get("raw_dependency_score")
                if isinstance(raw, (int, float)) and math.isfinite(raw):
                    state = "support" if raw <= boundary else "conflict"
            elif family == "crispr_codependency":
                effects = values.get("raw_effects_by_gene", {})
                finite = [value for value in effects.values() if isinstance(value, (int, float)) and math.isfinite(value)]
                if len(finite) == 2:
                    below = sum(value <= boundary for value in finite)
                    state = "support" if below == 2 else "neutral" if below == 1 else "conflict"
        counts[str(state)] += 1
    return {key: counts.get(key, 0) for key in ("support", "neutral", "conflict", "unknown", "not_applicable")}


def summarize_result(control: Mapping[str, Any], result: Mapping[str, Any]) -> dict[str, Any]:
    """Keep complete control-level method outputs while avoiding a second full evidence export."""
    statistic = result.get("qualification_statistic", {})
    packets = result.get("model_packets", [])
    lineage_counts = Counter(str(packet.get("lineage")) for packet in packets if packet.get("lineage"))
    assay = Counter()
    for packet in packets:
        for key, value in packet.items():
            if key.endswith("assay_state"):
                assay[f"{key}:{value}"] += 1
    return {
        "control_id": control["control_id"], "control_kind": control["kind"], "fixture_software_contract": control["kind"] in _FIXTURE_KINDS,
        "relationship": control["relationship"], "route_family": result.get("route_family"), "route_status": result.get("route_status"),
        "route_qualified": result.get("route_qualified"), "relationship_id": result.get("relationship_id"),
        "denominator": result.get("denominator", {}), "qualification_statistic": statistic,
        "model_support_counts": _packet_state_counts(result), "assay_missingness_counts": dict(sorted(assay.items())),
        "missingness": result.get("missingness", {}), "lineage_distribution": dict(sorted(lineage_counts.items())),
        "lineage_assay_bias": {"largest_lineage_models": max(lineage_counts.values(), default=0), "lineages_with_packets": len(lineage_counts), "assay_availability_is_route_packet_reported": True},
        "leave_one_lineage_out": statistic.get("leave_one_lineage_out"),
        "source_ablation": {"state": "not_defined", "reason": "Each implemented route has one required governed assay stream for its primary statistic; replacing or omitting it would define a different route rather than a source ablation."},
        "limitations": result.get("limitations", []),
    }


def _conditional_rows(result: Mapping[str, Any], *, quantiles: tuple[float, float] = (0.2, 0.8), lineage_minimum: int = 15) -> list[dict[str, Any]]:
    """Reconstruct BR04/BR05 rows exactly, including RNA-middle rows for centring."""
    packets = result.get("model_packets", [])
    rna = result.get("relationship", {}).get("context_type") == "rna_high"
    candidate: list[dict[str, Any]] = []
    for packet in packets:
        value = packet.get("route_native_value", {}).get("raw_dependency_score")
        lineage = packet.get("lineage")
        if packet.get("problematic_model_state") != "not_problematic" or not isinstance(lineage, str) or not isinstance(value, (int, float)) or not math.isfinite(value):
            continue
        group = None
        if not rna:
            group = "context" if packet.get("mutation_context") == "qualifying_event_recorded" else "reference" if packet.get("mutation_context") == "no_qualifying_event_recorded" else None
            if group is None:
                continue
        candidate.append({"ModelID": str(packet.get("ModelID")), "lineage": lineage, "raw": float(value), "rna": packet.get("route_native_value", {}).get("raw_rna_log2tpm1"), "group": group})
    by_lineage: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in candidate:
        by_lineage[row["lineage"]].append(row)
    output: list[dict[str, Any]] = []
    for lineage, rows in sorted(by_lineage.items()):
        if rna:
            rows = [row for row in rows if isinstance(row["rna"], (int, float)) and math.isfinite(row["rna"])]
            if len(rows) < lineage_minimum:
                continue
            values = np.asarray(sorted(float(row["rna"]) for row in rows))
            low, high = np.quantile(values, quantiles, method="linear")
            if high <= low:
                continue
            for row in rows:
                row["group"] = "context" if row["rna"] >= high else "reference" if row["rna"] <= low else None
                if row["group"] is None:
                    row["group"] = "middle"
        if not rna and len(rows) < lineage_minimum:
            continue
        median = float(np.median([row["raw"] for row in rows]))
        output.extend({**row, "adjusted": row["raw"] - median} for row in sorted(rows, key=lambda item: item["ModelID"]))
    return output


def _bootstrap_shift(rows: list[dict[str, Any]], replicates: int) -> list[float]:
    strata: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        strata[(row["lineage"], row["group"])].append(row["adjusted"])
    generator = np.random.default_rng(25)
    output = []
    for _ in range(replicates):
        sampled = {"context": [], "reference": []}
        for (_, group), values in sorted(strata.items()):
            values = np.asarray(sorted(values))
            sampled[group].extend(generator.choice(values, size=len(values), replace=True))
        output.append(float(np.median(sampled["reference"]) - np.median(sampled["context"])))
    return output


def _conditional_policy_status(rows: list[dict[str, Any]], result: Mapping[str, Any], *, group_minimum: int, total_minimum: int, lineage_minimum: int, stability_minimum: float, replicates: int, quantiles: tuple[float, float]) -> dict[str, Any]:
    # rows are reconstructed first, then filtered to make every policy change visible and bounded.
    eligible = _conditional_rows(result, quantiles=quantiles, lineage_minimum=lineage_minimum)
    analysis = [row for row in eligible if row["group"] in {"context", "reference"}]
    context = [row for row in analysis if row["group"] == "context"]
    reference = [row for row in analysis if row["group"] == "reference"]
    shift = float(np.median([row["adjusted"] for row in reference]) - np.median([row["adjusted"] for row in context])) if context and reference else None
    lineages = sorted({row["lineage"] for row in eligible})
    gates = len(context) >= group_minimum and len(reference) >= group_minimum and len(analysis) >= total_minimum and len(lineages) >= 2 and shift is not None
    boot = _bootstrap_shift(analysis, replicates) if gates else []
    ci = [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))] if boot else [None, None]
    held = []
    for lineage in lineages:
        retained = [row for row in analysis if row["lineage"] != lineage]
        a, b = [row for row in retained if row["group"] == "context"], [row for row in retained if row["group"] == "reference"]
        value = float(np.median([row["adjusted"] for row in b]) - np.median([row["adjusted"] for row in a])) if a and b else None
        estimable = value is not None and len(a) >= group_minimum and len(b) >= group_minimum and len(retained) >= total_minimum and len({row["lineage"] for row in retained}) >= 2
        held.append((estimable, value))
    estimable = [value for ok, value in held if ok]
    stability = sum((value > 0) == (shift > 0) for value in estimable) / len(estimable) if estimable and shift not in (None, 0) else None
    essential = result.get("common_essential", {}).get("is_common_essential")
    if not gates:
        status = "unavailable"
    elif essential is not False:
        status = "exploratory"
    elif shift > 0 and ci[0] > 0 and stability is not None and stability >= stability_minimum:
        status = "supported"
    elif shift < 0 and ci[1] < 0 and stability is not None and stability >= stability_minimum:
        status = "contradicted"
    else:
        status = "exploratory"
    return {"policy_status": status, "shift": shift, "ci": ci, "context_n": len(context), "reference_n": len(reference), "total_n": len(analysis), "lineages": len(lineages), "loo_stability": stability}


def _codependency_policy_status(result: Mapping[str, Any], *, shared_minimum: int, lineage_minimum: int, stability_minimum: float, replicates: int) -> dict[str, Any]:
    packets = result.get("model_packets", [])
    rows = []
    for packet in packets:
        effects = packet.get("route_native_value", {}).get("raw_effects_by_gene", {})
        lineage = packet.get("lineage")
        values = [effects[key] for key in sorted(effects)]
        if packet.get("problematic_model_state") == "not_problematic" and isinstance(lineage, str) and len(values) == 2 and all(isinstance(value, (int, float)) and math.isfinite(value) for value in values):
            rows.append({"ModelID": str(packet.get("ModelID")), "lineage": lineage, "a": float(values[0]), "b": float(values[1])})
    by_lineage: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_lineage[row["lineage"]].append(row)
    rows = [row for lineage, items in sorted(by_lineage.items()) if len(items) >= lineage_minimum for row in sorted(items, key=lambda item: item["ModelID"])]
    for lineage in sorted({row["lineage"] for row in rows}):
        members = [row for row in rows if row["lineage"] == lineage]
        ma, mb = np.median([row["a"] for row in members]), np.median([row["b"] for row in members])
        for row in members:
            row["aa"], row["bb"] = row["a"] - ma, row["b"] - mb
    def rho(items: Iterable[dict[str, Any]]) -> float | None:
        values = list(items)
        if len(values) < 2 or np.ptp([row["aa"] for row in values]) == 0 or np.ptp([row["bb"] for row in values]) == 0:
            return None
        value = float(spearmanr([row["aa"] for row in values], [row["bb"] for row in values]).statistic)
        return value if math.isfinite(value) else None
    value, lineages = rho(rows), sorted({row["lineage"] for row in rows})
    gates = value is not None and len(rows) >= shared_minimum and len(lineages) >= 2
    generator, boot = np.random.default_rng(25), []
    if gates:
        by_lineage = {lineage: sorted([row for row in rows if row["lineage"] == lineage], key=lambda item: item["ModelID"]) for lineage in lineages}
        for _ in range(replicates):
            sample = [items[i] for items in by_lineage.values() for i in generator.integers(0, len(items), len(items))]
            candidate = rho(sample)
            if candidate is not None:
                boot.append(candidate)
    ci = [float(np.quantile(boot, .025)), float(np.quantile(boot, .975))] if len(boot) == replicates else [None, None]
    full_sign = 1 if value is not None and value > 0 else -1 if value is not None and value < 0 else 0
    loo = []
    for held_out in lineages:
        retained = [row for row in rows if row["lineage"] != held_out]
        retained_lineages = sorted({row["lineage"] for row in retained})
        holdout = rho(retained)
        estimable = len(retained) >= shared_minimum and len(retained_lineages) >= 2 and holdout is not None
        holdout_sign = 1 if holdout is not None and holdout > 0 else -1 if holdout is not None and holdout < 0 else 0
        loo.append({"held_out_lineage": held_out, "adjusted_spearman_rho": holdout if estimable else None, "estimable": estimable, "matches_full_sign": bool(estimable and full_sign != 0 and holdout_sign == full_sign), "denominator": {"shared_evaluable_models": len(retained), "eligible_lineages": len(retained_lineages)}})
    estimable_loo = [row for row in loo if row["estimable"]]
    stable = sum(row["matches_full_sign"] for row in estimable_loo) / len(estimable_loo) if estimable_loo else None
    common = result.get("common_essential", {})
    if not gates:
        status = "unavailable"
    elif common.get("source_is_common_essential") is not False or common.get("target_is_common_essential") is not False:
        status = "exploratory"
    elif (ci[0] > 0 or ci[1] < 0) and stable is not None and stable >= stability_minimum:
        status = "supported"
    else:
        status = "exploratory"
    return {"policy_status": status, "adjusted_rho": value, "ci": ci, "shared_n": len(rows), "lineages": len(lineages), "leave_one_lineage_out": loo, "loo_stability": stable}


def sensitivity_rows(control: Mapping[str, Any], result: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Run only the predeclared, one-factor-at-a-time policy/model sensitivity grid."""
    grid, family, output = fixed_sensitivity_grid(), result.get("route_family"), []
    for boundary in grid["dependency_support_boundary_model_only"]:
        output.append({"control_id": control["control_id"], "route_family": family, "sensitivity_type": "model_only_boundary", "parameter": "dependency_support_boundary", "value": boundary, "route_status": result.get("route_status"), "model_support_counts": _packet_state_counts(result, boundary), "note": "Model-only boundary change; route policy status is retained by design."})
    if family == "conditional_dependency":
        rna = control["relationship"].get("context_type") == "rna_high"
        scenarios = [("group_minimum", value) for value in grid["conditional_context_reference_minimum_n"]] + [("total_minimum", value) for value in grid["conditional_total_minimum_n"]] + [("lineage_minimum", value) for value in grid["lineage_minimum_n"]] + [("loo_stability", value) for value in grid["leave_one_lineage_out_stability_minimum"]] + [("bootstrap_replicates", value) for value in grid["bootstrap_replicates"]]
        if rna:
            scenarios += [("rna_quantiles", value) for value in grid["rna_quantile_pairs"]]
        for parameter, value in scenarios:
            kwargs = {"group_minimum": 10, "total_minimum": 30, "lineage_minimum": 15, "stability_minimum": .8, "replicates": 500, "quantiles": (.2, .8)}
            if parameter == "group_minimum": kwargs["group_minimum"] = int(value)
            elif parameter == "total_minimum": kwargs["total_minimum"] = int(value)
            elif parameter == "lineage_minimum": kwargs["lineage_minimum"] = int(value)
            elif parameter == "loo_stability": kwargs["stability_minimum"] = float(value)
            elif parameter == "bootstrap_replicates": kwargs["replicates"] = int(value)
            elif parameter == "rna_quantiles": kwargs["quantiles"] = tuple(value)
            output.append({"control_id": control["control_id"], "route_family": family, "sensitivity_type": "policy_status_rerun", "parameter": parameter, "value": value, **_conditional_policy_status([], result, **kwargs)})
    elif family == "crispr_codependency":
        scenarios = [("shared_minimum", value) for value in grid["codependency_shared_minimum_n"]] + [("lineage_minimum", value) for value in grid["lineage_minimum_n"]] + [("loo_stability", value) for value in grid["leave_one_lineage_out_stability_minimum"]] + [("bootstrap_replicates", value) for value in grid["bootstrap_replicates"]]
        for parameter, value in scenarios:
            kwargs = {"shared_minimum": 30, "lineage_minimum": 15, "stability_minimum": .8, "replicates": 500}
            if parameter == "shared_minimum": kwargs["shared_minimum"] = int(value)
            elif parameter == "lineage_minimum": kwargs["lineage_minimum"] = int(value)
            elif parameter == "loo_stability": kwargs["stability_minimum"] = float(value)
            else: kwargs["replicates"] = int(value)
            output.append({"control_id": control["control_id"], "route_family": family, "sensitivity_type": "policy_status_rerun", "parameter": parameter, "value": value, **_codependency_policy_status(result, **kwargs)})
    return output


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_bytes(to_json_safe(value)) + b"\n")


def _write_csv(path: Path, rows: list[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(to_json_safe(value), sort_keys=True, separators=(",", ":")) if isinstance(value, (dict, list)) else value for key, value in row.items()})


def evaluate_control(control_id: str, output_dir: Path = EVALUATION_DIRECTORY, *, snapshot_id: str = DEFAULT_SNAPSHOT_ID, data_dir: Path = DEFAULT_DATA_DIRECTORY) -> dict[str, Any]:
    """Evaluate one frozen control and persist only its deterministic BR09 summary.

    Splitting real controls is operational only: it lets a long governed Parquet read resume after
    an interrupted desktop session.  It cannot alter the manifest, grid, policy, or control list.
    """
    manifest = verify_execution_manifest(output_dir, snapshot_id=snapshot_id, data_dir=data_dir)
    manifest_sha256 = _sha256_file(output_dir / "execution_manifest.json")
    if manifest["analytic_identity"]["snapshot_id"] != snapshot_id:
        raise RuntimeError("requested snapshot differs from frozen BR09 execution manifest")
    snapshot_path = DEFAULT_SNAPSHOT_DIRECTORY / f"{snapshot_id}.json"
    if _sha256_file(snapshot_path) != manifest["analytic_identity"]["snapshot_sha256"]:
        raise RuntimeError("frozen snapshot bytes changed after manifest creation")
    assert_validation_import_boundary()
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    controls = {row["control_id"]: row for row in load_controls()["controls"]}
    if control_id not in controls:
        raise ValueError(f"control_id is not frozen: {control_id}")
    control = controls[control_id]
    if control["kind"] in _FIXTURE_KINDS:
        summary, sensitivity = {"control_id": control["control_id"], "control_kind": control["kind"], "fixture_software_contract": True, "execution_state": "exercised_by_deterministic_fixture_tests", "relationship": control["relationship"], "route_status": None, "route_qualified": None, "limitations": ["Fixture/software-contract control only; it is not reported as real biology."]}, []
        _write_json(output_dir / "full_control_results" / f"{control_id}.json", {"execution_manifest_sha256": manifest_sha256, "control": control, "execution_state": "fixture_software_contract_only", "result": None})
    else:
        relationship = control["relationship"]
        adapter = adapt_route_evidence(snapshot, relationship["route_family"], relationship["source_gene"], relationship["target_gene"], relationship.get("context_type"), data_dir=data_dir, cache_dir=output_dir / "adapter_cache")
        if relationship["route_family"] == "direct_fusion": result = evaluate_direct_fusion(adapter)
        elif relationship["route_family"] == "conditional_dependency": result = evaluate_conditional_dependency(adapter)
        else: result = evaluate_codependency(adapter)
        summary = summarize_result(control, result)
        summary["execution_state"] = "real_governed_processed_data_internal_evaluation"
        sensitivity = sensitivity_rows(control, result)
        _write_json(output_dir / "full_control_results" / f"{control_id}.json", {"execution_manifest_sha256": manifest_sha256, "control": control, "result": result})
    _write_json(output_dir / "control_parts" / f"{control_id}.json", {"execution_manifest_sha256": manifest_sha256, "summary": summary, "sensitivity": sensitivity})
    return {"summary": summary, "sensitivity": sensitivity}


def compile_evaluation(output_dir: Path = EVALUATION_DIRECTORY) -> dict[str, Any]:
    """Compile every frozen control part in registry order without rerunning an engine."""
    manifest = verify_execution_manifest(output_dir)
    manifest_sha256 = _sha256_file(output_dir / "execution_manifest.json")
    controls = load_controls()["controls"]
    summaries, sensitivity = [], []
    for control in controls:
        part = output_dir / "control_parts" / f"{control['control_id']}.json"
        if not part.exists():
            raise RuntimeError(f"frozen control has not been evaluated: {control['control_id']}")
        value = json.loads(part.read_text(encoding="utf-8"))
        if value.get("execution_manifest_sha256") != manifest_sha256:
            raise RuntimeError(f"control result was not produced under the frozen final manifest: {control['control_id']}")
        summaries.append(value["summary"])
        sensitivity.extend(value["sensitivity"])
    _write_json(output_dir / "control_results.json", {"schema_version": VALIDATION_SCHEMA_VERSION, "execution_manifest_sha256": manifest_sha256, "analytic_identity": manifest["analytic_identity"], "controls": summaries})
    _write_json(output_dir / "source_availability.json", {"execution_manifest_sha256": manifest_sha256, **manifest["source_availability"]})
    _write_csv(output_dir / "sensitivity_results.csv", sensitivity)
    _write_json(output_dir / "sensitivity_results.json", {"schema_version": VALIDATION_SCHEMA_VERSION, "execution_manifest_sha256": manifest_sha256, "rows": sensitivity})
    counts = []
    for summary in summaries:
        for state, count in summary.get("model_support_counts", {}).items():
            counts.append({"control_id": summary["control_id"], "state": state, "count": count})
    _write_csv(output_dir / "model_support_counts.csv", counts)
    by_control = {summary["control_id"]: summary for summary in summaries}
    source = manifest["source_availability"]
    verdicts = {
        "direct_fusion": {"verdict": "adopt", "software_implemented": True, "internally_supported": by_control["direct_fusion_bcr_abl1_positive"]["route_status"] == "supported", "internally_unsupported_control_retained": by_control["direct_fusion_bcr_muc1_unrelated"]["route_status"] == "contradicted", "independently_corroborated": False, "independent_gap": source["direct_fusion_independent_corroboration"]["reason"]},
        "conditional_dependency_mutation": {"verdict": "adopt", "external_corroboration": "deferred", "software_implemented": True, "internally_supported_controls": [item for item in ("conditional_brca1_parp1_mutation", "conditional_braf_map2k1_mutation") if by_control[item]["route_status"] == "supported"], "independently_corroborated": False, "independent_gap": source["sanger_project_score"]["reason"]},
        "conditional_dependency_rna_high": {"verdict": "adopt", "external_corroboration": "deferred", "software_implemented": True, "internally_supported_controls": [item for item in ("conditional_erbb2_erbb3_rna_high", "conditional_myc_max_rna_high") if by_control[item]["route_status"] == "supported"], "independently_corroborated": False, "independent_gap": "PRISM is unavailable/not interpretable because its governed metadata and response files are LFS pointers and no gene-to-compound mapping is available."},
        "crispr_codependency": {"verdict": "defer", "presentation": "exploratory", "software_implemented": True, "internally_supported": by_control["codependency_braf_map2k1"]["route_status"] == "supported", "independently_corroborated": False, "independent_gap": "No suitable governed independent corroboration was available; the family remains an exploratory signed association."},
    }
    if any(value["verdict"] not in {"adopt", "revise", "defer"} for value in verdicts.values()):
        raise RuntimeError("BR09 verdict must use the controlled adopt/revise/defer vocabulary")
    _write_json(output_dir / "route_family_verdicts.json", {"schema_version": VALIDATION_SCHEMA_VERSION, "execution_manifest_sha256": manifest_sha256, "verdicts": verdicts})
    return {"controls": summaries, "sensitivity": sensitivity, "manifest": manifest}


def run_evaluation(output_dir: Path = EVALUATION_DIRECTORY, *, snapshot_id: str = DEFAULT_SNAPSHOT_ID, data_dir: Path = DEFAULT_DATA_DIRECTORY) -> dict[str, Any]:
    """Convenience complete runner; ``evaluate_control`` is the resumable equivalent."""
    for control in load_controls()["controls"]:
        evaluate_control(control["control_id"], output_dir, snapshot_id=snapshot_id, data_dir=data_dir)
    return compile_evaluation(output_dir)
