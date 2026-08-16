"""BR02 read-only projection of a frozen research snapshot into route evidence.

This module deliberately stops before route evaluation.  It resolves only exact frozen
``gene_reference`` entries, reads the complete saved snapshot partition (not its displayed
top-N), and makes long evidence rows without changing native evidence, D, rank, tier or veto.
The existing public ``backend.research_service.adapter`` remains the bounded measurement reader;
this module adds only the route-specific model spine, gene-reference and event-pair handling.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd

_RESEARCH_SERVICE_DIR = str(Path(__file__).resolve().parents[1])
if _RESEARCH_SERVICE_DIR not in sys.path:
    sys.path.insert(0, _RESEARCH_SERVICE_DIR)

import adapter  # noqa: E402 -- sys.path must be set up first
from serialization import canonical_bytes, content_hash, to_json_safe  # noqa: E402
from .contracts import (
    ContextType,
    EvidenceLedgerRecord,
    EvidenceState,
    Relationship,
    ResolvedGene,
    RouteFamily,
    validate_ledger,
)
from .registry import load_route_registry, policy_bytes, route_input_sources


_PARTITION_BUCKETS = (
    "ranked_cell_lines",
    "ranked_beyond_top_n",
    "low_confidence_lines",
    "insufficient_evidence_lines",
    "disqualified_lines",
)
_DEFAULT_CACHE_DIR = Path(__file__).resolve().parents[1] / "runtime" / "routes"
_DEFAULT_COMMON_ESSENTIAL = Path(__file__).resolve().parents[3] / "scoring" / "resources" / "common_essential_genes.json"
ADAPTER_SCHEMA_VERSION = "br02-evidence-adapter-v2"
_LAYER_DETAILS = {
    "rna": ("processed_expression_rna", "expression_rna_state", "log2tpm1", "log2(TPM+1)"),
    "dependency": ("processed_dependency", "dependency_state", "dependency_score", "Chronos dependency score"),
    "mutations": ("processed_mutations", "mutations_state", None, "processed mutation event"),
    "fusions": ("processed_fusions", "fusions_state", None, "processed fusion event"),
}


def _resolve_data_dir(data_dir: str | Path | None) -> Path:
    """Use the extension adapter's documented processed-data boundary without a private import."""
    if data_dir is not None:
        return Path(data_dir)
    configured = os.environ.get("SCORING_SERVICE_DATA_DIR")
    if configured:
        return Path(configured)
    return Path(__file__).resolve().parents[3] / "data" / "processed"


@dataclass(frozen=True)
class GeneResolution:
    """Result of exact frozen-reference lookup, including the explicit failure state."""

    submitted: str
    state: EvidenceState
    gene: ResolvedGene | None
    candidates: tuple[ResolvedGene, ...] = ()

    @property
    def resolved(self) -> bool:
        return self.state is EvidenceState.OBSERVED and self.gene is not None

    def as_dict(self) -> dict[str, Any]:
        """Return strict JSON data without losing ambiguous frozen-reference candidates."""
        return {
            "submitted": self.submitted,
            "state": self.state.value,
            "gene": self.gene.as_dict() if self.gene else None,
            "candidates": [candidate.as_dict() for candidate in self.candidates],
        }


def _processed_path(data_dir: Path, filename: str) -> Path:
    """Prefer the immutable Parquet mirror, matching the public adapter's read convention."""
    parquet = data_dir / filename.replace(".csv", ".parquet")
    return parquet if parquet.exists() else data_dir / filename


def _read_table(path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    """Read a small metadata table or an already bounded event subset from frozen outputs."""
    if not path.exists():
        raise FileNotFoundError(f"No processed table at {path}")
    if path.suffix == ".parquet":
        return pd.read_parquet(path, columns=columns)
    return pd.read_csv(path, usecols=columns)


def _gene_reference(data_dir: Path) -> pd.DataFrame:
    frame = _read_table(_processed_path(data_dir, "gene_reference.csv"), ["ensembl_id", "symbol"])
    if frame["ensembl_id"].isna().any() or frame["symbol"].isna().any():
        raise ValueError("processed gene_reference has null canonical identifiers")
    return frame.astype({"ensembl_id": "string", "symbol": "string"})


def resolve_gene(submitted: str, data_dir: str | Path | None = None) -> GeneResolution:
    """Resolve an exact symbol or exact Ensembl ID with no aliases, normalization or guessing.

    A duplicate symbol is intentionally returned as ``unresolvable`` with all candidates rather
    than choosing one.  Callers can surface that state without converting it to missing assay
    evidence or silently selecting an isoform-like match.
    """
    if not isinstance(submitted, str) or not submitted.strip():
        raise ValueError("submitted gene must be a non-empty string")
    value = submitted.strip()
    base_dir = _resolve_data_dir(data_dir)
    reference = _gene_reference(base_dir)
    matches = reference.loc[(reference["symbol"] == value) | (reference["ensembl_id"] == value)]
    candidates = tuple(
        ResolvedGene(symbol=str(row.symbol), ensembl_id=str(row.ensembl_id))
        for row in matches.drop_duplicates(["symbol", "ensembl_id"]).sort_values(["ensembl_id", "symbol"]).itertuples(index=False)
    )
    unique_ids = {candidate.ensembl_id for candidate in candidates}
    if len(unique_ids) != 1:
        return GeneResolution(value, EvidenceState.UNRESOLVABLE, None, candidates)
    selected = candidates[0]
    return GeneResolution(value, EvidenceState.OBSERVED, selected, candidates)


def evidence_state_for(coverage_state: Any, has_row: bool, applicable: bool = True) -> EvidenceState:
    """Translate frozen coverage plus row presence without inventing a biological negative."""
    if not applicable:
        return EvidenceState.NOT_APPLICABLE
    if has_row:
        return EvidenceState.OBSERVED
    if coverage_state == "not_assayed" or coverage_state is None or pd.isna(coverage_state):
        return EvidenceState.NOT_ASSAYED
    # Both explicit measured_absent and a measured model with no requested-gene row mean no
    # retained value for this query; neither is coerced to a numeric zero.
    return EvidenceState.MEASURED_ABSENT


def _snapshot_models(snapshot: Mapping[str, Any]) -> list[str]:
    """Enumerate every immutable partition bucket in fixed bucket/model order."""
    models: list[str] = []
    seen: set[str] = set()
    native = snapshot.get("native_result") or {}
    for bucket in _PARTITION_BUCKETS:
        for row in native.get(bucket, []) or []:
            model_id = row.get("model_id") or row.get("ModelID")
            if isinstance(model_id, str) and model_id and model_id not in seen:
                seen.add(model_id)
                models.append(model_id)
    if not models:
        for row in snapshot.get("models", []) or []:
            model_id = row.get("ModelID")
            if isinstance(model_id, str) and model_id and model_id not in seen:
                seen.add(model_id)
                models.append(model_id)
    if not models:
        raise ValueError("snapshot has no immutable native model partition")
    return models


def _model_context(model_ids: list[str], data_dir: Path) -> pd.DataFrame:
    """Join only requested model metadata and coverage, retaining the snapshot's model order."""
    metadata, unknown = adapter.read_model_metadata(model_ids, data_dir=data_dir)
    if unknown:
        raise ValueError(f"snapshot ModelID values absent from frozen cell_lines: {unknown}")
    coverage_path = _processed_path(data_dir, "coverage.csv")
    coverage = _read_table(coverage_path)
    coverage = coverage.loc[coverage["ModelID"].isin(model_ids)].copy()
    if coverage["ModelID"].duplicated().any():
        raise ValueError("processed coverage has duplicate ModelID rows")
    result = pd.DataFrame({"ModelID": model_ids}).merge(metadata, on="ModelID", how="left", validate="one_to_one")
    result = result.merge(coverage, on="ModelID", how="left", validate="one_to_one")
    return result


def _event_rows(data_dir: Path, filename: str, model_ids: Iterable[str], gene_ids: Iterable[str]) -> pd.DataFrame:
    """Read event tables for the requested model/gene set while preserving every event row.

    Public ``adapter.read_measurements`` already handles single gene-key event reads.  Fusion
    pairs additionally need ``partner_ensembl_id`` so a direct exact query cannot lose calls
    where the requested gene is stored on the partner side.
    """
    path = _processed_path(data_dir, filename)
    models, genes = list(dict.fromkeys(model_ids)), list(dict.fromkeys(gene_ids))
    if not models or not genes:
        return pd.DataFrame()

    if path.suffix == ".parquet":
        if filename == "fusions.csv":
            # Parquet filters accept disjunctive normal form: requested model AND either side
            # of the pair. This avoids materialising the full fusion table for a two-gene query.
            filters = [
                [("ModelID", "in", models), ("ensembl_id", "in", genes)],
                [("ModelID", "in", models), ("partner_ensembl_id", "in", genes)],
            ]
        else:
            filters = [("ModelID", "in", models), ("ensembl_id", "in", genes)]
        return pd.read_parquet(path, filters=filters).reset_index(drop=True)

    if not path.exists():
        raise FileNotFoundError(f"No processed table at {path}")
    matches = []
    for chunk in pd.read_csv(path, chunksize=1_000_000):
        mask = chunk["ModelID"].isin(models) & chunk["ensembl_id"].isin(genes)
        if filename == "fusions.csv":
            mask = chunk["ModelID"].isin(models) & (
                chunk["ensembl_id"].isin(genes) | chunk["partner_ensembl_id"].isin(genes)
            )
        matches.append(chunk.loc[mask])
    return pd.concat(matches, ignore_index=True) if matches else pd.DataFrame()


def _json_record(row: Mapping[str, Any]) -> dict[str, Any]:
    """Preserve source columns/raw values while making pandas nulls valid strict JSON nulls."""
    return to_json_safe(dict(row))


def _source_stamp(paths: Iterable[Path]) -> str:
    """Hash immutable file identities, not full multi-gigabyte frozen tables at query time."""
    entries = []
    for path in sorted(paths, key=lambda item: str(item)):
        if not path.exists():
            raise FileNotFoundError(f"required frozen route source is absent: {path}")
        stat = path.stat()
        entries.append({"name": path.name, "size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    return content_hash(entries)


def _common_essential(path: Path, genes: Iterable[ResolvedGene]) -> dict[str, dict[str, Any]]:
    """Attach per-gene selectivity warnings, never an aggregate score or route conclusion."""
    if not path.exists():
        raise FileNotFoundError(f"common-essential policy resource is absent: {path}")
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    pan_essential = payload.get("pan_essential", {})
    return {
        gene.ensembl_id: {
            "source_id": "common_essential_gene_policy",
            "source_version": str(payload.get("version", "frozen-v7-resource")),
            "symbol": gene.symbol,
            "ensembl_id": gene.ensembl_id,
            "is_common_essential": gene.ensembl_id in pan_essential,
            "fraction": pan_essential.get(gene.ensembl_id),
            "limitation": "A common-essential flag is a selectivity warning, not route evidence or a score.",
        }
        for gene in sorted(genes, key=lambda item: item.ensembl_id)
    }


def _ledger_record(
    query_id: str,
    relationship: Relationship,
    policy_version: str,
    source_id: str,
    source_version: str,
    data_id: str,
    model_id: str,
    gene_id: str,
    raw_value: Any,
    unit: str,
    state: EvidenceState,
    source_record_id: str | None,
    denominator: Mapping[str, int],
) -> EvidenceLedgerRecord:
    return EvidenceLedgerRecord(
        research_query_id=query_id,
        route_family=relationship.route_family,
        model_id=model_id,
        ensembl_id=gene_id,
        relationship=relationship,
        route_policy_version=policy_version,
        source_id=source_id,
        source_version=source_version,
        data_id=data_id,
        raw_value=raw_value,
        unit=unit,
        transformation="none; frozen processed value retained",
        denominator=denominator,
        evidence_state=state,
        qc_state="not_evaluated_by_BR02",
        limitation="BR02 is a read-only evidence adapter; it does not qualify, score or rank a route.",
        source_record_id=source_record_id,
    )


def _measurement_records(
    query_id: str,
    relationship: Relationship,
    policy_version: str,
    context: pd.DataFrame,
    frames: Mapping[str, pd.DataFrame],
    registry: Mapping[str, Any],
) -> list[EvidenceLedgerRecord]:
    """Build one deterministic long row per retained value or per explicit missing state."""
    registrations = {item["source_id"]: item for item in registry["sources"]}
    records: list[EvidenceLedgerRecord] = []
    route_layers = {
        RouteFamily.DIRECT_FUSION: ("fusions",),
        RouteFamily.CONDITIONAL_DEPENDENCY: ("rna", "dependency", "mutations"),
        RouteFamily.CRISPR_CODEPENDENCY: ("dependency",),
    }[relationship.route_family]
    genes_for_layer = {
        "rna": (relationship.source_gene,),
        "dependency": (relationship.target_gene,) if relationship.route_family is RouteFamily.CONDITIONAL_DEPENDENCY else (relationship.source_gene, relationship.target_gene),
        "mutations": (relationship.source_gene,),
        "fusions": (relationship.source_gene,),
    }
    for layer in route_layers:
        source_id, coverage_column, value_column, unit = _LAYER_DETAILS[layer]
        registration = registrations[source_id]
        frame = frames[layer]
        denominator = {
            "snapshot_models": len(context),
            "layer_measured_models": int(context[coverage_column].isin(["measured", "measured_absent"]).sum()),
        }
        for gene in genes_for_layer[layer]:
            for model in context.itertuples(index=False):
                model_id = str(model.ModelID)
                coverage_state = getattr(model, coverage_column, None)
                matches = frame.loc[(frame["ModelID"] == model_id) & (frame["ensembl_id"] == gene.ensembl_id)]
                if layer == "fusions" and "partner_ensembl_id" in frame:
                    matches = frame.loc[
                        (frame["ModelID"] == model_id)
                        & ((frame["ensembl_id"] == relationship.source_gene.ensembl_id) | (frame["partner_ensembl_id"] == relationship.source_gene.ensembl_id))
                    ]
                if matches.empty:
                    records.append(_ledger_record(
                        query_id, relationship, policy_version, source_id, registration["source_version"],
                        f"{source_id}:{model_id}:{gene.ensembl_id}", model_id, gene.ensembl_id,
                        None, unit, evidence_state_for(coverage_state, False), None, denominator,
                    ))
                    continue
                # Source row position is not identity. Sort strict source content and ordinalise
                # exact duplicates so a CSV/Parquet representation or row reordering is stable.
                ordered_events = []
                for event in matches.to_dict("records"):
                    raw = event[value_column] if value_column else _json_record(event)
                    event_hash = hashlib.sha256(canonical_bytes(_json_record(event))).hexdigest()
                    ordered_events.append((event_hash, canonical_bytes(_json_record(event)), raw))
                occurrence_by_hash: dict[str, int] = {}
                for event_hash, _, raw in sorted(ordered_events, key=lambda item: (item[0], item[1])):
                    occurrence = occurrence_by_hash.get(event_hash, 0)
                    occurrence_by_hash[event_hash] = occurrence + 1
                    source_record_id = f"{event_hash}:{occurrence}"
                    records.append(_ledger_record(
                        query_id, relationship, policy_version, source_id, registration["source_version"],
                        f"{source_id}:{source_record_id}", model_id, gene.ensembl_id,
                        raw, unit, EvidenceState.OBSERVED, source_record_id, denominator,
                    ))
    return records


def _cache_path(cache_dir: Path, query_hash: str, data_hash: str, policy_hash: str) -> Path:
    key = hashlib.sha256(f"{query_hash}:{data_hash}:{policy_hash}".encode("utf-8")).hexdigest()
    return cache_dir / f"{key}.json"


def adapt_route_evidence(
    snapshot: Mapping[str, Any],
    route_family: RouteFamily | str,
    source_gene: str,
    target_gene: str,
    context_type: ContextType | str | None = None,
    *,
    data_dir: str | Path | None = None,
    cache_dir: str | Path | None = None,
    common_essential_path: str | Path | None = None,
) -> dict[str, Any]:
    """Return a deterministic, cached BR02 evidence projection for one declared relationship.

    Unresolved/ambiguous identity returns an explicit ``unresolvable`` adapter result and no
    fabricated relationship/ledger.  A resolved call makes no biological determination; later
    BR03--BR06 route engines consume the preserved rows.
    """
    original_snapshot_bytes = canonical_bytes(snapshot)
    family = RouteFamily(route_family)
    base_dir = _resolve_data_dir(data_dir)
    source_resolution = resolve_gene(source_gene, base_dir)
    target_resolution = resolve_gene(target_gene, base_dir)
    query_id = snapshot.get("query_id")
    if not isinstance(query_id, str) or not query_id:
        raise ValueError("snapshot requires a non-empty query_id")
    if not source_resolution.resolved or not target_resolution.resolved:
        return {
            "schema_version": ADAPTER_SCHEMA_VERSION,
            "research_query_id": query_id,
            "route_family": family.value,
            "evidence_state": EvidenceState.UNRESOLVABLE.value,
            "gene_resolution": {"source": source_resolution.as_dict(), "target": target_resolution.as_dict()},
            "ledger": [],
            "limitation": "An unresolved or ambiguous frozen gene-reference entry is not guessed.",
        }

    relationship = Relationship(
        route_family=family,
        source_gene=source_resolution.gene,
        target_gene=target_resolution.gene,
        context_type=ContextType(context_type) if context_type is not None else None,
    )
    registry = load_route_registry()
    policy_version = registry["policy_version"]
    # This call is a hard policy gate and documents exactly which registered sources are consumed.
    registered_sources = route_input_sources(family, registry)
    model_ids = _snapshot_models(snapshot)
    context = _model_context(model_ids, base_dir)
    source_files = [
        _processed_path(base_dir, "cell_lines.csv"), _processed_path(base_dir, "gene_reference.csv"),
        _processed_path(base_dir, "coverage.csv"),
    ]
    layers = {
        RouteFamily.DIRECT_FUSION: ("fusions",),
        RouteFamily.CONDITIONAL_DEPENDENCY: ("rna", "dependency", "mutations"),
        RouteFamily.CRISPR_CODEPENDENCY: ("dependency",),
    }[family]
    source_files.extend(_processed_path(base_dir, adapter.LAYER_FILES[layer]) for layer in layers)
    common_essential = None
    if family in {RouteFamily.CONDITIONAL_DEPENDENCY, RouteFamily.CRISPR_CODEPENDENCY}:
        common_essential_file = Path(common_essential_path) if common_essential_path else _DEFAULT_COMMON_ESSENTIAL
        source_files.append(common_essential_file)
        essential_genes = (
            (relationship.target_gene,)
            if family is RouteFamily.CONDITIONAL_DEPENDENCY
            else (relationship.source_gene, relationship.target_gene)
        )
        common_essential = _common_essential(common_essential_file, essential_genes)
    query_hash = content_hash({"query_id": query_id, "relationship": relationship.as_dict(), "models": model_ids})
    data_hash = _source_stamp(source_files)
    policy_hash = hashlib.sha256(policy_bytes()).hexdigest()
    selected_cache_dir = Path(cache_dir) if cache_dir is not None else _DEFAULT_CACHE_DIR
    cache_path = _cache_path(selected_cache_dir, query_hash, data_hash, policy_hash)
    if cache_path.exists():
        with cache_path.open(encoding="utf-8") as handle:
            cached = json.load(handle)
        if cached.get("schema_version") == ADAPTER_SCHEMA_VERSION:
            cached["cache"] = {"status": "hit", "path": str(cache_path), "query_hash": query_hash, "data_hash": data_hash, "policy_hash": policy_hash}
            if canonical_bytes(snapshot) != original_snapshot_bytes:
                raise RuntimeError("BR02 adapter attempted to mutate the immutable native snapshot")
            return cached

    # The public adapter performs bounded Parquet/CSV reads for all normal gene-keyed layers.
    frames: dict[str, pd.DataFrame] = {}
    normal_layers = [layer for layer in layers if layer != "fusions"]
    if normal_layers:
        frames.update(adapter.read_measurements(model_ids, [relationship.source_gene.ensembl_id, relationship.target_gene.ensembl_id], normal_layers, data_dir=base_dir))
    if "fusions" in layers:
        frames["fusions"] = _event_rows(base_dir, "fusions.csv", model_ids, [relationship.source_gene.ensembl_id, relationship.target_gene.ensembl_id])
    records = _measurement_records(query_id, relationship, policy_version, context, frames, registry)
    validate_ledger(records)
    ledger = sorted((record.as_dict() for record in records), key=lambda row: (row["ModelID"], row["ensembl_id"], row["source_id"], row["evidence_id"]))
    result = {
        "schema_version": ADAPTER_SCHEMA_VERSION,
        "research_query_id": query_id,
        "route_family": family.value,
        "route_policy_version": policy_version,
        "relationship": relationship.as_dict(),
        "model_spine": [
            {"ModelID": row.ModelID, "lineage": getattr(row, "lineage", None), "is_problematic": getattr(row, "is_problematic", None)}
            for row in context.itertuples(index=False)
        ],
        "gene_resolution": {"source": source_resolution.as_dict(), "target": target_resolution.as_dict()},
        "route_input_sources": [source.source_id for source in registered_sources],
        "common_essential": common_essential,
        "ledger": ledger,
        "compatibility_manifest": {
            "snapshot_partition": "complete immutable native partition; displayed top-N is not used as a cohort",
            "gene_resolution": "exact symbol or exact Ensembl match in frozen gene_reference only; aliases are not guessed",
            "copy_number": "relative copy number is disclosure-only and no absolute-copy-number threshold is applied",
            "dependency_qc": "BR02 retains frozen dependency_score values and coverage state; no route-specific dependency QC flag is inferred from absent columns",
            "native_behavior": "native D, rank, tier, veto, filters and partitions are read-only and not recalculated",
            "event_multiplicity": "each retained mutation/fusion source row becomes a separate ledger record",
        },
        "cache": {"status": "miss", "path": str(cache_path), "query_hash": query_hash, "data_hash": data_hash, "policy_hash": policy_hash},
    }
    selected_cache_dir.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(to_json_safe(result), handle, sort_keys=True, separators=(",", ":"), allow_nan=False)
        handle.write("\n")
    if canonical_bytes(snapshot) != original_snapshot_bytes:
        raise RuntimeError("BR02 adapter attempted to mutate the immutable native snapshot")
    return result
