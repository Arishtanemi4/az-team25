"""Local orchestration for the optional S09 research API.

This module is the only bridge from a completed native ranking to the extension stages.  It
persists a trusted, deep-copied native export and then lazily attaches additive research packets.
It never accepts a score payload from an HTTP client and never calls the native scorer itself.
"""

import copy
import hashlib
import threading
from pathlib import Path

import adapter
import alternatives
import compare
import context
import graph_context
import profiles
import research_export
import snapshot
from backend.research_service.routes import router as route_router
from backend.research_service.routes import serialization as route_serialization
from backend.research_service.routes.automatic_plan import (
    create_automatic_plan,
    plan_identity,
    validate_persisted_plan,
)
from backend.research_service.routes.contracts import Relationship, ResolvedGene
from backend.research_service.routes.registry import automatic_policy_bytes, load_controls, load_route_registry


class ResearchNotFoundError(LookupError):
    """A syntactically valid query ID has no locally persisted research snapshot."""


class ResearchStateError(RuntimeError):
    """The requested derived operation needs a prior, explicit research selection."""


_PARTITION_BUCKETS = (
    "ranked_cell_lines",
    "ranked_beyond_top_n",
    "low_confidence_lines",
    "insufficient_evidence_lines",
    "disqualified_lines",
)
_AUTOMATIC_EXECUTION_SCHEMA_VERSION = "automatic-route-execution-v1"
_AUTOMATIC_ERROR_LIMIT = 240


def _records(metadata):
    """Make a detached JSON-shaped model spine without importing pandas in this layer."""
    if hasattr(metadata, "to_dict"):
        return copy.deepcopy(metadata.to_dict("records"))
    if isinstance(metadata, list):
        return copy.deepcopy(metadata)
    raise ValueError("Ranking callback model metadata must be a table or list of records")


def _canonical_query(native_result):
    """Derive canonical C1 query identifiers from native per-gene evidence.

    S02 preserves the original user tokens in ``native_result.query``.  C5 must instead pin the
    resolved Ensembl IDs.  We derive only that extension-owned view, leaving native evidence
    untouched and retaining native filters verbatim.
    """
    seen, genes = {"inclusion": set(), "exclusion": set()}, {"inclusion": [], "exclusion": []}
    for bucket in _PARTITION_BUCKETS:
        for line in native_result.get(bucket, []) or []:
            for gene in line.get("per_gene", []) or []:
                role = gene.get("role")
                ensembl_id = gene.get("ensembl_id")
                if role in genes and isinstance(ensembl_id, str) and ensembl_id not in seen[role]:
                    seen[role].add(ensembl_id)
                    genes[role].append(ensembl_id)
    native_query = native_result.get("query") or {}
    return {
        "inclusion_genes": genes["inclusion"] or copy.deepcopy(native_query.get("inclusion_genes", [])),
        "exclusion_genes": genes["exclusion"] or copy.deepcopy(native_query.get("exclusion_genes", [])),
        "filters": copy.deepcopy(native_query.get("filters", {})),
        "identifier_basis": "native_per_gene.ensembl_id" if any(genes.values()) else "native_query",
    }


def _native_model_ids(native_result):
    return {
        line.get("model_id")
        for bucket in _PARTITION_BUCKETS
        for line in native_result.get(bucket, []) or []
        if line.get("model_id")
    }


class ResearchService:
    """Snapshot-backed operations exposed through the opt-in local research router."""

    def __init__(self, snapshot_dir=None, data_dir=None, adapter_module=adapter,
                 route_adapter_function=None, route_evaluators=None):
        self.snapshot_dir = Path(snapshot_dir) if snapshot_dir is not None else None
        self.data_dir = data_dir
        self.adapter = adapter_module
        self.route_adapter_function = route_adapter_function
        self.route_evaluators = route_evaluators
        self._automatic_lock = threading.RLock()

    def capture_native_result(self, native_result, model_metadata):
        """Create the trusted S03 snapshot after native export has completely succeeded."""
        manifest = {
            "baseline_branch": "v7/validate",
            "native_export_schema": native_result.get("schema_version"),
            "extension_method": "research-api-v1",
            "model_metadata_source": "ranking_service.cell_lines_df",
        }
        record = snapshot.create_snapshot(native_result, manifest, _records(model_metadata))
        record["query"] = _canonical_query(native_result)
        # S08 export must reproduce the method used when alternatives are later attached.
        record["method_policy"] = profiles.load_similarity_policy()
        snapshot.save_snapshot(record, self.snapshot_dir)
        return {"research_query_id": record["query_id"], "research_status": "ready"}

    def _load(self, query_id):
        # Validate at this boundary as well as in snapshot.load_snapshot.  This keeps the
        # traversal guarantee true even when a test or future storage implementation replaces
        # the loader behind this service.
        snapshot._validate_query_id(query_id)
        try:
            return snapshot.load_snapshot(query_id, self.snapshot_dir)
        except snapshot.SnapshotNotFoundError as exc:
            raise ResearchNotFoundError("Research result not found; run ranking again to create it.") from exc

    def _save(self, record):
        snapshot.save_snapshot(record, self.snapshot_dir)

    @staticmethod
    def _require_model(record, model_id):
        if model_id not in _native_model_ids(record.get("native_result", {})):
            raise ValueError(f"Unknown model_id for this research query: {model_id!r}")

    @staticmethod
    def _gene_roles(record):
        query = record.get("query") or {}
        ids, roles = [], {}
        for role, field in (("inclusion", "inclusion_genes"), ("exclusion", "exclusion_genes")):
            for gene_id in query.get(field, []) or []:
                if gene_id not in roles:
                    ids.append(gene_id)
                    roles[gene_id] = role
        return ids, roles

    @staticmethod
    def _model_metadata(record, model_id):
        """Return the one snapshot-pinned metadata row, never a name-based lookup or guess."""
        for row in record.get("models", []) or []:
            if isinstance(row, dict) and row.get("ModelID") == model_id:
                return copy.deepcopy(row)
        return None

    def get_context(self, query_id, model_id):
        record = self._load(query_id)
        self._require_model(record, model_id)
        attachments = record.get("context") if isinstance(record.get("context"), dict) else {}
        measured = attachments.setdefault("measured_by_model", {})
        changed = False
        if model_id not in measured:
            measured[model_id] = context.build_measured_context(
                record, model_id, self.adapter, data_dir=self.data_dir
            )
            changed = True
        if "graph_context" not in attachments:
            gene_ids, roles = self._gene_roles(record)
            attachments["graph_context"] = graph_context.build_graph_context(gene_ids, roles)
            changed = True
        record["context"] = attachments
        if changed:
            self._save(record)
        return {
            "query_id": record["query_id"],
            "model_id": model_id,
            "model_metadata": self._model_metadata(record, model_id),
            "measured_context": copy.deepcopy(measured[model_id]),
            "graph_context": copy.deepcopy(attachments["graph_context"]),
        }

    def get_alternatives(self, query_id, model_id):
        record = self._load(query_id)
        self._require_model(record, model_id)
        policy = copy.deepcopy(record.get("method_policy") or profiles.load_similarity_policy())
        fixed_profiles = profiles.extract_fixed_cohort_profiles(
            record, model_id, self.adapter, data_dir=self.data_dir
        )
        result = alternatives.find_alternatives(record, model_id, fixed_profiles, policy)
        result["policy"] = copy.deepcopy(policy)
        record["method_policy"] = policy
        record["alternatives"] = result
        self._save(record)
        return copy.deepcopy(result)

    def compare(self, query_id, model_ids):
        record = self._load(query_id)
        result = compare.compare_models(record, model_ids)
        attachments = record.get("context") if isinstance(record.get("context"), dict) else {}
        attachments["last_comparison_ids"] = list(model_ids)
        record["context"] = attachments
        self._save(record)
        return result

    def export_comparison(self, query_id, model_ids, format):
        record = self._load(query_id)
        # Validates the same selection as the export and records it for the read-only GET route.
        result = research_export.export_research(record, model_ids, format)
        attachments = record.get("context") if isinstance(record.get("context"), dict) else {}
        attachments["last_comparison_ids"] = list(model_ids)
        record["context"] = attachments
        self._save(record)
        return result

    def export_latest(self, query_id, format):
        record = self._load(query_id)
        attachments = record.get("context") if isinstance(record.get("context"), dict) else {}
        model_ids = attachments.get("last_comparison_ids")
        if not isinstance(model_ids, list):
            raise ResearchStateError(
                "Select two or three models with /research/{query_id}/compare before exporting."
            )
        return research_export.export_research(record, model_ids, format)

    @staticmethod
    def _native_line(record, model_id):
        """Find one original native packet without adapting, scoring, or reading a table."""
        native = record.get("native_result", {})
        for bucket in _PARTITION_BUCKETS:
            for line in native.get(bucket, []) or []:
                if line.get("model_id") == model_id or line.get("ModelID") == model_id:
                    return bucket, line
        raise ValueError(f"Unknown model_id for this research query: {model_id!r}")

    def native_molecular_fallback(self, record, model_id):
        """Expose already-saved native molecular evidence as context, never a second score."""
        bucket, line = self._native_line(record, model_id)
        has_routes = bool(route_serialization.attachments(record.get("routes"))["by_relationship_id"])
        fallback = {
            "statement": (
                "native multi-omic suitability evidence; relationship-specific route evidence is shown separately and does not change native D or rank"
                if has_routes else
                "native multi-omic suitability evidence; no relationship-specific route was declared"
            ),
            "model_id": model_id,
            "partition_bucket": bucket,
            "native": {
                "D": copy.deepcopy(line.get("D")),
                "rank": copy.deepcopy(line.get("rank")),
                "confidence_tier": copy.deepcopy(line.get("confidence_tier")),
                "veto": copy.deepcopy(line.get("veto")),
                "per_gene": copy.deepcopy(line.get("per_gene", [])),
            },
        }
        context_attachment = record.get("context")
        if isinstance(context_attachment, dict) and "graph_context" in context_attachment:
            fallback["graph_pathway_context"] = copy.deepcopy(context_attachment["graph_context"])
        return fallback

    def _evaluate_route_record(self, record, relationship_type, source_gene, target_gene, context_type=None):
        """Evaluate one declared relationship against an already-loaded snapshot record."""
        original_native = copy.deepcopy(record.get("native_result"))
        kwargs = {"data_dir": self.data_dir}
        if self.route_adapter_function is not None:
            kwargs["adapter_function"] = self.route_adapter_function
        family, evidence = route_router.adapt_declared_route(
            record, relationship_type, source_gene, target_gene, context_type, **kwargs
        )
        relationship = evidence.get("relationship")
        policy_version = evidence.get("route_policy_version")
        if not isinstance(relationship, dict) or not isinstance(policy_version, str):
            raise ValueError("resolved BR02 route evidence requires relationship and policy version")
        stored_relationship_id = Relationship(
            family,
            ResolvedGene(**relationship["source_gene"]),
            ResolvedGene(**relationship["target_gene"]),
            relationship.get("context_type"),
        ).relationship_id(policy_version)
        existing = route_serialization.attachments(record.get("routes"))["by_relationship_id"].get(stored_relationship_id)
        if existing is not None:
            if original_native != record.get("native_result"):
                raise RuntimeError("route adaptation attempted to mutate immutable native_result")
            entry, idempotent = copy.deepcopy(existing), True
        else:
            result = route_router.evaluate_adapter_evidence(family, evidence, self.route_evaluators)
            attachments, entry, idempotent = route_serialization.put_attachment(record.get("routes"), result)
            record["routes"] = attachments
        if original_native != record.get("native_result"):
            raise RuntimeError("route evaluation attempted to mutate immutable native_result")
        return {
            "query_id": record["query_id"],
            "idempotent": idempotent,
            "request_orientation": {
                "type": relationship_type,
                "source_gene": source_gene,
                "target_gene": target_gene,
                "context_type": context_type,
            },
            "stored_request_orientation": copy.deepcopy(entry["request_orientation"]),
            "route": copy.deepcopy(entry["result"]),
        }

    def evaluate_route(self, query_id, relationship_type, source_gene, target_gene, context_type=None):
        """The manual research/debug evaluation trigger retained from BR07."""
        record = self._load(query_id)
        result = self._evaluate_route_record(record, relationship_type, source_gene, target_gene, context_type)
        if not result["idempotent"]:
            self._save(record)
        return result

    @staticmethod
    def _refresh_automatic_plan_state(plan):
        """Recount task state after one bounded transition; route outcome is not a rank input."""
        tasks = plan["tasks"]
        pending = sum(task["status"] == "pending" for task in tasks)
        completed = sum(task["status"] == "completed" for task in tasks)
        failed = sum(task["status"] == "failed" for task in tasks)
        plan["state"].update({
            "status": "completed" if pending == 0 else "in_progress",
            "pending_task_count": pending,
            "completed_task_count": completed,
            "failed_task_count": failed,
        })

    @staticmethod
    def _automatic_not_computed(record):
        """Represent absence truthfully without creating a plan from a read-only GET."""
        return {
            "query_id": record["query_id"],
            "automatic_routes_status": "not_computed",
            "reason": "Automatic route execution has not been started for this saved query.",
            "plan": None,
            "results": [],
            "conflicts": [],
            "limitations": [],
        }

    @staticmethod
    def _automatic_execution(record):
        """Load strict saved state only; callers choose whether a missing plan may be created."""
        value = record.get("automatic_routes")
        if value is None:
            return None
        if not isinstance(value, dict) or value.get("schema_version") != _AUTOMATIC_EXECUTION_SCHEMA_VERSION:
            raise ValueError("snapshot automatic routes attachment has an unsupported schema")
        plan = value.get("plan")
        validate_persisted_plan(plan)
        if value.get("policy_version") != plan.get("policy_version") or value.get("policy_identity") != plan.get("policy_identity"):
            raise ValueError("automatic execution policy identity does not match its saved plan")
        results = value.get("results_by_task_id")
        if not isinstance(results, dict):
            raise ValueError("automatic execution requires results_by_task_id")
        task_ids = {task["task_id"] for task in plan["tasks"]}
        if not set(results).issubset(task_ids):
            raise ValueError("automatic execution result references an unknown task ID")
        for task in plan["tasks"]:
            if task["status"] in {"completed", "failed"} and task["task_id"] not in results:
                raise ValueError("automatic execution completed task has no persisted result")
            if task["status"] == "pending" and task["task_id"] in results:
                raise ValueError("automatic execution pending task already has a result")
        return value

    @staticmethod
    def _automatic_response(record, execution):
        """Return a detached progress payload whose task ledger stays available for export."""
        plan = execution["plan"]
        return {
            "query_id": record["query_id"],
            "automatic_routes_status": plan["state"]["status"],
            "policy_version": plan["policy_version"],
            "policy_identity": plan["policy_identity"],
            "plan": copy.deepcopy(plan),
            "results": [copy.deepcopy(execution["results_by_task_id"][task["task_id"]])
                        for task in plan["tasks"] if task["task_id"] in execution["results_by_task_id"]],
            "conflicts": copy.deepcopy(execution.get("conflicts", [])),
            "limitations": copy.deepcopy(execution.get("limitations", [])),
        }

    def _new_automatic_execution(self, record):
        """Persist the AUTO02 plan once, pinning the exact policy bytes that selected its tasks."""
        plan = create_automatic_plan(record)
        validate_persisted_plan(plan)
        return {
            "schema_version": _AUTOMATIC_EXECUTION_SCHEMA_VERSION,
            "policy_version": plan["policy_version"],
            "policy_identity": plan["policy_identity"],
            "plan": plan,
            "results_by_task_id": {},
            "conflicts": [],
            "limitations": [
                "Automatic route evidence is supporting or conflicting context only and does not change native D, rank, tier, veto, candidate membership, or order."
            ],
        }

    @staticmethod
    def _automatic_error_result(task, exc):
        """Bound an evaluator failure as saved unavailable evidence rather than failing native output."""
        message = str(exc).replace("\n", " ").strip()[:_AUTOMATIC_ERROR_LIMIT]
        return {
            "task_id": task["task_id"],
            "status": "error",
            "route_status": "unavailable",
            "route_qualified": False,
            "relationship": {
                "route_family": task["route_family"],
                "source_ensembl_id": task["source_ensembl_id"],
                "target_ensembl_id": task["target_ensembl_id"],
                "context_type": task["context_type"],
            },
            "error": {"type": type(exc).__name__, "message": message or "route evaluation failed"},
            "limitations": ["The automatic route evaluator failed; native result remains unchanged."],
        }

    def get_automatic_routes(self, query_id):
        """Read persisted automatic progress only; this endpoint never creates or evaluates work."""
        record = self._load(query_id)
        execution = self._automatic_execution(record)
        return self._automatic_not_computed(record) if execution is None else self._automatic_response(record, execution)

    def advance_automatic_routes(self, query_id):
        """Create/resume a deterministic plan and evaluate at most one still-pending task."""
        with self._automatic_lock:
            record = self._load(query_id)
            original_native = copy.deepcopy(record.get("native_result"))
            execution = self._automatic_execution(record)
            if execution is None:
                execution = self._new_automatic_execution(record)
                record["automatic_routes"] = execution
                self._refresh_automatic_plan_state(execution["plan"])
            else:
                current_plan = create_automatic_plan(record)
                if plan_identity(execution["plan"]) != plan_identity(current_plan):
                    raise ValueError("saved automatic plan does not match the current canonical query and policy")
                current_policy_identity = hashlib.sha256(automatic_policy_bytes()).hexdigest()
                if execution.get("policy_identity") != current_policy_identity:
                    raise ValueError("saved automatic plan policy identity does not match the installed policy")

            task = next((item for item in execution["plan"]["tasks"] if item["status"] == "pending"), None)
            if task is not None:
                try:
                    route = self._evaluate_route_record(
                        record, task["route_family"], task["source_ensembl_id"],
                        task["target_ensembl_id"], task["context_type"],
                    )
                    result = {
                        "task_id": task["task_id"], "status": "completed",
                        "idempotent_route_attachment": route["idempotent"],
                        "relationship_id": route["route"].get("relationship_id"),
                        "route": route["route"],
                    }
                    task["status"] = "completed"
                    if route["route"].get("route_status") in {"contradicted", "conflicting"}:
                        execution["conflicts"].append({
                            "task_id": task["task_id"], "relationship_id": result["relationship_id"],
                            "route_status": route["route"].get("route_status"),
                        })
                except Exception as exc:
                    result = self._automatic_error_result(task, exc)
                    task["status"] = "failed"
                execution["results_by_task_id"][task["task_id"]] = result
                self._refresh_automatic_plan_state(execution["plan"])

            record["automatic_routes"] = execution
            if original_native != record.get("native_result"):
                raise RuntimeError("automatic route execution attempted to mutate immutable native_result")
            self._save(record)
            return self._automatic_response(record, execution)

    def get_routes(self, query_id):
        """Return bounded route cards from storage; no adapter/engine/data-reader call occurs."""
        record = self._load(query_id)
        attachments = route_serialization.attachments(record.get("routes"))
        return {
            "query_id": record["query_id"],
            "routes_status": attachments["status"],
            "routes": route_serialization.summaries(attachments),
            "native_molecular_fallback": {
                "statement": "native multi-omic suitability evidence; no relationship-specific route was declared",
                "scope": "Request a model route packet to receive its saved native D/rank/tier/veto and per-gene evidence.",
            },
        }

    def get_model_routes(self, query_id, model_id):
        """Read matching stored packets plus the requested model's saved native context."""
        record = self._load(query_id)
        self._require_model(record, model_id)
        attachments = route_serialization.attachments(record.get("routes"))
        packets = []
        for relationship_id, entry in attachments["by_relationship_id"].items():
            route = entry.get("result", {})
            packet = next(
                (copy.deepcopy(item) for item in route.get("model_packets", []) if item.get("ModelID") == model_id),
                None,
            )
            packets.append({
                "relationship_id": relationship_id,
                "relationship": copy.deepcopy(entry.get("request_orientation")),
                "route_status": route.get("route_status"),
                "route_qualified": route.get("route_qualified"),
                "model_packet": packet,
            })
        return {
            "query_id": record["query_id"],
            "model_id": model_id,
            "routes_status": attachments["status"],
            "routes": packets,
            "native_molecular_fallback": self.native_molecular_fallback(record, model_id),
        }

    def export_routes(self, query_id, format):
        """Export attached evidence only; this is deliberately not an evaluation trigger."""
        record = self._load(query_id)
        return research_export.export_routes(
            record, format, route_registry=load_route_registry(), route_controls=load_controls()
        )
