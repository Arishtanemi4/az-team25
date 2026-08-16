"""Load and validate BR01's versioned, fixed route policy registry.

Only sources explicitly registered as ``route_input`` can be returned to a route evaluator.  In
particular, validation-only sources are rejected before any later route stage can consume them.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping

from .contracts import RouteFamily, RouteInvocation, SourceRole


POLICY_DIRECTORY = Path(__file__).resolve().parents[1] / "route_policies"
ROUTE_REGISTRY_PATH = POLICY_DIRECTORY / "route_registry.json"
CONTROL_QUERIES_PATH = POLICY_DIRECTORY / "control_queries.json"
VALIDATION_SOURCES_PATH = POLICY_DIRECTORY / "validation_sources.json"
AUTOMATIC_ROUTE_POLICY_PATH = POLICY_DIRECTORY / "automatic_route_policy.json"


_AUTOMATIC_TASK_BATTERY = (
    ("direct_fusion", "direct_fusion", None, "unordered"),
    ("crispr_codependency", "crispr_codependency", None, "unordered"),
    ("conditional_mutation_first_to_second", "conditional_dependency", "mutation", "first_to_second"),
    ("conditional_mutation_second_to_first", "conditional_dependency", "mutation", "second_to_first"),
    ("conditional_rna_high_first_to_second", "conditional_dependency", "rna_high", "first_to_second"),
    ("conditional_rna_high_second_to_first", "conditional_dependency", "rna_high", "second_to_first"),
)


@dataclass(frozen=True)
class SourceRegistration:
    source_id: str
    source_version: str
    role: SourceRole
    allowed_route_families: tuple[RouteFamily, ...]
    purpose: str


def _load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _nonempty_string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string")
    return value.strip()


def _parse_source(item: Mapping[str, Any]) -> SourceRegistration:
    families = tuple(RouteFamily(family) for family in item.get("allowed_route_families", []))
    role = SourceRole(item["role"])
    if role is SourceRole.ROUTE_INPUT and not families:
        raise ValueError(f"route_input {item.get('source_id')} requires allowed_route_families")
    if role is not SourceRole.ROUTE_INPUT and families:
        raise ValueError(f"{role.value} source {item.get('source_id')} cannot declare route families")
    return SourceRegistration(
        source_id=_nonempty_string(item["source_id"], "source_id"),
        source_version=_nonempty_string(item["source_version"], "source_version"),
        role=role,
        allowed_route_families=families,
        purpose=_nonempty_string(item["purpose"], "purpose"),
    )


def _validate_method_policy(registry: Mapping[str, Any]) -> None:
    policy = registry.get("method_policy")
    if not isinstance(policy, Mapping):
        raise ValueError("route registry requires method_policy")
    expected = {
        "random_seed": 25,
        "bootstrap_replicates": 500,
        "bootstrap_percentiles": [0.025, 0.975],
        "lineage_minimum_n": 15,
        "context_minimum_n": 10,
        "reference_minimum_n": 10,
        "total_or_shared_minimum_n": 30,
        "eligible_lineages_minimum": 2,
        "leave_one_lineage_out_sign_consistency_minimum": 0.80,
        "model_dependency_support_boundary": -0.5,
    }
    for name, value in expected.items():
        if policy.get(name) != value:
            raise ValueError(f"method_policy {name} must equal locked value {value!r}")


def load_route_registry(path: Path | None = None) -> dict[str, Any]:
    """Return the complete registry only if every immutable BR01 policy is internally valid."""
    registry = _load_json(path or ROUTE_REGISTRY_PATH)
    if registry.get("schema_version") != "biological-routes-policy-v1":
        raise ValueError("unsupported route registry schema_version")
    _nonempty_string(registry.get("policy_version"), "policy_version")
    _validate_method_policy(registry)
    sources = [_parse_source(item) for item in registry.get("sources", [])]
    source_ids = [source.source_id for source in sources]
    if not sources or len(source_ids) != len(set(source_ids)):
        raise ValueError("route registry sources must be non-empty with unique source_id values")
    route_sources = registry.get("route_sources")
    if not isinstance(route_sources, Mapping):
        raise ValueError("route registry requires route_sources")
    source_lookup = {source.source_id: source for source in sources}
    for family_name, source_ids_for_family in route_sources.items():
        family = RouteFamily(family_name)
        if not isinstance(source_ids_for_family, list) or not source_ids_for_family:
            raise ValueError(f"route_sources {family.value} must be a non-empty list")
        for source_id in source_ids_for_family:
            source = source_lookup.get(source_id)
            if source is None:
                raise ValueError(f"route source {source_id} is not registered")
            assert_route_input(source, family)
    return registry


def load_controls(path: Path | None = None) -> dict[str, Any]:
    """Load pre-registered, deterministic BR09 controls without evaluating them."""
    controls = _load_json(path or CONTROL_QUERIES_PATH)
    if controls.get("schema_version") != "biological-route-controls-v1":
        raise ValueError("unsupported control query schema_version")
    control_rows = controls.get("controls")
    if not isinstance(control_rows, list) or not control_rows:
        raise ValueError("control_queries requires a non-empty controls list")
    ids = []
    for control in control_rows:
        control_id = _nonempty_string(control.get("control_id"), "control_id")
        ids.append(control_id)
        family = RouteFamily(control.get("relationship", {}).get("route_family"))
        context = control["relationship"].get("context_type")
        if family is RouteFamily.CONDITIONAL_DEPENDENCY and context not in {"mutation", "rna_high"}:
            raise ValueError(f"conditional control {control_id} requires mutation or rna_high context")
        if family is not RouteFamily.CONDITIONAL_DEPENDENCY and context is not None:
            raise ValueError(f"symmetric control {control_id} cannot carry context_type")
    if len(ids) != len(set(ids)):
        raise ValueError("control_id values must be unique")
    return controls


def load_validation_sources(path: Path | None = None) -> dict[str, Any]:
    """Load governed corroboration sources and reject any role other than validation_only."""
    policy = _load_json(path or VALIDATION_SOURCES_PATH)
    if policy.get("schema_version") != "biological-route-validation-sources-v1":
        raise ValueError("unsupported validation source schema_version")
    sources = policy.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError("validation_sources requires a non-empty sources list")
    for source in sources:
        if SourceRole(source.get("role")) is not SourceRole.VALIDATION_ONLY:
            raise ValueError("validation_sources may only declare validation_only roles")
    return policy


def _validate_automatic_task_battery(tasks: Any) -> None:
    """Require the full predeclared six-task pair battery, in its public execution order."""
    if not isinstance(tasks, list) or len(tasks) != len(_AUTOMATIC_TASK_BATTERY):
        raise ValueError("automatic route policy requires exactly six fixed tasks per unordered pair")
    actual = []
    for task in tasks:
        if not isinstance(task, Mapping):
            raise ValueError("automatic route policy task must be an object")
        task_id = _nonempty_string(task.get("task_id"), "automatic task_id")
        family = RouteFamily(task.get("route_family"))
        context = task.get("context_type")
        if context is not None and context not in {"mutation", "rna_high"}:
            raise ValueError(f"automatic task {task_id} has unsupported context_type")
        direction = _nonempty_string(task.get("pair_direction"), "automatic pair_direction")
        actual.append((task_id, family.value, context, direction))
    if tuple(actual) != _AUTOMATIC_TASK_BATTERY:
        raise ValueError("automatic route policy task battery must equal the fixed six-task battery")


def load_automatic_route_policy(path: Path | None = None) -> dict[str, Any]:
    """Load AUTO01's fail-closed, result-independent background route-selection policy.

    ``manual_research`` remains the existing explicit-declaration path.  This loader deliberately
    does not plan or evaluate a route: AUTO02 will consume this validated declaration only.
    """
    policy = _load_json(path or AUTOMATIC_ROUTE_POLICY_PATH)
    if policy.get("schema_version") != "automatic-biological-route-policy-v1":
        raise ValueError("unsupported automatic route policy schema_version")
    _nonempty_string(policy.get("policy_version"), "automatic policy_version")
    if policy.get("invocation") != RouteInvocation.AUTOMATIC_BACKGROUND.value:
        raise ValueError("automatic route policy invocation must be automatic_background")
    if policy.get("manual_research_invocation") != RouteInvocation.MANUAL_RESEARCH.value:
        raise ValueError("automatic route policy must preserve manual_research invocation")
    if policy.get("selection_timing") != "before_route_outcomes":
        raise ValueError("automatic route policy must select tasks before route outcomes")
    if policy.get("pair_ordering") != (
        "Sort canonical Ensembl identifiers ascending and generate unordered pairs lexicographically."
    ):
        raise ValueError("automatic route policy requires deterministic canonical pair ordering")
    if policy.get("maximum_unordered_pairs") != 10:
        raise ValueError("automatic route policy maximum_unordered_pairs must equal 10")
    omission = policy.get("omission")
    if not isinstance(omission, Mapping) or omission.get("reason") != "automatic_pair_limit":
        raise ValueError("automatic route policy requires automatic_pair_limit omission rule")
    if omission.get("record_in_research_export") is not True:
        raise ValueError("automatic route policy must retain omitted pairs in the research export")
    if policy.get("single_gene_pair_route_status") != "not_applicable":
        raise ValueError("automatic route policy must mark single-gene pair routes not_applicable")
    if policy.get("invent_target_gene") is not False:
        raise ValueError("automatic route policy must prohibit invented target genes")
    if policy.get("exclusion_gene_pair_hypotheses") is not False:
        raise ValueError("automatic route policy must prohibit exclusion-gene pair hypotheses")
    _validate_automatic_task_battery(policy.get("fixed_task_battery"))
    boundary = policy.get("native_ranking_boundary")
    expected_boundary = ["D", "rank", "tier", "veto", "candidate_membership", "candidate_order"]
    if not isinstance(boundary, Mapping) or boundary.get("unchanged_fields") != expected_boundary:
        raise ValueError("automatic route policy must lock the native ranking boundary")
    if boundary.get("biological_routes_may_change_native_result") is not False:
        raise ValueError("automatic route policy must prohibit biological changes to native results")
    return policy


def assert_route_input(source: SourceRegistration, route_family: RouteFamily) -> None:
    """The hard gate preventing validation/context/ranking sources from becoming route inputs."""
    if source.role is not SourceRole.ROUTE_INPUT:
        raise ValueError(
            f"source {source.source_id} has role {source.role.value}; only route_input is consumable"
        )
    if route_family not in source.allowed_route_families:
        raise ValueError(f"source {source.source_id} is not registered for {route_family.value}")


def route_input_sources(route_family: RouteFamily, registry: Mapping[str, Any] | None = None) -> tuple[SourceRegistration, ...]:
    """Return only policy-permitted route inputs in registry order; no validation source can enter."""
    registry = dict(registry) if registry is not None else load_route_registry()
    parsed_sources = {_parse_source(item).source_id: _parse_source(item) for item in registry["sources"]}
    requested_ids = registry["route_sources"][RouteFamily(route_family).value]
    result = tuple(parsed_sources[source_id] for source_id in requested_ids)
    for source in result:
        assert_route_input(source, RouteFamily(route_family))
    return result


def policy_bytes(path: Path | None = None) -> bytes:
    """Canonical bytes for deterministic policy loading and future attachment provenance."""
    return json.dumps(load_route_registry(path), sort_keys=True, separators=(",", ":"), allow_nan=False).encode(
        "utf-8"
    )


def automatic_policy_bytes(path: Path | None = None) -> bytes:
    """Canonical AUTO01 bytes for persisted plan/export provenance in later stages."""
    return json.dumps(
        load_automatic_route_policy(path), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
