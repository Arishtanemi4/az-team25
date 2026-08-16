"""Explicit BR07 dispatch for immutable snapshot route evidence.

There is intentionally no "best route" selection here.  A POST supplies one declared
relationship, BR02 projects evidence for that declaration, and exactly its registered evaluator
is called once.  GET handlers consume the persisted result through :mod:`service` only.
"""

from __future__ import annotations

from typing import Any, Callable, Mapping

from .codependency import evaluate_codependency
from .conditional_dependency import evaluate_conditional_dependency
from .contracts import ContextType, EvidenceState, RouteFamily
from .direct_fusion import evaluate_direct_fusion
from .evidence_adapter import adapt_route_evidence


def validate_declaration(
    relationship_type: str, source_gene: str, target_gene: str, context_type: str | None
) -> tuple[RouteFamily, ContextType | None]:
    """Reject malformed declarations before BR02 can touch a data reader."""
    try:
        family = RouteFamily(relationship_type)
    except (TypeError, ValueError) as exc:
        raise ValueError("type must be direct_fusion, conditional_dependency, or crispr_codependency") from exc
    if not isinstance(source_gene, str) or not source_gene.strip():
        raise ValueError("source_gene must be a non-empty string")
    if not isinstance(target_gene, str) or not target_gene.strip():
        raise ValueError("target_gene must be a non-empty string")
    if source_gene.strip() == target_gene.strip():
        raise ValueError("source_gene and target_gene must be distinct")
    if family is RouteFamily.CONDITIONAL_DEPENDENCY:
        try:
            context = ContextType(context_type)
        except (TypeError, ValueError) as exc:
            raise ValueError("conditional_dependency requires context_type mutation or rna_high") from exc
        return family, context
    if context_type is not None:
        raise ValueError(f"{family.value} does not accept a context_type")
    return family, None


def adapt_declared_route(
    snapshot: Mapping[str, Any],
    relationship_type: str,
    source_gene: str,
    target_gene: str,
    context_type: str | None = None,
    *,
    data_dir: str | None = None,
    adapter_function: Callable[..., dict[str, Any]] = adapt_route_evidence,
    evaluators: Mapping[RouteFamily, Callable[[Mapping[str, Any]], dict[str, Any]]] | None = None,
) -> tuple[RouteFamily, dict[str, Any]]:
    """Resolve/adapt one declaration so callers can check an existing attachment first."""
    family, context = validate_declaration(relationship_type, source_gene, target_gene, context_type)
    evidence = adapter_function(
        snapshot, family, source_gene.strip(), target_gene.strip(),
        context.value if context else None, data_dir=data_dir,
    )
    if evidence.get("evidence_state") == EvidenceState.UNRESOLVABLE.value:
        raise ValueError("source_gene or target_gene is unresolved or ambiguous in frozen gene_reference")
    return family, evidence


def evaluate_adapter_evidence(
    family: RouteFamily, evidence: Mapping[str, Any],
    evaluators: Mapping[RouteFamily, Callable[[Mapping[str, Any]], dict[str, Any]]] | None = None,
) -> dict[str, Any]:
    """Dispatch exactly one fixed evaluator over an already-created BR02 projection."""
    dispatch = evaluators or {
        RouteFamily.DIRECT_FUSION: evaluate_direct_fusion,
        RouteFamily.CONDITIONAL_DEPENDENCY: evaluate_conditional_dependency,
        RouteFamily.CRISPR_CODEPENDENCY: evaluate_codependency,
    }
    evaluator = dispatch.get(family)
    if evaluator is None:
        raise ValueError(f"No evaluator registered for {family.value}")
    result = evaluator(evidence)
    if result.get("route_family") != family.value:
        raise RuntimeError("route evaluator returned a mismatched route family")
    return result


def evaluate_declared_route(
    snapshot: Mapping[str, Any], relationship_type: str, source_gene: str, target_gene: str,
    context_type: str | None = None, *, data_dir: str | None = None,
    adapter_function: Callable[..., dict[str, Any]] = adapt_route_evidence,
    evaluators: Mapping[RouteFamily, Callable[[Mapping[str, Any]], dict[str, Any]]] | None = None,
) -> dict[str, Any]:
    """Convenience one-shot API retained for direct router tests and non-persistent callers."""
    family, evidence = adapt_declared_route(
        snapshot, relationship_type, source_gene, target_gene, context_type,
        data_dir=data_dir, adapter_function=adapter_function,
    )
    return evaluate_adapter_evidence(family, evidence, evaluators)
