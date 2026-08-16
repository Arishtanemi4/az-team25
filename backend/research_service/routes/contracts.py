"""BR01 immutable contracts for biological evidence routes.

These objects validate declared, already-resolved relationships and individual long-ledger
evidence rows.  They do not read data, evaluate a route, calculate a score, or import validation
code.  Native D, rank, tier and veto values remain outside this package.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
import hashlib
import json
import math
from typing import Any, Mapping


class _StringEnum(str, Enum):
    """String-valued enums make policy JSON and strict serialization unambiguous."""

    def __str__(self) -> str:
        return self.value


class RouteFamily(_StringEnum):
    DIRECT_FUSION = "direct_fusion"
    CONDITIONAL_DEPENDENCY = "conditional_dependency"
    CRISPR_CODEPENDENCY = "crispr_codependency"


class ContextType(_StringEnum):
    MUTATION = "mutation"
    RNA_HIGH = "rna_high"


class EvidenceState(_StringEnum):
    OBSERVED = "observed"
    MEASURED_ABSENT = "measured_absent"
    NOT_ASSAYED = "not_assayed"
    UNRESOLVABLE = "unresolvable"
    NOT_APPLICABLE = "not_applicable"


class RouteStatus(_StringEnum):
    SUPPORTED = "supported"
    CONTRADICTED = "contradicted"
    EXPLORATORY = "exploratory"
    UNAVAILABLE = "unavailable"


class ModelSupportState(_StringEnum):
    SUPPORT = "support"
    NEUTRAL = "neutral"
    CONFLICT = "conflict"
    UNKNOWN = "unknown"
    NOT_APPLICABLE = "not_applicable"


class SourceRole(_StringEnum):
    RANKING_INPUT = "ranking_input"
    ROUTE_INPUT = "route_input"
    CONTEXT_ONLY = "context_only"
    VALIDATION_ONLY = "validation_only"


class RouteInvocation(_StringEnum):
    """The declared route entry point, kept separate from route-family evidence semantics."""

    AUTOMATIC_BACKGROUND = "automatic_background"
    MANUAL_RESEARCH = "manual_research"


def _canonical_json(value: Any) -> bytes:
    """Encode only finite JSON values so deterministic IDs cannot conceal NaN/Infinity."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _require_nonempty(value: str, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty string")
    return value.strip()


@dataclass(frozen=True)
class ResolvedGene:
    """A canonical gene-reference result; callers may not provide a symbol-only guess."""

    symbol: str
    ensembl_id: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "symbol", _require_nonempty(self.symbol, "symbol"))
        object.__setattr__(self, "ensembl_id", _require_nonempty(self.ensembl_id, "ensembl_id"))

    def as_dict(self) -> dict[str, str]:
        return {"symbol": self.symbol, "ensembl_id": self.ensembl_id}


@dataclass(frozen=True)
class Relationship:
    """A user-declared relationship after both genes resolved through gene_reference.csv."""

    route_family: RouteFamily
    source_gene: ResolvedGene
    target_gene: ResolvedGene
    context_type: ContextType | None = None

    def __post_init__(self) -> None:
        family = RouteFamily(self.route_family)
        context = ContextType(self.context_type) if self.context_type is not None else None
        object.__setattr__(self, "route_family", family)
        object.__setattr__(self, "context_type", context)

        if self.source_gene.ensembl_id == self.target_gene.ensembl_id:
            raise ValueError("source_gene and target_gene must resolve to distinct Ensembl IDs")
        if family is RouteFamily.CONDITIONAL_DEPENDENCY and context is None:
            raise ValueError("conditional_dependency requires context_type mutation or rna_high")
        if family is not RouteFamily.CONDITIONAL_DEPENDENCY and context is not None:
            raise ValueError(f"{family.value} does not accept a context_type")

    @property
    def is_symmetric(self) -> bool:
        return self.route_family in {
            RouteFamily.DIRECT_FUSION,
            RouteFamily.CRISPR_CODEPENDENCY,
        }

    @property
    def canonical_gene_pair(self) -> tuple[str, str]:
        """The unordered identity used by fusion/codependency while request orientation is retained."""
        return tuple(sorted((self.source_gene.ensembl_id, self.target_gene.ensembl_id)))

    def as_dict(self) -> dict[str, Any]:
        return {
            "route_family": self.route_family.value,
            "source_gene": self.source_gene.as_dict(),
            "target_gene": self.target_gene.as_dict(),
            "context_type": self.context_type.value if self.context_type else None,
            "canonical_gene_pair": list(self.canonical_gene_pair),
        }

    def relationship_id(self, policy_version: str) -> str:
        """Stable attachment key: symmetric pairs canonicalize, conditional direction does not."""
        policy_version = _require_nonempty(policy_version, "policy_version")
        if self.is_symmetric:
            identity = {
                "route_family": self.route_family.value,
                "canonical_gene_pair": self.canonical_gene_pair,
                "context_type": None,
            }
        else:
            identity = {
                "route_family": self.route_family.value,
                "source_ensembl_id": self.source_gene.ensembl_id,
                "target_ensembl_id": self.target_gene.ensembl_id,
                "context_type": self.context_type.value if self.context_type else None,
            }
        payload = {"policy_version": policy_version, "relationship": identity}
        return hashlib.sha256(_canonical_json(payload)).hexdigest()


def _json_value_is_finite(value: Any) -> bool:
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, Mapping):
        return all(isinstance(key, str) and _json_value_is_finite(item) for key, item in value.items())
    if isinstance(value, (list, tuple)):
        return all(_json_value_is_finite(item) for item in value)
    return value is None or isinstance(value, (str, int, bool))


@dataclass(frozen=True)
class EvidenceLedgerRecord:
    """One evidence item in the route CSV/JSON long ledger.

    The tuple `(research_query_id, route_family, model_id, ensembl_id, evidence_id)` is its
    immutable key.  `evidence_id` is derived from the complete row contents, so event multiplicity
    is retained rather than silently reduced to a model/gene summary.
    """

    research_query_id: str
    route_family: RouteFamily
    model_id: str
    ensembl_id: str
    relationship: Relationship
    route_policy_version: str
    source_id: str
    source_version: str
    data_id: str
    raw_value: Any
    unit: str
    transformation: str
    denominator: Mapping[str, Any]
    evidence_state: EvidenceState
    qc_state: str
    limitation: str
    source_record_id: str | None = None
    evidence_id: str = field(init=False)

    def __post_init__(self) -> None:
        for field_name in (
            "research_query_id",
            "model_id",
            "ensembl_id",
            "route_policy_version",
            "source_id",
            "source_version",
            "data_id",
            "unit",
            "transformation",
            "qc_state",
            "limitation",
        ):
            object.__setattr__(self, field_name, _require_nonempty(getattr(self, field_name), field_name))
        if self.source_record_id is not None:
            object.__setattr__(
                self, "source_record_id", _require_nonempty(self.source_record_id, "source_record_id")
            )
        object.__setattr__(self, "route_family", RouteFamily(self.route_family))
        object.__setattr__(self, "evidence_state", EvidenceState(self.evidence_state))
        if self.relationship.route_family is not self.route_family:
            raise ValueError("ledger route_family must match relationship route_family")
        if not isinstance(self.denominator, Mapping) or not self.denominator:
            raise ValueError("denominator must be a non-empty mapping with explicit cohort counts")
        if not _json_value_is_finite(self.raw_value) or not _json_value_is_finite(self.denominator):
            raise ValueError("ledger raw_value and denominator must contain finite JSON values")
        object.__setattr__(self, "evidence_id", self._make_evidence_id())

    def _id_payload(self) -> dict[str, Any]:
        return {
            "research_query_id": self.research_query_id,
            "route_family": self.route_family.value,
            # Preserve the processed/native spine's public column spelling in every export row.
            "ModelID": self.model_id,
            "ensembl_id": self.ensembl_id,
            "relationship": self.relationship.as_dict(),
            "route_policy_version": self.route_policy_version,
            "source_id": self.source_id,
            "source_version": self.source_version,
            "data_id": self.data_id,
            "raw_value": self.raw_value,
            "unit": self.unit,
            "transformation": self.transformation,
            "denominator": dict(self.denominator),
            "evidence_state": self.evidence_state.value,
            "qc_state": self.qc_state,
            "limitation": self.limitation,
            "source_record_id": self.source_record_id,
        }

    def _make_evidence_id(self) -> str:
        return hashlib.sha256(_canonical_json(self._id_payload())).hexdigest()

    @property
    def key(self) -> tuple[str, str, str, str, str]:
        return (
            self.research_query_id,
            self.route_family.value,
            self.model_id,
            self.ensembl_id,
            self.evidence_id,
        )

    def as_dict(self) -> dict[str, Any]:
        row = self._id_payload()
        row["relationship_id"] = self.relationship.relationship_id(self.route_policy_version)
        row["evidence_id"] = self.evidence_id
        return row


def validate_ledger(records: list[EvidenceLedgerRecord]) -> None:
    """Reject duplicate long-ledger keys before later export can overwrite evidence events."""
    keys = [record.key for record in records]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate evidence ledger key")
