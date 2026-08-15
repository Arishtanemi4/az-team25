"""Response models for /filters endpoints."""

from pydantic import BaseModel


class LineageOptions(BaseModel):
    """Sorted list of distinct lineage values, for the lineage filter dropdown."""

    lineages: list[str]


class PrimaryDiseaseOptions(BaseModel):
    """Sorted list of distinct primary_disease values, optionally scoped to a lineage."""

    primary_diseases: list[str]
