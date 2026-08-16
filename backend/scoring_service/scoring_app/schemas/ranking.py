"""Request model for POST /rank. The response is intentionally NOT re-declared field-by-field
here -- it is exactly scoring/export.py::export_query_result's output contract, and duplicating
that as nested Pydantic models would create a second copy of the contract that could silently
drift from the one scoring/ actually produces. See app/services/ranking_service.py.
"""

from pydantic import BaseModel, Field, model_validator


class RankRequest(BaseModel):
    inclusion_genes: list[str] = Field(
        default_factory=list, description="Gene symbols or Ensembl IDs that must be present/active"
    )
    exclusion_genes: list[str] = Field(
        default_factory=list, description="Gene symbols or Ensembl IDs that must be absent/inactive"
    )
    lineage: str | None = Field(default=None, description="Restrict candidates to one lineage")
    primary_disease: str | None = Field(
        default=None, description="Restrict candidates to one primary_disease"
    )
    exclude_problematic: bool = Field(
        default=False, description="Remove cell lines flagged as problematic before scoring"
    )
    msi_high: bool | None = Field(
        default=None, description="Restrict candidates by MSI-high status when specified"
    )
    ploidy_min: float | None = Field(
        default=None, description="Minimum supported ploidy value"
    )
    ploidy_max: float | None = Field(
        default=None, description="Maximum supported ploidy value"
    )
    require_metabolomics: bool = Field(
        default=False, description="Keep only lines with metabolomics availability"
    )
    require_mirna: bool = Field(
        default=False, description="Keep only lines with miRNA availability"
    )
    top_k: int = Field(default=10, ge=1, le=100, description="How many ranked lines to return")

    @model_validator(mode="after")
    def validate_research_query(self):
        # The native scorer can inspect exclusions, but confidence/eligibility are defined around
        # at least one desired molecular target. Do not present exclusion-only output as a ranking.
        if not self.inclusion_genes:
            raise ValueError(
                "At least one inclusion gene is required; exclusion genes are optional criteria."
            )
        if self.ploidy_min is not None and self.ploidy_max is not None:
            if self.ploidy_min > self.ploidy_max:
                raise ValueError("ploidy_min must be less than or equal to ploidy_max.")
        return self
