"""Response models for /data/eda/*."""

from pydantic import BaseModel


class CoverageResponse(BaseModel):
    coverage_state_counts: dict[str, dict[str, int]]


class LineageCoverageRow(BaseModel):
    lineage: str
    layer: str
    n: int
    n_measured: int
    measured_fraction: float
    flagged: bool


class LineageCoverageResponse(BaseModel):
    min_lineage_n: int
    rows: list[LineageCoverageRow]


class ConcordanceAvailable(BaseModel):
    available: bool = True
    generated_utc: str
    rna_cross_source: dict
    rna_protein: dict


class ConcordanceUnavailable(BaseModel):
    available: bool = False
    reason: str
    build_command: str
