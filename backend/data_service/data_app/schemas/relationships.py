"""Response models for /data/relationships/*."""

from pydantic import BaseModel


class GeneCellLinePoint(BaseModel):
    model_id: str
    cell_line_name: str | None
    lineage: str | None
    value: float


class GeneCellLineResponse(BaseModel):
    ensembl_id: str
    symbol: str | None
    layer: str
    unit: str
    n_measured: int
    n_models_total: int
    points: list[GeneCellLinePoint]


class CellLineDiseaseGroup(BaseModel):
    lineage: str
    primary_disease: str
    n_models: int


class CellLineDiseaseResponse(BaseModel):
    n_models_total: int
    cells: list[CellLineDiseaseGroup]


class GeneDiseaseGroup(BaseModel):
    primary_disease: str
    n_measured: int
    n_models: int
    median: float | None
    iqr_low: float | None
    iqr_high: float | None
    flagged: bool


class GeneDiseaseResponse(BaseModel):
    ensembl_id: str
    symbol: str | None
    layer: str
    min_n: int
    groups: list[GeneDiseaseGroup]
