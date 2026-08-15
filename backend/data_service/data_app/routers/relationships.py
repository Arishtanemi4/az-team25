"""GET /data/relationships/* -- gene<->cell-line, cell-line<->disease, gene<->disease
(PRODUCT_SURFACE.md SS3.2). Every response states its denominator; a group below `min_n` is
returned flagged, never dropped."""

from fastapi import APIRouter, Query

from data_app.schemas.relationships import (
    CellLineDiseaseResponse, GeneCellLineResponse, GeneDiseaseResponse,
)
from data_app.services import relationships_service

router = APIRouter(prefix="/data/relationships", tags=["relationships"])


@router.get("/gene-cell-line", response_model=GeneCellLineResponse)
def gene_cell_line(
    ensembl_id: str = Query(...),
    layer: str = Query("rna"),
    limit: int = Query(2000, ge=1, le=2000),
):
    return relationships_service.gene_cell_line(ensembl_id, layer, limit)


@router.get("/cell-line-disease", response_model=CellLineDiseaseResponse)
def cell_line_disease(lineage: str | None = Query(None)):
    return relationships_service.cell_line_disease(lineage)


@router.get("/gene-disease", response_model=GeneDiseaseResponse)
def gene_disease(
    ensembl_id: str = Query(...),
    layer: str = Query("rna"),
    min_n: int = Query(5, ge=1),
):
    return relationships_service.gene_disease(ensembl_id, layer, min_n)
