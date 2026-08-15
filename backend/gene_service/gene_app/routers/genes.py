"""GET /genes/search -- the endpoint behind the search+dropdown gene input on the frontend."""

from fastapi import APIRouter, Query, Request

from gene_app.schemas.genes import GeneMatch

router = APIRouter(prefix="/genes", tags=["genes"])


@router.get("/search", response_model=list[GeneMatch])
def search_genes(
    request: Request,
    q: str = Query(
        ..., min_length=0, description="Partial gene symbol or Ensembl ID; empty returns a default alphabetical subset"
    ),
    limit: int = Query(20, ge=1, le=100),
):
    service = request.app.state.gene_reference_service
    return service.search_genes(q, limit=limit)
