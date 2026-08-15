"""GET /filters/lineages and /filters/primary-diseases -- populate the biological/disease-context
filter dropdowns without ever asking the user to pick an algorithm (root _.md SS1 item 5)."""

from fastapi import APIRouter, Query, Request

from gene_app.schemas.filters import LineageOptions, PrimaryDiseaseOptions

router = APIRouter(prefix="/filters", tags=["filters"])


@router.get("/lineages", response_model=LineageOptions)
def get_lineages(request: Request):
    service = request.app.state.gene_reference_service
    return LineageOptions(lineages=service.list_lineages())


@router.get("/primary-diseases", response_model=PrimaryDiseaseOptions)
def get_primary_diseases(request: Request, lineage: str | None = Query(default=None)):
    service = request.app.state.gene_reference_service
    return PrimaryDiseaseOptions(primary_diseases=service.list_primary_diseases(lineage=lineage))
