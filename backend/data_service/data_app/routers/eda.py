"""GET /data/eda/* -- per-layer coverage, lineage x layer coverage bias, and precomputed
concordance aggregates (PRODUCT_SURFACE.md SS3.3)."""

from fastapi import APIRouter

from data_app.schemas.eda import CoverageResponse, LineageCoverageResponse
from data_app.services import eda_service

router = APIRouter(prefix="/data/eda", tags=["eda"])


@router.get("/coverage", response_model=CoverageResponse)
def coverage():
    return eda_service.coverage()


@router.get("/lineage-coverage", response_model=LineageCoverageResponse)
def lineage_coverage():
    return eda_service.lineage_coverage()


@router.get("/concordance")
def concordance():
    # available:false is a first-class response shape (PRODUCT_SURFACE.md SS3.3), so this
    # endpoint deliberately has no single response_model -- Pydantic would force one shape.
    return eda_service.concordance()
