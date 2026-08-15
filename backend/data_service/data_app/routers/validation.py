"""GET /data/validation/* -- the 4-study external validation manifest and each study's own
result CSV, reproduced as-is (PRODUCT_SURFACE.md SS3.4)."""

from fastapi import APIRouter

from data_app.schemas.validation import ValidationStudiesResponse, ValidationStudyResponse
from data_app.services import validation_service

router = APIRouter(prefix="/data/validation", tags=["validation"])


@router.get("/studies", response_model=ValidationStudiesResponse)
def studies():
    return validation_service.list_studies()


@router.get("/{study}", response_model=ValidationStudyResponse)
def study(study: str):
    return validation_service.get_study(study)
