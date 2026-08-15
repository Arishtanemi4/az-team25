"""Response models for /data/validation/*."""

from pydantic import BaseModel


class ValidationStudiesResponse(BaseModel):
    studies: list[dict]


class ValidationStudyResponse(BaseModel):
    study: str
    total_rows: int
    returned_rows: int
    columns: list[str]
    rows: list[dict]
