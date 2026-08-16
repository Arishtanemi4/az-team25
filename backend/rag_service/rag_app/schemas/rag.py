"""Request models for rag_service's four endpoints. Responses are intentionally NOT re-declared
as pydantic models here -- each is exactly what the corresponding rag/ entry point already
returns (a plain JSON-serializable dict), and duplicating that shape as nested pydantic models
would create a second copy of the contract that could silently drift. See
app/services/rag_service.py.
"""

from typing import Any

from pydantic import BaseModel, Field


class NarrateRequest(BaseModel):
    # Exactly the /rank response from scoring_service, passed through verbatim -- narrator.py's
    # evidence_record parameter.
    evidence_record: dict[str, Any]
    top_k_context: int = Field(default=5, ge=1, le=20)


class MethodologyQuestionRequest(BaseModel):
    question: str = Field(min_length=1)


class QueryExpansionRequest(BaseModel):
    inclusion_genes: list[str] = Field(min_length=1)


class LiteratureQuestionRequest(BaseModel):
    question: str = Field(min_length=1)
