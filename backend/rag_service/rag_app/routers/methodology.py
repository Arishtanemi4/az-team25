"""POST /methodology/ask -- answers a researcher's question about the scoring methodology from
project documentation (rag/methodology_agent.py)."""

from typing import Any

from fastapi import APIRouter, HTTPException, Request

from rag_app.schemas.rag import MethodologyQuestionRequest

router = APIRouter(tags=["methodology"])


@router.post("/methodology/ask")
def ask(request: Request, body: MethodologyQuestionRequest) -> dict[str, Any]:
    service = request.app.state.rag_service
    try:
        return service.answer_methodology_question(body.question)
    except RuntimeError as exc:
        # e.g. NVIDIA_API_KEY is not set -- a service misconfiguration, not a bad request.
        raise HTTPException(status_code=500, detail=str(exc)) from exc
