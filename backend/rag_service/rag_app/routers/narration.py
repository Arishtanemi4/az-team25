"""POST /narrate -- turns one query's ranked result set into a grounded, verifier-checked
narrative (rag/narrator.py)."""

from typing import Any

from fastapi import APIRouter, HTTPException, Request

from rag_app.schemas.rag import NarrateRequest

router = APIRouter(tags=["narration"])


@router.post("/narrate")
def narrate(request: Request, body: NarrateRequest) -> dict[str, Any]:
    service = request.app.state.rag_service
    try:
        return service.narrate_result(body.evidence_record, top_k=body.top_k_context)
    except RuntimeError as exc:
        # narrate()'s own documented failure mode: two attempts both failed schema validation or
        # the grounding check -- a genuine backend failure, not a bad request.
        raise HTTPException(status_code=502, detail=str(exc)) from exc
