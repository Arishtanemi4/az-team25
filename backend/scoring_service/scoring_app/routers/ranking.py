"""POST /rank -- streams Server-Sent Events: one `stage` event per pipeline checkpoint, then one
`result` event with the full output contract (or one `error` event on a failure discovered after
streaming has begun)."""

import json
from typing import Any

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import StreamingResponse

from scoring_app.schemas.ranking import RankRequest
from scoring_app.services.ranking_service import TOTAL_STAGES

router = APIRouter(tags=["ranking"])


def _sse_event(event: str, data: dict[str, Any]) -> str:
    return f"event: {event}\ndata: {json.dumps(data)}\n\n"


@router.post("/rank")
def rank(request: Request, body: RankRequest) -> StreamingResponse:
    service = request.app.state.ranking_service
    try:
        # Resolved synchronously, before any SSE framing starts, so an ambiguous/unresolved
        # gene still surfaces as a normal 422 -- the status code can't change once streaming
        # headers are sent.
        service.check_genes(body.inclusion_genes, body.exclusion_genes)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    def event_stream():
        stage_index = 0
        try:
            for kind, payload in service.rank_stream(
                body.inclusion_genes,
                body.exclusion_genes,
                lineage=body.lineage,
                primary_disease=body.primary_disease,
                exclude_problematic=body.exclude_problematic,
                msi_high=body.msi_high,
                ploidy_min=body.ploidy_min,
                ploidy_max=body.ploidy_max,
                require_metabolomics=body.require_metabolomics,
                require_mirna=body.require_mirna,
                top_k=body.top_k,
            ):
                if kind == "stage":
                    stage_index += 1
                    yield _sse_event("stage", {"label": payload, "index": stage_index, "total": TOTAL_STAGES})
                else:  # kind == "result"
                    yield _sse_event("result", payload)
        except ValueError as exc:
            yield _sse_event("error", {"message": str(exc)})
        except Exception:
            yield _sse_event("error", {"message": "Ranking failed unexpectedly. Check scoring_service logs."})

    return StreamingResponse(event_stream(), media_type="text/event-stream")
