"""POST /expand -- suggests genes related to the researcher's inclusion genes, with knowledge-
graph provenance (rag/query_expansion_agent.py). Suggestions are proposals only; nothing here
changes a query or a score."""

from typing import Any

from fastapi import APIRouter, HTTPException, Request

from rag_app.lib import rag_path  # noqa: F401 -- import order matters: puts rag/ on sys.path
                                    # before the QueryExpansionFailedError import below can succeed.

from query_expansion_agent import QueryExpansionFailedError

from rag_app.schemas.rag import QueryExpansionRequest

router = APIRouter(tags=["query_expansion"])


@router.post("/expand")
def expand(request: Request, body: QueryExpansionRequest) -> dict[str, Any]:
    service = request.app.state.rag_service
    try:
        return service.expand_query(body.inclusion_genes)
    except QueryExpansionFailedError as exc:
        # The model exhausted its tool-call turns without a parseable final answer -- a genuine
        # agent failure, not a bad request. Caught before the generic RuntimeError clause below,
        # since QueryExpansionFailedError is itself a RuntimeError subclass.
        raise HTTPException(status_code=502, detail=str(exc)) from exc
    except RuntimeError as exc:
        # e.g. NVIDIA_API_KEY is not set -- a service misconfiguration, not a bad request.
        raise HTTPException(status_code=500, detail=str(exc)) from exc
