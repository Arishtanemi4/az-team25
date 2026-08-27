"""POST /literature/search -- finds quote-verified PubMed evidence for a gene/biology claim
(rag/literature_agent.py). Every finding is a verbatim quote plus a PMID, or it is dropped."""

from typing import Any

from fastapi import APIRouter, HTTPException, Request

from rag_app.lib import rag_path  # noqa: F401 -- import order matters: puts rag/ on sys.path
                                    # before the LiteratureAgentFailedError import below can succeed.

from literature_agent import LiteratureAgentFailedError
from provider import ProviderRateLimitedError

from rag_app.schemas.rag import LiteratureQuestionRequest

router = APIRouter(tags=["literature"])


@router.post("/literature/search")
def search(request: Request, body: LiteratureQuestionRequest) -> dict[str, Any]:
    service = request.app.state.rag_service
    try:
        return service.find_literature(body.question)
    except LiteratureAgentFailedError as exc:
        # An unparseable final response or exhausted turns -- never raised for "found nothing",
        # which is a normal, honest result (CONSTRAINTS.md R4). Caught before the generic
        # RuntimeError clause below, since this is itself a RuntimeError subclass.
        raise HTTPException(
            status_code=502,
            detail=f"The literature agent couldn't produce a usable answer ({exc}). This looks "
                   "like a limitation of the current AI model on this question, not a bug -- try "
                   "rephrasing your question, or try again in a moment.",
        ) from exc
    except ProviderRateLimitedError as exc:
        raise HTTPException(status_code=429, detail=str(exc)) from exc
    except RuntimeError as exc:
        # e.g. NVIDIA_API_KEY is not set -- a service misconfiguration, not a bad request.
        raise HTTPException(status_code=500, detail=str(exc)) from exc
