"""Opt-in FastAPI routes for local, snapshot-backed research extensions (S09)."""

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import Response
from pydantic import BaseModel, ConfigDict

import service
from serialization import to_json_safe

router = APIRouter(prefix="/research", tags=["research"])


class ComparisonRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    model_ids: list[str]


class ComparisonExportRequest(ComparisonRequest):
    format: str = "csv"


class RelationshipDeclaration(BaseModel):
    """Only an explicit relationship declaration can start an experimental route evaluation."""
    model_config = ConfigDict(extra="forbid")
    type: str
    source_gene: str
    target_gene: str
    context_type: str | None = None


class RouteEvaluationRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    relationship: RelationshipDeclaration


class AutomaticAdvanceRequest(BaseModel):
    """No client route choice is accepted for the fixed automatic battery."""
    model_config = ConfigDict(extra="forbid")


def _service(request):
    return request.app.state.research_service


def _call(operation):
    try:
        result = operation()
    except service.ResearchNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except service.ResearchStateError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    # Downloads (_download) already return a Response with pre-formatted bytes; every other
    # endpoint returns a plain dict straight from the service layer, which can carry NaN from
    # real graph/expression data -- to_json_safe is the module's own public JSON-safety boundary
    # (serialization.py), applied once here rather than scattered across each service method.
    return result if isinstance(result, Response) else to_json_safe(result)


def _download(result):
    return Response(content=result["content"], media_type=result["media_type"])


@router.get("/queries/{query_id}/models/{model_id}/context")
def get_context(query_id: str, model_id: str, request: Request):
    return _call(lambda: _service(request).get_context(query_id, model_id))


@router.get("/queries/{query_id}/models/{model_id}/alternatives")
def get_alternatives(query_id: str, model_id: str, request: Request):
    return _call(lambda: _service(request).get_alternatives(query_id, model_id))


@router.post("/queries/{query_id}/compare")
def compare_models(query_id: str, body: ComparisonRequest, request: Request):
    return _call(lambda: _service(request).compare(query_id, body.model_ids))


@router.get("/queries/{query_id}/export")
def export_latest(query_id: str, request: Request, format: str = "json"):
    return _call(lambda: _download(_service(request).export_latest(query_id, format)))


@router.post("/queries/{query_id}/export-comparison")
def export_comparison(query_id: str, body: ComparisonExportRequest, request: Request):
    return _call(lambda: _download(_service(request).export_comparison(query_id, body.model_ids, body.format)))


@router.post("/queries/{query_id}/routes/evaluate")
def evaluate_route(query_id: str, body: RouteEvaluationRequest, request: Request):
    relationship = body.relationship
    return _call(lambda: _service(request).evaluate_route(
        query_id, relationship.type, relationship.source_gene, relationship.target_gene,
        relationship.context_type,
    ))


@router.get("/queries/{query_id}/routes")
def get_routes(query_id: str, request: Request):
    return _call(lambda: _service(request).get_routes(query_id))


@router.get("/queries/{query_id}/routes/automatic")
def get_automatic_routes(query_id: str, request: Request):
    """Read automatic-route state only; it cannot create a plan or run an evaluator."""
    return _call(lambda: _service(request).get_automatic_routes(query_id))


@router.post("/queries/{query_id}/routes/automatic/advance")
def advance_automatic_routes(query_id: str, request: Request, body: AutomaticAdvanceRequest | None = None):
    """Advance the fixed automatic plan by at most one persisted pending task."""
    return _call(lambda: _service(request).advance_automatic_routes(query_id))


@router.get("/queries/{query_id}/models/{model_id}/routes")
def get_model_routes(query_id: str, model_id: str, request: Request):
    return _call(lambda: _service(request).get_model_routes(query_id, model_id))


@router.get("/queries/{query_id}/routes/export")
def export_routes(query_id: str, request: Request, format: str = "json"):
    return _call(lambda: _download(_service(request).export_routes(query_id, format)))
