"""GET /graph/neighborhood -- seed-gene k-hop neighbourhood from the knowledge graph
(PRODUCT_SURFACE.md SS3.5). Explanation only; never scores, never feeds a query (F4)."""

from fastapi import APIRouter, Query

from rag_app.services import graph_service

router = APIRouter(tags=["graph"])


@router.get("/graph/neighborhood")
def neighborhood(
    gene: str = Query(...),
    hops: int = Query(1, ge=1, le=2),
    max_neighbors: int = Query(50, ge=1, le=200),
):
    # available:false is a first-class response shape, same precedent as
    # data_app/routers/eda.py::concordance -- no single response_model on purpose.
    return graph_service.get_neighborhood(gene, hops, max_neighbors)
