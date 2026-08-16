"""Backs GET /graph/neighborhood -- a seed-gene k-hop neighbourhood from the knowledge graph
(PRODUCT_SURFACE.md SS3.5). Wraps rag/kg_tools.py::kg_neighbors rather than re-implementing graph
traversal; reshapes its neighbour-with-nested-edges shape into flat nodes[]/edges[] lists a
Cytoscape rendering can consume directly. Explanation only -- never scores, never calls
scoring/score_panel (F4)."""

from rag_app.lib import rag_path  # noqa: F401 -- import order matters: puts rag/ on sys.path
                                    # before the `import kg_tools` below can succeed.

import kg_tools

BUILD_COMMAND = "python rag/build_knowledge_graph.py"


def _json_safe_weight(weight):
    """BioGRID edges carry a NaN weight, not None (kg_tools.py's own _best_edge_weight already
    treats NaN as absent via `w == w`). Starlette's default JSONResponse serializes with
    allow_nan=False, so a raw NaN float crashes the response -- normalise it to None, the same
    "no weight recorded" meaning kg_tools already gives it."""
    if isinstance(weight, float) and weight != weight:
        return None
    return weight


def get_neighborhood(gene: str, hops: int, max_neighbors: int) -> dict:
    # kg_tools.kg_neighbors lazily loads and caches the graph on first call (its own
    # module-level _get_graph()) -- loading it eagerly here at service startup would crash the
    # whole rag_service if the derived Parquet pair is absent, which is exactly the
    # blank-canvas-by-crash outcome F6 forbids. A missing Parquet pair surfaces as
    # FileNotFoundError from pandas/pyarrow inside knowledge_graph.load_graph.
    try:
        result = kg_tools.kg_neighbors(gene, k=hops, max_neighbors=max_neighbors)
    except FileNotFoundError:
        return {
            "available": False, "seed": gene, "resolved_ensembl_id": None,
            "nodes": [], "edges": [], "truncated": False,
            "reason": "the knowledge graph has not been built on this machine",
            "build_command": BUILD_COMMAND,
        }

    resolution = result["resolution"]
    if resolution["status"] != "resolved":
        return {
            "available": True, "seed": gene, "resolved_ensembl_id": None,
            "nodes": [], "edges": [], "truncated": False,
            "reason": f"'{gene}' did not resolve to a gene ({resolution['status']})",
        }

    seed_node_id = f"gene:{resolution['ensembl_id']}"
    nodes = [{"node_id": seed_node_id, "node_type": "gene", "display_name": gene, "hop_distance": 0}]
    edges = []
    for neighbor in result["neighbors"]:
        nodes.append({k: neighbor[k] for k in ("node_id", "node_type", "display_name", "hop_distance")})
        for e in neighbor["edges"]:
            edges.append({
                "source": e["connected_to"], "target": neighbor["node_id"],
                "edge_type": e["edge_type"], "source_db": e["source_db"],
                "weight": _json_safe_weight(e["weight"]),
            })
    return {
        "available": True, "seed": gene, "resolved_ensembl_id": resolution["ensembl_id"],
        "nodes": nodes, "edges": edges, "truncated": result["truncated"], "reason": None,
    }
