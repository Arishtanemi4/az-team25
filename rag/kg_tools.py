import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT / "scoring") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scoring"))

import knowledge_graph
import score as scoring_score  # scoring/score.py -- resolve_genes only, see module docstring

GENE_REFERENCE_PATH = REPO_ROOT / "data" / "processed" / "gene_reference.csv"

_graph_cache = None
_gene_reference_cache = None


def _get_graph():
    global _graph_cache
    if _graph_cache is None:
        _graph_cache = knowledge_graph.load_graph()
    return _graph_cache


def _get_gene_reference():
    global _gene_reference_cache
    if _gene_reference_cache is None:
        _gene_reference_cache = pd.read_csv(GENE_REFERENCE_PATH)
    return _gene_reference_cache


def resolve_symbol(token):
    resolved, ambiguous, unresolved = scoring_score.resolve_genes([token], _get_gene_reference())
    if token in resolved:
        return {"token": token, "ensembl_id": resolved[token], "status": "resolved"}
    if token in ambiguous:
        return {"token": token, "ensembl_id": None, "status": "ambiguous"}
    return {"token": token, "ensembl_id": None, "status": "unresolved"}


def _node_id_for_gene(token):
    resolution = resolve_symbol(token)
    if resolution["status"] != "resolved":
        return None, resolution
    return f"gene:{resolution['ensembl_id']}", resolution


DEFAULT_MAX_NEIGHBORS = 25


def _best_edge_weight(neighbor):
    weights = [e.get("weight") for e in neighbor.get("edges", [])]
    numeric = [w for w in weights if isinstance(w, (int, float)) and w == w]  # w == w excludes NaN
    return max(numeric) if numeric else float("-inf")


def kg_neighbors(gene, k=1, max_neighbors=DEFAULT_MAX_NEIGHBORS):
    node_id, resolution = _node_id_for_gene(gene)
    if node_id is None:
        return {"resolution": resolution, "neighbors": [], "total_found": 0, "truncated": False}

    neighbors = knowledge_graph.get_neighborhood(_get_graph(), node_id, hops=k)
    total_found = len(neighbors)
    if max_neighbors is not None and total_found > max_neighbors:
        neighbors = sorted(neighbors, key=_best_edge_weight, reverse=True)[:max_neighbors]
        truncated = True
    else:
        truncated = False
    return {
        "resolution": resolution, "neighbors": neighbors,
        "total_found": total_found, "truncated": truncated,
    }


def shared_pathways(genes):
    pathway_sets = []
    pathway_details = {}
    for gene in genes:
        node_id, _ = _node_id_for_gene(gene)
        if node_id is None:
            return {"shared": [], "note": f"'{gene}' did not resolve to a gene in the graph"}
        neighbors = knowledge_graph.get_neighborhood(_get_graph(), node_id, hops=1)
        pathway_ids = {n["node_id"] for n in neighbors if n["node_type"] == "pathway"}
        pathway_sets.append(pathway_ids)
        for n in neighbors:
            if n["node_type"] == "pathway":
                pathway_details[n["node_id"]] = n["display_name"]

    if not pathway_sets:
        return {"shared": [], "note": "no genes supplied"}
    shared_ids = set.intersection(*pathway_sets)
    return {"shared": [{"node_id": pid, "display_name": pathway_details[pid]} for pid in shared_ids]}


def interaction_evidence(gene_a, gene_b):
    node_a, resolution_a = _node_id_for_gene(gene_a)
    node_b, resolution_b = _node_id_for_gene(gene_b)
    if node_a is None or node_b is None:
        return {"resolutions": [resolution_a, resolution_b], "edges": []}

    graph = _get_graph()
    if not graph.has_edge(node_a, node_b):
        return {"resolutions": [resolution_a, resolution_b], "edges": []}

    edges = [
        {"edge_type": attrs.get("edge_type"), "source_db": attrs.get("source_db"), "weight": attrs.get("weight")}
        for attrs in graph.get_edge_data(node_a, node_b).values()
    ]
    return {"resolutions": [resolution_a, resolution_b], "edges": edges}
