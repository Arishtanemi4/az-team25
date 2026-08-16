"""Source-labelled graph and pathway-membership context (extensions/plans/CONTRACTS.md C4).

`build_graph_context(resolved_genes, roles, config)` answers, generically, "what does the graph
say about these genes" -- for every requested gene, including exclusion-role genes and a
single-gene query. It never reads or writes an S03 snapshot or an S04 measured-context packet;
it is combined with the measured packet by a later stage (S08), without altering either.

Two data sources, tried in this order, never both required:
1. The full built graph (rag/knowledge_graph.py's public load_graph/get_neighborhood, explicit
   paths only -- Reactome + STRING + BioGRID + fusion + co-dependency edges).
2. reactome_membership.py's pathway-membership-only fallback, parsed from the raw Reactome
   files directly -- explicitly partial, no gene-gene edges.
Neither available -> every gene gets `claim_type: "missing_context"` and the software is still
exercised by fixtures; a real R4 demonstration needs one of the two sources present on disk.

C4's governing rule, restated here because it binds this module specifically: membership is not
activation, enrichment, causality, or sensitivity. Every edge/pathway record below is a generic
reference fact about the two entities it connects, never a claim about any specific cell line --
which is also why this module never touches a model_id.
"""

import sys
from pathlib import Path

_RESEARCH_DIR = str(Path(__file__).resolve().parent)
if _RESEARCH_DIR not in sys.path:
    sys.path.insert(0, _RESEARCH_DIR)

_RAG_DIR = str(Path(__file__).resolve().parents[2] / "rag")
if _RAG_DIR not in sys.path:
    sys.path.insert(0, _RAG_DIR)

import reactome_membership  # noqa: E402 -- sys.path must be set up first
import knowledge_graph as rag_knowledge_graph  # noqa: E402 -- rag/knowledge_graph.py, public API only

MAX_DISPLAYED_EDGES = 10
MAX_DISPLAYED_PATHWAYS = 10
_GENE_NODE_PREFIX = "gene:"
_PATHWAY_NODE_PREFIX = "pathway:"


def _bounded(items, cap):
    """C4: 'UI cap 10 edges/10 pathways per gene with total/truncation counts; export all
    retrieved direct records.' The full list is always returned as `total_count`; only the
    *displayed* slice is capped -- nothing retrieved is ever discarded from the record."""
    return {
        "displayed": items[:cap],
        "total_count": len(items),
        "truncated": len(items) > cap,
    }


def _try_load_full_graph(config):
    """Primary path. Returns (graph, None) on success, (None, reason) on any failure -- a
    missing or corrupt graph is reported, never silently treated as 'no gene has context'."""
    config = config or {}
    nodes_path = config.get("graph_nodes_path", rag_knowledge_graph.NODES_PATH)
    edges_path = config.get("graph_edges_path", rag_knowledge_graph.EDGES_PATH)
    try:
        graph = rag_knowledge_graph.load_graph(nodes_path=nodes_path, edges_path=edges_path)
    except FileNotFoundError as exc:
        return None, f"Full graph Parquet not found: {exc}"
    except Exception as exc:  # noqa: BLE001 -- reported, never raised past this function
        return None, f"Full graph failed to load: {exc}"
    return graph, None


def _graph_gene_node_id(ensembl_id):
    """The frozen graph's public node-table schema namespaces canonical Ensembl IDs.

    Callers of this extension supply canonical IDs (for example ``ENSG00000146648``), while
    ``rag/build_knowledge_graph.py::build_nodes_table`` persists those as
    ``gene:ENSG00000146648``.  Keep that storage detail at this one boundary rather than
    leaking it into the S01--S04 identifier contracts.
    """
    return ensembl_id if str(ensembl_id).startswith(_GENE_NODE_PREFIX) else f"{_GENE_NODE_PREFIX}{ensembl_id}"


def _source_entity_id(node_id, node_type):
    """Return the source identifier while retaining the exact graph node ID separately."""
    node_id = str(node_id)
    if node_type == "gene" and node_id.startswith(_GENE_NODE_PREFIX):
        return node_id[len(_GENE_NODE_PREFIX):]
    if node_type == "pathway" and node_id.startswith(_PATHWAY_NODE_PREFIX):
        return node_id[len(_PATHWAY_NODE_PREFIX):]
    return node_id


def _flat_edges_for_gene(graph, ensembl_id):
    """Wraps rag/knowledge_graph.py's public get_neighborhood(graph, node_id, hops=1) -- 1-hop
    neighbours ARE the direct edges C4 asks for; no separate traversal is written here. Flattens
    from 'one record per neighbour node' to 'one record per (neighbour, edge)' pair, since a
    MultiGraph can carry more than one edge to the same neighbour (e.g. both a STRING and a
    BioGRID edge) and C4 requires those preserved as independent records, not collapsed.
    Sorted by (source_db, edge_type, neighbour node_id) -- C4's 'stable ordering: source, type,
    neighbour ID.'"""
    graph_gene_node_id = _graph_gene_node_id(ensembl_id)
    if graph_gene_node_id not in graph:
        return None  # distinct from "measured, zero edges" -- the gene isn't a node at all
    neighbours = rag_knowledge_graph.get_neighborhood(graph, graph_gene_node_id, hops=1)
    flat = []
    for neighbour in neighbours:
        # get_neighborhood deliberately gives the public direct-neighbour traversal. Its compact
        # edge summary omits evidence_detail, so recover every parallel edge's full persisted
        # attributes from the already loaded public graph object -- no private builder import or
        # second traversal is involved. This preserves C4 provenance details as well as source.
        parallel_edges = graph.get_edge_data(graph_gene_node_id, neighbour["node_id"], default={})
        for edge in parallel_edges.values():
            flat.append(
                {
                    "neighbour_id": _source_entity_id(neighbour["node_id"], neighbour["node_type"]),
                    "neighbour_node_id": neighbour["node_id"],
                    "neighbour_type": neighbour["node_type"],
                    "display_name": neighbour["display_name"],
                    "edge_type": edge.get("edge_type"),
                    "source_db": edge.get("source_db"),
                    "weight": edge.get("weight"),
                    "evidence_detail": edge.get("evidence_detail"),
                }
            )
    flat.sort(
        key=lambda e: (
            e["source_db"] or "",
            e["edge_type"] or "",
            e["neighbour_id"],
            e["neighbour_node_id"],
            str(e["evidence_detail"] or ""),
        )
    )
    return flat


def _gene_context_from_full_graph(graph, ensembl_id, role):
    """Returns (entry, full_pathway_ids) -- the second element is the *untruncated* set of
    pathway neighbour IDs, kept only long enough for _shared_pathways to use below; it is never
    part of the returned entry itself."""
    flat_edges = _flat_edges_for_gene(graph, ensembl_id)
    if flat_edges is None:
        entry = {
            "role": role,
            "claim_type": "missing_context",
            "reason": "Gene is not present in the graph (no node for this ensembl_id).",
            "direct_edges": _bounded([], MAX_DISPLAYED_EDGES),
            "pathway_membership": _bounded([], MAX_DISPLAYED_PATHWAYS),
        }
        return entry, set()
    pathway_edges = [e for e in flat_edges if e["neighbour_type"] == "pathway"]
    entry = {
        "role": role,
        "claim_type": "reference_association",
        "reason": None,
        "direct_edges": _bounded(flat_edges, MAX_DISPLAYED_EDGES),
        "pathway_membership": _bounded(pathway_edges, MAX_DISPLAYED_PATHWAYS),
    }
    return entry, {e["neighbour_id"] for e in pathway_edges}


def _gene_context_from_membership_only(memberships_by_gene, ensembl_id, role):
    """Returns (entry, full_pathway_ids), same convention as the full-graph builder above."""
    memberships = memberships_by_gene.get(ensembl_id, [])
    pathway_edges = [
        {
            "neighbour_id": m["reactome_id"],
            # Membership-only mode has no materialised graph node table. Keep the field present
            # for a stable edge-record shape, but do not manufacture a graph-node identifier.
            "neighbour_node_id": None,
            "neighbour_type": "pathway",
            "display_name": m["pathway_name"],
            "edge_type": "member_of_pathway",
            "source_db": "Reactome",
            "weight": None,
            "evidence_detail": None,
        }
        for m in memberships
    ]
    entry = {
        "role": role,
        "claim_type": "reference_association" if pathway_edges else "missing_context",
        "reason": None if pathway_edges else "No Reactome pathway membership found for this gene.",
        "direct_edges": _bounded([], MAX_DISPLAYED_EDGES),  # membership-only: no gene-gene edges
        "pathway_membership": _bounded(pathway_edges, MAX_DISPLAYED_PATHWAYS),
    }
    return entry, {e["neighbour_id"] for e in pathway_edges}


def _shared_pathways(full_pathway_ids_by_gene):
    """Pairwise membership intersections -- C4: 'Shared pathways are membership intersections,
    not a significance test.' Computed over the FULL retrieved pathway set for each gene (the
    untruncated sets passed in here), never the display-capped slice, so a shared pathway just
    past the top-10 cap is never missed."""
    gene_ids = sorted(full_pathway_ids_by_gene)
    shared = {}
    for i, gene_a in enumerate(gene_ids):
        for gene_b in gene_ids[i + 1:]:
            common = sorted(full_pathway_ids_by_gene[gene_a] & full_pathway_ids_by_gene[gene_b])
            if common:
                shared[f"{gene_a}|{gene_b}"] = {"pathway_ids": common, "count": len(common)}
    return shared


def build_graph_context(resolved_genes, roles, config=None):
    """Assembles graph context for every gene in `resolved_genes` (a list of canonical
    ensembl_ids -- C1), including exclusion-role genes and a single-gene query. `roles` is a
    `{ensembl_id: "inclusion"|"exclusion"}` mapping covering the same genes. `config` is an
    optional dict of path overrides (`graph_nodes_path`, `graph_edges_path`, `reactome_dir`,
    `reactome_cache_dir`); every key defaults to the same paths the baseline itself uses.

    No LLM and no live network are ever used (C4, S05's own scope note)."""
    config = config or {}
    graph, graph_failure_reason = _try_load_full_graph(config)

    memberships_by_gene = None
    membership_failure_reason = None
    if graph is None:
        try:
            memberships_by_gene = reactome_membership.gene_pathway_memberships(
                resolved_genes,
                reactome_dir=config.get("reactome_dir"),
                cache_dir=config.get("reactome_cache_dir"),
            )
        except Exception as exc:  # noqa: BLE001 -- corrupt local fallback data is unavailable, not no context
            membership_failure_reason = str(exc)

    if graph is not None:
        mode, source_reason, scope = "full_graph", None, (
            "Reactome + STRING + BioGRID + project-derived fusion/co-dependency edges"
        )
    elif memberships_by_gene is not None:
        mode, source_reason, scope = "reactome_membership_only", None, "Reactome pathway membership only"
    else:
        mode, source_reason, scope = "unavailable", (
            f"Full graph unavailable ({graph_failure_reason}); Reactome membership fallback "
            f"also unavailable ({membership_failure_reason})."
        ), None

    genes_context = {}
    full_pathway_ids_by_gene = {}
    for ensembl_id in resolved_genes:
        role = roles.get(ensembl_id, "unknown")
        if mode == "full_graph":
            entry, full_pathway_ids = _gene_context_from_full_graph(graph, ensembl_id, role)
        elif mode == "reactome_membership_only":
            entry, full_pathway_ids = _gene_context_from_membership_only(
                memberships_by_gene, ensembl_id, role
            )
        else:
            entry, full_pathway_ids = {
                "role": role,
                "claim_type": "missing_context",
                "reason": "No graph source available (see graph_source.reason).",
                "direct_edges": _bounded([], MAX_DISPLAYED_EDGES),
                "pathway_membership": _bounded([], MAX_DISPLAYED_PATHWAYS),
            }, set()
        genes_context[ensembl_id] = entry
        full_pathway_ids_by_gene[ensembl_id] = full_pathway_ids

    shared_pathways = _shared_pathways(full_pathway_ids_by_gene)

    limitations = []
    if mode == "reactome_membership_only":
        limitations.append({
            "claim_type": "limitation",
            "message": (
                "Graph context for this query used the Reactome-membership-only fallback: "
                "pathway membership is available, but gene-gene interaction/association edges "
                "(STRING, BioGRID, fusion, co-dependency) are not, since the full built graph "
                "was unavailable."
            ),
        })
    elif mode == "unavailable":
        limitations.append({"claim_type": "limitation", "message": source_reason})

    return {
        "schema_version": "graph-context-v1",
        "graph_source": {"mode": mode, "reason": source_reason, "scope": scope},
        "genes": genes_context,
        "shared_pathways": shared_pathways,
        "limitations": limitations,
    }
