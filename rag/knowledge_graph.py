from pathlib import Path

import networkx as nx
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
NODES_PATH = REPO_ROOT / "data" / "external" / "knowledge_graph" / "derived" / "knowledge_graph_nodes.parquet"
EDGES_PATH = REPO_ROOT / "data" / "external" / "knowledge_graph" / "derived" / "knowledge_graph_edges.parquet"


def load_graph(nodes_path=NODES_PATH, edges_path=EDGES_PATH):
    nodes_df = pd.read_parquet(nodes_path)
    edges_df = pd.read_parquet(edges_path)

    graph = nx.MultiGraph()
    for row in nodes_df.itertuples():
        graph.add_node(
            row.node_id,
            node_type=row.node_type,
            source_db=row.source_db,
            display_name=row.display_name,
        )
    for row in edges_df.itertuples():
        if row.source_node_id not in graph or row.target_node_id not in graph:
            continue  # a node absent from nodes_df is a build defect, not something to guess past
        graph.add_edge(
            row.source_node_id,
            row.target_node_id,
            edge_type=row.edge_type,
            source_db=row.source_db,
            weight=row.weight,
            evidence_detail=row.evidence_detail,
        )
    return graph


def get_neighborhood(graph, node_id, hops=1):
    if node_id not in graph:
        return []

    ego = nx.ego_graph(graph, node_id, radius=hops)
    distances = nx.single_source_shortest_path_length(graph, node_id, cutoff=hops)

    records = []
    for neighbor_id in ego.nodes:
        if neighbor_id == node_id:
            continue
        node_attrs = graph.nodes[neighbor_id]
        neighbor_distance = distances[neighbor_id]
        # Only edges that step exactly one hop closer to the seed -- e.g. for a 1-hop neighbour,
        # only its direct edge(s) to the seed itself. `ego_graph`'s induced subgraph also includes
        # *lateral* edges between two neighbours at the same distance (two of EGFR's 1-hop
        # neighbours can easily also be directly connected to each other); attributing those to
        # both neighbours would inflate "how many neighbours are STRING-connected to the seed"
        # far past the real direct-edge count -- caught via a real EGFR STRING-degree sanity check
        # during V6-7 (503 real direct edges vs. an unfiltered count in the thousands).
        inbound_edges = [
            {
                "connected_to": other,
                "edge_type": edge_attrs.get("edge_type"),
                "source_db": edge_attrs.get("source_db"),
                "weight": edge_attrs.get("weight"),
            }
            for _, other, edge_attrs in ego.edges(neighbor_id, data=True)
            if distances.get(other, hops + 1) == neighbor_distance - 1
        ]
        records.append(
            {
                "node_id": neighbor_id,
                "node_type": node_attrs.get("node_type"),
                "display_name": node_attrs.get("display_name"),
                "hop_distance": neighbor_distance,
                "edges": inbound_edges,
            }
        )

    return records
