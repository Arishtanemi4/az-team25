// Mirrors backend/rag_service/rag_app/services/graph_service.py's response shape
// (PRODUCT_SURFACE.md SS3.5). One flat schema, not an available/unavailable union -- `reason`
// is always present (null when not applicable), matching the artefact's own single-line
// contract: {available, seed, resolved_ensembl_id, nodes, edges, truncated, reason}.

export interface GraphNode {
  node_id: string;
  node_type: string; // "gene" | "pathway" today -- rag/build_knowledge_graph.py's own vocabulary
  display_name: string;
  hop_distance: number;
}

export interface GraphEdge {
  source: string;
  target: string;
  // exact strings rag/build_knowledge_graph.py writes:
  // gene_in_pathway | pathway_parent_of | protein_interaction | curated_interaction | fusion_partner | codependency
  edge_type: string;
  // reactome | string_v12 | biogrid | project_fusions | project_dependency
  source_db: string;
  weight: number | null;
}

export interface GraphNeighborhoodResponse {
  available: boolean;
  seed: string;
  resolved_ensembl_id: string | null;
  nodes: GraphNode[];
  edges: GraphEdge[];
  truncated: boolean;
  reason: string | null;
  build_command?: string; // present only when available is false
}
