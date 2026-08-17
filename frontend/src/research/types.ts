import type { LayerDetail, WarningEntry } from "../types/scoring";

export interface EventPacket {
  status: string;
  reason?: string | null;
  rows?: Array<Record<string, unknown>>;
  total_count?: number;
  truncation_note?: string;
}

export interface ContextGene {
  ensembl_id: string;
  symbol: string;
  role: string;
  gene_class: string;
  d_gene: number | null;
  layers: Record<string, LayerDetail & { unit_description?: string; caveats?: Array<{ type: string; message: string }>; events?: EventPacket }>;
  narrative: string | null;
  evidence_caveats: Array<{ type: string; message: string }>;
}

export interface MeasuredContext {
  schema_version: "measured-context-v1" | string;
  query_id: string;
  model_id: string;
  partition_bucket: string;
  scoring_snapshot: {
    D: number | null;
    confidence_tier: string;
    veto: unknown;
    hpa_agreement: boolean | null;
    geo_agreement: boolean | null;
    warnings: WarningEntry[];
  };
  genes: ContextGene[];
  model_context: {
    is_problematic: boolean | null;
    tumour_representativeness_note?: string | null;
    hallmark_tags_note?: string | null;
    supplementary_assay_availability?: Record<string, { status: string; reason?: string | null; available?: boolean | null }>;
  };
  boundary_statement?: string | null;
}

export interface GraphEdge {
  neighbour_id: string;
  neighbour_node_id?: string | null;
  neighbour_type: string;
  display_name?: string | null;
  edge_type?: string | null;
  source_db?: string | null;
  weight?: number | null;
  evidence_detail?: string | null;
  url?: string | null;
}

export interface BoundedGraphEntries {
  displayed: GraphEdge[];
  total_count: number;
  truncated: boolean;
}

export interface GraphGeneContext {
  role: string;
  claim_type: string;
  reason?: string | null;
  direct_edges: BoundedGraphEntries;
  pathway_membership: BoundedGraphEntries;
}

export interface GraphContext {
  schema_version: "graph-context-v1" | string;
  graph_source: { mode: string; reason?: string | null; scope?: string | null };
  genes: Record<string, GraphGeneContext>;
  shared_pathways: Record<string, { pathway_ids: string[]; count: number }>;
  limitations: Array<{ claim_type?: string; message: string }>;
}

export interface ContextResponse {
  query_id: string;
  model_id: string;
  model_metadata: Record<string, unknown> | null;
  measured_context: MeasuredContext;
  graph_context: GraphContext;
}

export interface AlternativeViewScore {
  status: string;
  score?: number;
  weight?: number;
  joint_gene_fraction?: number;
  reason?: string | null;
}

export interface AlternativeCandidate {
  model_id: string;
  rank: number;
  S: number;
  native_D: number | null;
  native_confidence_tier: string;
  view_scores: Record<string, AlternativeViewScore>;
  missing_views: string[];
  available_weight: number;
  available_weight_fraction: number;
  small_cohort_caveats: Array<{ message: string; view_id?: string; gene_id?: string; observed_n?: number }>;
  one_gene_caveat: boolean;
  leave_one_view_out_S_range: { minimum: number | null; maximum: number | null; by_omitted_view: Record<string, number | null> };
}

export interface AlternativesResponse {
  schema_version: string;
  query_id: string;
  anchor_id: string;
  method: { version?: string; label: string; is_probability: boolean; boundary: string };
  reference_cohort: { count?: number; model_ids?: string[]; hash?: string | null };
  candidate_counts: { considered: number; qualifying: number; returned: number };
  status: string;
  reason?: string | null;
  top_alternatives: AlternativeCandidate[];
  limitations: string[];
}

export interface ComparisonRow {
  row_id: string;
  section: string;
  label: string;
  values: Record<string, unknown>;
  gene_id?: string;
  gene_role?: string;
  layer?: string;
  unit?: string;
  provenance?: string;
}

export interface ComparisonResponse {
  schema_version: string;
  query_id: string;
  selected_model_ids: string[];
  rows: ComparisonRow[];
  limitations: string[];
  boundary: string;
}
