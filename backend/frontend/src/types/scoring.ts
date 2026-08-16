// Mirrors scoring/export.py::export_query_result. The backend intentionally has no parallel
// response schema, so this client contract names native fields and renders unknown future states
// explicitly instead of treating them as a score or a safe default.

export type LayerName = "rna" | "protein" | "fusion" | "copy_number" | "mutation" | "dependency";
export type KnownLayerState = "measured" | "measured_absent" | "not_assayed" | "non_detected" | "excluded_pan_essential";
export type LayerState = KnownLayerState | (string & {});
export type KnownConfidenceTier = "High" | "Moderate" | "Low" | "Insufficient";
export type ConfidenceTier = KnownConfidenceTier | (string & {});
export type GeneRole = "inclusion" | "exclusion" | (string & {});

export interface LayerDetail {
  value: number | string | null;
  unit: string;
  d: number | null;
  state: LayerState;
  confidence_factor?: number;
  // gene_class: only present on the mutation/copy_number layers. pan_essential: only present on
  // the dependency layer. Both optional here because all six layers share this one shape.
  gene_class?: string;
  pan_essential?: boolean;
}

export interface PerGeneResult {
  ensembl_id: string;
  symbol: string;
  role: GeneRole;
  d_gene: number | null;
  layers: Record<LayerName, LayerDetail>;
  missing_layers?: LayerName[];
  narrative: string | null;
  // True exactly when scoring/score.py zeroed the fusion layer's weight for this (exclusion-role)
  // gene's measured_absent fusion reading (V6-2c, CONSTRAINTS.md S3/U3). Optional so older
  // exported JSON fixtures/snapshots without the field still type-check.
  fusion_weight_zeroed?: boolean;
}

export interface VetoInfo {
  ensembl_id: string;
  symbol: string;
  reason: "exclusion_expressed" | "inclusion_below_floor" | (string & {});
}

export interface MissingEvidenceEntry {
  ensembl_id: string;
  symbol?: string;
  role?: GeneRole;
  layer?: LayerName | string | null;
  severity: "cannot_certify_absence" | "untested" | "gene_role_unknown" | "measured_unresolvable" | (string & {});
}

export interface WarningEntry {
  type: string;
  message: string;
  ensembl_id?: string;
  symbol?: string;
}

export interface FusionContextEntry {
  ensembl_id: string;
  symbol: string;
  events: number;
}

export interface SimilarLine {
  model_id: string;
  correlation: number;
  D: number | null;
}

export interface RankedCellLine {
  model_id: string;
  // Populated by backend/scoring_service/.../ranking_service.py from data/processed/cell_lines.csv,
  // joined on ModelID -- display metadata only, never a scored input (V7-1c). Absent/None when
  // the ModelID has no name on record; render the bare ModelID in that case.
  cell_line_name?: string | null;
  D: number | null;
  confidence_tier: ConfidenceTier;
  hpa_agreement: boolean | null;
  geo_agreement: boolean | null;
  veto: VetoInfo | null;
  rho_bar: Record<string, number | null>;
  m_eff: number;
  m: number;
  warnings: WarningEntry[];
  per_gene: PerGeneResult[];
  missing_evidence: MissingEvidenceEntry[];
  is_problematic: boolean;
  tumour_representativeness: null;
  tumour_representativeness_note: string;
  hallmark_tags: null;
  hallmark_tags_note: string;
  fusion_context: FusionContextEntry[];
  similar_lines: SimilarLine[];
}

export interface CandidateCounts {
  evaluated: number;
  eligible: number;
  displayed: number;
  ranked_beyond_top_n: number;
  low_confidence: number;
  insufficient: number;
  disqualified: number;
}

export interface DiagnosticInfo {
  message: string;
  most_disqualifying_gene: [string, string] | null;
  veto_counts: Record<string, number>;
  veto_counts_by_gene?: Array<{ ensembl_id: string; symbol: string | null; count: number }>;
}

export interface RankResponse {
  schema_version: string;
  query: {
    inclusion_genes: string[];
    exclusion_genes: string[];
    filters: Record<string, string | number | boolean | null>;
  };
  ranked_cell_lines: RankedCellLine[];
  ranked_beyond_top_n: RankedCellLine[];
  low_confidence_lines: RankedCellLine[];
  insufficient_evidence_lines: RankedCellLine[];
  disqualified_lines: RankedCellLine[];
  total_ranked: number;
  candidate_counts: CandidateCounts;
  diagnostic: DiagnosticInfo | null;
  boundary_statement: string;
  // Present only when the S09 local extension API is enabled; native ranking remains usable when
  // this is absent or unavailable.
  research_query_id?: string | null;
  research_status?: "ready" | "unavailable" | string;
  research_reason?: string;
  total_candidates_scored?: number;
}

export interface GeneMatch {
  ensembl_id: string;
  symbol: string;
}
