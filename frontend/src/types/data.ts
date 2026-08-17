// Mirrors backend/data_service/data_app/schemas/tables.py field-for-field. `role` and
// `cardinality` stay `string` rather than a frontend union -- they are server-controlled
// vocabulary (table_registry.py), the same reasoning scoring.ts's `(string & {})` escape
// hatch already uses for state strings.

export interface ColumnInfo {
  name: string;
  role: string; // join_key | value | state | context
  dtype: string;
}

export interface JoinInfo {
  to: string;
  on: string; // ModelID | ensembl_id -- the only two automatic join keys
  cardinality: string;
}

export interface TableInfo {
  name: string;
  grain: string;
  row_count: number;
  columns: ColumnInfo[];
  joins: JoinInfo[];
}

export interface SchemaResponse {
  tables: TableInfo[];
}

export interface TableSummary {
  name: string;
  row_count: number;
}

export interface TablesResponse {
  tables: TableSummary[];
}

export interface PreviewResponse {
  table: string;
  grain: string;
  total_rows: number;
  returned_rows: number;
  columns: string[];
  rows: Record<string, unknown>[];
}

// --- V7-6: mirrors backend/data_service/data_app/schemas/{relationships,eda}.py field-for-field ---

export type RelationshipLayer = "rna" | "rna_hpa" | "rna_geo" | "protein" | "dependency" | "copy_number";

export interface GeneCellLinePoint {
  model_id: string;
  cell_line_name: string | null;
  lineage: string | null;
  value: number;
}

export interface GeneCellLineResponse {
  ensembl_id: string;
  symbol: string | null;
  layer: string;
  unit: string;
  n_measured: number;
  n_models_total: number;
  points: GeneCellLinePoint[];
}

export interface CellLineDiseaseGroup {
  lineage: string;
  primary_disease: string;
  n_models: number;
}

export interface CellLineDiseaseResponse {
  n_models_total: number;
  cells: CellLineDiseaseGroup[];
}

export interface GeneDiseaseGroup {
  primary_disease: string;
  n_measured: number;
  n_models: number;
  median: number | null;
  iqr_low: number | null;
  iqr_high: number | null;
  flagged: boolean;
}

export interface GeneDiseaseResponse {
  ensembl_id: string;
  symbol: string | null;
  layer: string;
  min_n: number;
  groups: GeneDiseaseGroup[];
}

export interface CoverageResponse {
  coverage_state_counts: Record<string, Record<string, number>>;
}

export interface LineageCoverageRow {
  lineage: string;
  layer: string;
  n: number;
  n_measured: number;
  measured_fraction: number;
  flagged: boolean;
}

export interface LineageCoverageResponse {
  min_lineage_n: number;
  rows: LineageCoverageRow[];
}

// The single leaf shape build_eda_aggregates.py writes for every source pair -- all fields are
// null/0 when a pair's model intersection was empty.
export interface ConcordanceSummary {
  n_models_scored: number;
  n_models_in_intersection: number;
  median_rho: number | null;
  iqr_low: number | null;
  iqr_high: number | null;
  median_n_genes_used: number | null;
  min_genes_floor: number;
}

export interface EdaConcordanceAvailable {
  available: true;
  generated_utc: string;
  rna_cross_source: Record<string, ConcordanceSummary>;
  rna_protein: Record<string, ConcordanceSummary>;
}

export interface EdaConcordanceUnavailable {
  available: false;
  reason: string;
  build_command: string;
}

// /data/eda/concordance deliberately has no single FastAPI response_model (PRODUCT_SURFACE.md
// SS3.3) -- the client must branch on `available` before reading either shape further.
export type EdaConcordanceResponse = EdaConcordanceAvailable | EdaConcordanceUnavailable;

// --- V7-7: mirrors backend/data_service/data_app/schemas/validation.py field-for-field ---

// One row of validation/study_comparability.csv, the 4-row manifest read first.
export interface ValidationStudySummary {
  study: string;
  mode: string;
  comparable_to_ranking: string; // "partial" | "no" today -- server's own vocabulary, not a closed union
  reason: string;
}

export interface ValidationStudiesResponse {
  studies: ValidationStudySummary[];
}

// rows stays a generic bag -- the 5 studies' own CSVs have genuinely different columns
// (validation/*.csv), and the backend schema is deliberately list[dict], not five separate models.
export interface ValidationStudyResponse {
  study: string;
  total_rows: number;
  returned_rows: number;
  columns: string[];
  rows: Record<string, unknown>[];
}
