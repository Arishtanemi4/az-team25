// Mirrors the JSON contracts backend/rag_service's four routers actually return (each router
// returns rag/'s own entry-point dict verbatim -- see backend/rag_service/app/schemas/rag.py).
// Kept in sync by hand -- there is deliberately no second backend schema to drift from, same
// convention as types/scoring.ts.

export interface VerificationResult {
  passed: boolean;
  numeric_failures: string[];
  citation_failures: string[];
}

// --- Result narrator (POST /narrate) ---

export interface GeneExplanation {
  ensembl_id: string;
  symbol: string;
  explanation: string;
}

export interface CellLineExplanation {
  model_id: string;
  summary: string;
  gene_explanations: GeneExplanation[];
  limitations: string;
}

export interface NarrationResult {
  overview: string;
  cell_lines: CellLineExplanation[];
  boundary_statement: string;
}

export interface NarrationResponse {
  narration: NarrationResult;
  text: string;
  verification: VerificationResult;
}

// --- Methodology Q&A agent (POST /methodology/ask) ---

export interface MethodologyAnswer {
  answer: string;
  context: string;
  abstained: boolean;
  verification: VerificationResult | null;
}

// --- Query-expansion agent (POST /expand) ---

export interface ExpansionSuggestion {
  gene: string;
  reason: string;
  source_db: string | null;
  edge_type: string | null;
  score: number | null;
}

export interface ExpansionResponse {
  suggestions: ExpansionSuggestion[];
  missing_provenance: ExpansionSuggestion[];
  verification: VerificationResult;
}

// --- Literature agent (POST /literature/search) ---

export interface LiteratureFinding {
  pmid: string;
  quote: string;
  relevance: string;
}

export interface LiteratureResponse {
  findings: LiteratureFinding[];
  n_reported_by_model: number;
  n_quote_confirmed: number;
}
