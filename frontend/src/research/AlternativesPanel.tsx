import { useEffect, useState } from "react";
import { getAlternatives, ResearchRequestError } from "./client";
import type { AlternativesResponse } from "./types";
import { formatScore } from "../lib/formatScore";

interface AlternativesPanelProps {
  queryId: string | null | undefined;
  modelId: string;
  onAddToCompare: (modelId: string) => void;
}

function AlternativesContent({ result, onAddToCompare }: { result: AlternativesResponse; onAddToCompare: (modelId: string) => void }) {
  if (result.status !== "ok") return <p className="context-unavailable">{result.reason ?? "Query-profile alternatives are unavailable for this model."}</p>;
  return (
    <div className="alternatives-content">
      <p className="alternative-boundary"><strong>query-profile similarity; not experimental equivalence.</strong> S is a separate heuristic and does not replace native D or confidence tier.</p>
      <p className="context-muted">Reference cohort: {result.reference_cohort.count ?? "not reported"}; candidates considered: {result.candidate_counts.considered}; qualifying: {result.candidate_counts.qualifying}.</p>
      {result.top_alternatives.length === 0 ? <p className="context-unavailable">No qualifying same-lineage alternative was returned. No legacy RNA-only fallback is substituted.</p> : (
        <ol className="alternative-list">
          {result.top_alternatives.map((candidate) => (
            <li key={candidate.model_id}>
              <div className="alternative-header"><strong>#{candidate.rank} <code>{candidate.model_id}</code></strong><span>S = {formatScore(candidate.S)}</span><span>native D = {formatScore(candidate.native_D)}; {candidate.native_confidence_tier}</span><button type="button" onClick={() => onAddToCompare(candidate.model_id)}>Add to comparison</button></div>
              <p>Available view weight: {(candidate.available_weight_fraction * 100).toFixed(0)}%; leave-one-view-out S range: {candidate.leave_one_view_out_S_range.minimum === null ? "not available" : `${formatScore(candidate.leave_one_view_out_S_range.minimum)}–${formatScore(candidate.leave_one_view_out_S_range.maximum)}`}</p>
              <ul className="alternative-view-list">{Object.entries(candidate.view_scores).map(([view, score]) => <li key={view}>{view}: {score.status === "available" ? `S=${formatScore(score.score)}, joint coverage=${((score.joint_gene_fraction ?? 0) * 100).toFixed(0)}%` : `unavailable (${score.reason ?? "no overlap"})`}</li>)}</ul>
              {candidate.one_gene_caveat && <p className="context-muted">One-gene query: this similarity is intentionally narrow in scope.</p>}
              {candidate.small_cohort_caveats.map((caveat, index) => <p className="context-muted" key={index}>{caveat.message}</p>)}
            </li>
          ))}
        </ol>
      )}
      {result.limitations.length > 0 && <ul className="context-caveats">{result.limitations.map((limitation, index) => <li key={index}>{limitation}</li>)}</ul>}
    </div>
  );
}

export function AlternativesPanel({ queryId, modelId, onAddToCompare }: AlternativesPanelProps) {
  const [open, setOpen] = useState(false);
  const [result, setResult] = useState<AlternativesResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => { setOpen(false); setResult(null); setError(null); }, [queryId, modelId]);
  useEffect(() => {
    if (!open || !queryId) return;
    const controller = new AbortController();
    setLoading(true); setError(null);
    getAlternatives(queryId, modelId, controller.signal)
      .then((payload) => { if (payload.query_id === queryId && payload.anchor_id === modelId) setResult(payload); })
      .catch((requestError: unknown) => { if (!controller.signal.aborted) setError(requestError instanceof ResearchRequestError ? requestError.message : "Alternatives could not be loaded."); })
      .finally(() => { if (!controller.signal.aborted) setLoading(false); });
    return () => controller.abort();
  }, [open, queryId, modelId]);

  if (!queryId) return <section className="alternatives-panel"><h3>Query-profile alternatives</h3><p className="context-unavailable">A trusted research snapshot is required for alternatives.</p></section>;
  return <section className="alternatives-panel"><button type="button" className="context-toggle" onClick={() => setOpen((value) => !value)} aria-expanded={open}>{open ? "Hide query-profile alternatives" : "Show query-profile alternatives"}</button>{open && <>{loading && <p className="context-muted" role="status">Loading fixed-cohort alternatives...</p>}{error && <p className="context-unavailable" role="alert">{error}</p>}{result && !loading && <AlternativesContent result={result} onAddToCompare={onAddToCompare} />}</>}</section>;
}
