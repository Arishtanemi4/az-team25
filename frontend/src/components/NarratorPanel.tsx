import { useEffect, useRef, useState } from "react";
import { narrateResult, RagRequestError } from "../api/ragServiceClient";
import type { NarrationResponse } from "../types/rag";
import type { RankResponse } from "../types/scoring";
import { formatCellLineLabel } from "../lib/formatCellLineLabel";

interface NarratorPanelProps {
  evidenceRecord: RankResponse;
}

// The narration response's own CellLineExplanation carries only model_id, no name (types/rag.ts) --
// the name is looked up here from the same evidence record already passed in, rather than adding
// a second field the rag/ service would have to keep in sync (V7-1c).
function buildCellLineNameLookup(evidenceRecord: RankResponse): Map<string, string | null | undefined> {
  const lookup = new Map<string, string | null | undefined>();
  for (const partition of [
    evidenceRecord.ranked_cell_lines,
    evidenceRecord.ranked_beyond_top_n,
    evidenceRecord.low_confidence_lines,
    evidenceRecord.insufficient_evidence_lines,
    evidenceRecord.disqualified_lines,
  ]) {
    for (const line of partition) lookup.set(line.model_id, line.cell_line_name);
  }
  return lookup;
}

// One narration call covers the WHOLE result set, not one per row -- narrator.narrate()'s own
// design rule (rag/narrator.py) is one generation call per result set, so this is button-
// triggered here rather than fired automatically when results render or a row expands.
export function NarratorPanel({ evidenceRecord }: NarratorPanelProps) {
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<NarrationResponse | null>(null);
  const cellLineNames = buildCellLineNameLookup(evidenceRecord);
  const abortControllerRef = useRef<AbortController | null>(null);

  // NarratorPanel only renders while a ranking result exists (AiAssistantWidget.tsx) --
  // RankPage clears that result the moment a new query starts, unmounting this panel mid-request
  // if narration is still in flight. Without this, the fetch kept running after unmount: a
  // wasted, rate-limited LLM call whose eventual response tried to setState on a gone component.
  useEffect(() => {
    return () => {
      abortControllerRef.current?.abort();
    };
  }, []);

  async function handleGenerate() {
    const controller = new AbortController();
    abortControllerRef.current = controller;
    setIsLoading(true);
    setError(null);
    try {
      setResult(await narrateResult(evidenceRecord, 5, controller.signal));
    } catch (err) {
      if (err instanceof DOMException && err.name === "AbortError") return;
      setError(err instanceof RagRequestError ? err.message : "Narration request failed. Is rag_service running?");
    } finally {
      abortControllerRef.current = null;
      setIsLoading(false);
    }
  }

  function handleCancel() {
    abortControllerRef.current?.abort();
  }

  return (
    <div className="narrator-panel">
      <button type="button" onClick={handleGenerate} disabled={isLoading}>
        {isLoading ? "Generating narrative..." : "Generate AI narrative"}
      </button>
      {isLoading && (
        <span className="narrator-loading-note">
          This calls an LLM over the full result set and can take 20-60 seconds (longer if a
          repair pass is needed).{" "}
          <button type="button" onClick={handleCancel}>Cancel</button>
        </span>
      )}
      {error && <p className="error-notice">{error}</p>}
      {result && (
        <div className="narrator-output">
          <p className="ai-disclaimer">
            AI-generated narrative, verified against the evidence record above -- distinct from
            the deterministic per-gene explanations in each result's evidence table.
          </p>
          <p>{result.narration.overview}</p>
          {result.narration.cell_lines.map((line) => (
            <div key={line.model_id} className="narrator-cell-line">
              <strong>{formatCellLineLabel(line.model_id, cellLineNames.get(line.model_id))}</strong>
              <p>{line.summary}</p>
              <ul className="narrator-gene-explanations">
                {line.gene_explanations.map((gene) => (
                  <li key={gene.ensembl_id}>
                    <strong>{gene.symbol}</strong>: {gene.explanation}
                  </li>
                ))}
              </ul>
              <p className="muted">Limitations: {line.limitations}</p>
            </div>
          ))}
          <p className="boundary-statement">{result.narration.boundary_statement}</p>
        </div>
      )}
    </div>
  );
}
