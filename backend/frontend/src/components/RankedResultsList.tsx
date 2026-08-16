import { useRef, useState } from "react";
import type { CandidateCounts, ConfidenceTier, RankedCellLine } from "../types/scoring";
import { CellLineDetail } from "./CellLineDetail";
import { ResultFlowchartModal } from "./ResultFlowchartModal";
import { formatScore } from "../lib/formatScore";
import { formatCellLineLabel } from "../lib/formatCellLineLabel";

interface RankedResultsListProps {
  lines: RankedCellLine[];
  candidateCounts: CandidateCounts;
  researchQueryId?: string | null;
  researchStatus?: string;
  researchReason?: string;
  selectedModelIds?: string[];
  onToggleCompare?: (modelId: string) => void;
  onAddToCompare?: (modelId: string) => void;
}

const TIER_CLASS: Record<string, string> = {
  High: "tier-high",
  Moderate: "tier-moderate",
  Low: "tier-low",
  Insufficient: "tier-insufficient",
};

function confidenceLabel(tier: ConfidenceTier): string {
  return TIER_CLASS[tier] ? `${tier} confidence` : `Unknown confidence state: ${tier}`;
}

function confidenceClass(tier: ConfidenceTier): string {
  return TIER_CLASS[tier] ?? "tier-unknown";
}

export function RankedResultsList({ lines, candidateCounts, researchQueryId, researchStatus, researchReason, selectedModelIds = [], onToggleCompare, onAddToCompare }: RankedResultsListProps) {
  const [expandedId, setExpandedId] = useState<string | null>(null);
  const [flowchartLine, setFlowchartLine] = useState<RankedCellLine | null>(null);
  const flowchartTriggerRef = useRef<HTMLButtonElement | null>(null);

  return (
    <div className="ranked-results">
      <p className="results-summary">
        Displayed recommendations: <strong>{candidateCounts.displayed}</strong>; eligible High- or Moderate-confidence
        candidates: <strong>{candidateCounts.eligible}</strong>; evaluated candidates: <strong>{candidateCounts.evaluated}</strong>.
      </p>
      <p className="results-order-note">Exact native-score ties retain stable ModelID order; coverage does not break ties.</p>
      <p className="score-precision-note">Scores are shown to four decimal places of arithmetic precision -- this reflects what the calculation computed, not how accurately it was measured.</p>
      <ol className="ranked-list">
        {lines.map((line, index) => (
          <li key={line.model_id} className="ranked-row">
            <button
              type="button"
              className="ranked-row-header"
              onClick={() => setExpandedId(expandedId === line.model_id ? null : line.model_id)}
              aria-expanded={expandedId === line.model_id}
            >
              <span className="rank-index">#{index + 1}</span>
              <span className="model-id">{formatCellLineLabel(line.model_id, line.cell_line_name)}</span>
              <span className="d-score">D = {formatScore(line.D)}</span>
              <span className={`confidence-badge ${confidenceClass(line.confidence_tier)}`}>
                {confidenceLabel(line.confidence_tier)}
              </span>
              {line.is_problematic && <span className="problematic-badge">flagged: identity/contamination</span>}
            </button>
            <button type="button" className="compare-toggle" disabled={!researchQueryId} onClick={() => onToggleCompare?.(line.model_id)} aria-pressed={selectedModelIds.includes(line.model_id)}>
              {selectedModelIds.includes(line.model_id) ? "Remove from comparison" : "Compare"}
            </button>
            <button
              type="button"
              className="flowchart-trigger"
              onClick={(event) => {
                flowchartTriggerRef.current = event.currentTarget;
                setFlowchartLine(line);
              }}
            >
              Show result path
            </button>
            {expandedId === line.model_id && <CellLineDetail line={line} researchQueryId={researchQueryId} researchStatus={researchStatus} researchReason={researchReason} onAddToCompare={onAddToCompare} />}
          </li>
        ))}
      </ol>
      {flowchartLine && (
        <ResultFlowchartModal
          line={flowchartLine}
          partitionBucket="ranked"
          onClose={() => {
            setFlowchartLine(null);
            flowchartTriggerRef.current?.focus();
          }}
        />
      )}
    </div>
  );
}
