import { useRef, useState, type ReactNode } from "react";
import type { ConfidenceTier, RankedCellLine } from "../types/scoring";
import { CellLineDetail } from "../components/CellLineDetail";
import { ResultFlowchartModal } from "../components/ResultFlowchartModal";
import { formatScore } from "../lib/formatScore";
import { formatCellLineLabel } from "../lib/formatCellLineLabel";

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

interface ResultCategoriesProps {
  rankedBeyondTopN: RankedCellLine[];
  lowConfidence: RankedCellLine[];
  disqualified: RankedCellLine[];
  insufficient: ReactNode;
  researchQueryId?: string | null;
  researchStatus?: string;
  researchReason?: string;
  selectedModelIds?: string[];
  onToggleCompare?: (modelId: string) => void;
  onAddToCompare?: (modelId: string) => void;
}

interface CategoryProps {
  title: string;
  explanation: string;
  lines: RankedCellLine[];
  bucket: string;
  researchQueryId?: string | null;
  researchStatus?: string;
  researchReason?: string;
  selectedModelIds: string[];
  onToggleCompare?: (modelId: string) => void;
  onAddToCompare?: (modelId: string) => void;
}

function InspectableCategory({ title, explanation, lines, bucket, researchQueryId, researchStatus, researchReason, selectedModelIds, onToggleCompare, onAddToCompare }: CategoryProps) {
  const [expanded, setExpanded] = useState(false);
  const [expandedId, setExpandedId] = useState<string | null>(null);
  const [flowchartLine, setFlowchartLine] = useState<RankedCellLine | null>(null);
  const flowchartTriggerRef = useRef<HTMLButtonElement | null>(null);

  if (lines.length === 0) return null;
  return (
    <section className="result-category">
      <button type="button" className="category-toggle" onClick={() => setExpanded((value) => !value)} aria-expanded={expanded}>
        {expanded ? "Hide" : "Show"} {lines.length} {title}
      </button>
      <p className="category-explanation">{explanation}</p>
      {expanded && (
        <ul className="unranked-list">
          {lines.map((line) => (
            <li key={line.model_id} className="unranked-row">
              <button type="button" className="unranked-row-header" onClick={() => setExpandedId(expandedId === line.model_id ? null : line.model_id)} aria-expanded={expandedId === line.model_id}>
                <span className="model-id">{formatCellLineLabel(line.model_id, line.cell_line_name)}</span>
                <span className="d-score">D = {formatScore(line.D)}</span>
                <span className={`confidence-badge ${confidenceClass(line.confidence_tier)}`}>{confidenceLabel(line.confidence_tier)}</span>
                {line.veto && <span className="veto-badge">native veto</span>}
              </button>
              <button type="button" className="compare-toggle" disabled={!researchQueryId} onClick={() => onToggleCompare?.(line.model_id)} aria-pressed={selectedModelIds.includes(line.model_id)}>{selectedModelIds.includes(line.model_id) ? "Remove from comparison" : "Compare"}</button>
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
        </ul>
      )}
      {flowchartLine && (
        <ResultFlowchartModal
          line={flowchartLine}
          partitionBucket={bucket}
          onClose={() => {
            setFlowchartLine(null);
            flowchartTriggerRef.current?.focus();
          }}
        />
      )}
    </section>
  );
}

// These are partitions of the complete native output, deliberately shown outside the ordered
// recommendation list. A reader can inspect their evidence without mistaking them for ranks.
export function ResultCategories({ rankedBeyondTopN, lowConfidence, disqualified, insufficient, researchQueryId, researchStatus, researchReason, selectedModelIds = [], onToggleCompare, onAddToCompare }: ResultCategoriesProps) {
  return (
    <div className="result-categories">
      <InspectableCategory
        title="eligible candidates beyond the displayed top N"
        explanation="These candidates passed the native eligibility rule but were not displayed because of your top-N limit. They are not numbered here."
        lines={rankedBeyondTopN}
        bucket="ranked_beyond_top_n"
        researchQueryId={researchQueryId}
        researchStatus={researchStatus}
        researchReason={researchReason}
        selectedModelIds={selectedModelIds}
        onToggleCompare={onToggleCompare}
        onAddToCompare={onAddToCompare}
      />
      <InspectableCategory
        title="Low-confidence candidates"
        explanation="These have native evidence but do not meet the High/Moderate confidence recommendation threshold. They are not recommendations."
        lines={lowConfidence}
        bucket="low_confidence"
        researchQueryId={researchQueryId}
        researchStatus={researchStatus}
        researchReason={researchReason}
        selectedModelIds={selectedModelIds}
        onToggleCompare={onToggleCompare}
        onAddToCompare={onAddToCompare}
      />
      <InspectableCategory
        title="disqualified candidates"
        explanation="These failed a native veto or floor rule. Their evidence is retained for inspection, not offered as an alternative recommendation."
        lines={disqualified}
        bucket="disqualified"
        researchQueryId={researchQueryId}
        researchStatus={researchStatus}
        researchReason={researchReason}
        selectedModelIds={selectedModelIds}
        onToggleCompare={onToggleCompare}
        onAddToCompare={onAddToCompare}
      />
      {insufficient}
    </div>
  );
}
