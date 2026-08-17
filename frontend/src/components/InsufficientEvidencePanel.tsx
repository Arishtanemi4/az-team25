import { useRef, useState } from "react";
import type { RankedCellLine } from "../types/scoring";
import { CellLineDetail } from "./CellLineDetail";
import { ResultFlowchartModal } from "./ResultFlowchartModal";
import { formatCellLineLabel } from "../lib/formatCellLineLabel";

interface InsufficientEvidencePanelProps {
  lines: RankedCellLine[];
  researchQueryId?: string | null;
  researchStatus?: string;
  researchReason?: string;
  selectedModelIds?: string[];
  onToggleCompare?: (modelId: string) => void;
  onAddToCompare?: (modelId: string) => void;
}

// Lines with insufficient evidence are never ranked last -- they're listed separately, since
// unmeasured is not the same as unsuitable (ALGORITHM_SPEC.md edge case 4).
export function InsufficientEvidencePanel({ lines, researchQueryId, researchStatus, researchReason, selectedModelIds = [], onToggleCompare, onAddToCompare }: InsufficientEvidencePanelProps) {
  const [isExpanded, setIsExpanded] = useState(false);
  const [expandedId, setExpandedId] = useState<string | null>(null);
  const [flowchartLine, setFlowchartLine] = useState<RankedCellLine | null>(null);
  const flowchartTriggerRef = useRef<HTMLButtonElement | null>(null);

  if (lines.length === 0) return null;

  return (
    <div className="insufficient-evidence-panel">
      <button type="button" className="text-button" onClick={() => setIsExpanded((v) => !v)}>
        {isExpanded ? "Hide" : "Show"} {lines.length} cell line(s) with insufficient evidence (not ranked)
      </button>
      {isExpanded && (
        <ul>
          {lines.map((line) => (
            <li key={line.model_id} className="insufficient-line">
              <button type="button" className="text-button" onClick={() => setExpandedId(expandedId === line.model_id ? null : line.model_id)} aria-expanded={expandedId === line.model_id}>
                {formatCellLineLabel(line.model_id, line.cell_line_name)} — {line.missing_evidence.length} missing-evidence entry/entries; unmeasured is not unsuitable.
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
          partitionBucket="insufficient_evidence"
          onClose={() => {
            setFlowchartLine(null);
            flowchartTriggerRef.current?.focus();
          }}
        />
      )}
    </div>
  );
}
