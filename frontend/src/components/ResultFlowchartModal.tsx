import { useEffect, useRef } from "react";
import type { RankedCellLine } from "../types/scoring";
import { formatCellLineLabel } from "../lib/formatCellLineLabel";
import { ResultFlowchartDiagram } from "./ResultFlowchartDiagram";
import "./resultFlowchart.css";

interface ResultFlowchartModalProps {
  line: RankedCellLine;
  partitionBucket: string;
  onClose: () => void;
}

// V7-8: the result path flowchart, opened from a ranked row. Read-only trace view of one
// already-fetched result -- extends CellLineDetail.tsx's 4-step stepper rather than replacing
// it (that stepper stays exactly as it is). No portal infrastructure exists anywhere in this
// codebase, so this is a plain fixed-position overlay, consistent with AiAssistantWidget.tsx.
export function ResultFlowchartModal({ line, partitionBucket, onClose }: ResultFlowchartModalProps) {
  const dialogRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    function handleKeyDown(event: KeyboardEvent) {
      if (event.key === "Escape") onClose();
    }
    document.addEventListener("keydown", handleKeyDown);
    return () => document.removeEventListener("keydown", handleKeyDown);
  }, [onClose]);

  useEffect(() => {
    dialogRef.current?.focus();
  }, []);

  return (
    <div className="flowchart-modal-backdrop" onClick={onClose}>
      <div
        className="flowchart-modal"
        role="dialog"
        aria-modal="true"
        aria-label={`Result path for ${formatCellLineLabel(line.model_id, line.cell_line_name)}`}
        tabIndex={-1}
        ref={dialogRef}
        onClick={(event) => event.stopPropagation()}
      >
        <div className="flowchart-modal-header">
          <h2>{formatCellLineLabel(line.model_id, line.cell_line_name)} -- result path</h2>
          <button type="button" className="flowchart-modal-close" aria-label="Close" onClick={onClose}>
            ×
          </button>
        </div>
        <p className="flowchart-modal-note">
          Every node below is read from this result's own record, or explicitly marked "not traced"
          when the record does not carry it. Nothing here is re-scored or inferred.
        </p>
        <ResultFlowchartDiagram line={line} partitionBucket={partitionBucket} />
      </div>
    </div>
  );
}
