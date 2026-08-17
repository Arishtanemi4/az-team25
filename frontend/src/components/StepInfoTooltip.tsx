import { useState } from "react";

interface StepInfoTooltipProps {
  title: string;
  formula: string;
  biology: string;
}

// The (i) icon on each pipeline step. Click-to-toggle (not just hover) so it also works on
// touch devices; the popover content is static text from content/pipelineExplainers.ts, not
// fetched per query, since the formula/rationale never changes between queries.
export function StepInfoTooltip({ title, formula, biology }: StepInfoTooltipProps) {
  const [isOpen, setIsOpen] = useState(false);

  return (
    <span className="info-tooltip-wrapper">
      <button
        type="button"
        className="info-icon"
        aria-label={`How ${title} is computed`}
        aria-expanded={isOpen}
        onClick={() => setIsOpen((v) => !v)}
        onMouseEnter={() => setIsOpen(true)}
        onMouseLeave={() => setIsOpen(false)}
      >
        i
      </button>
      {isOpen && (
        <div className="info-popover" role="tooltip">
          <p className="info-popover-formula">{formula}</p>
          <p className="info-popover-biology">{biology}</p>
        </div>
      )}
    </span>
  );
}
