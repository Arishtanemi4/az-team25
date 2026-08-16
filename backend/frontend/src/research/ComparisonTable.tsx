import type { ComparisonResponse } from "./types";
import { formatScore } from "../lib/formatScore";

interface ComparisonTableProps {
  comparison: ComparisonResponse;
  partitionByModel: Record<string, string>;
}

function displayValue(value: unknown): string {
  if (value === null || value === undefined) return "not reported";
  if (typeof value === "number") return formatScore(value);
  if (typeof value === "object") return JSON.stringify(value);
  return String(value);
}

function partitionLabel(partition: string | undefined): string {
  if (partition === "ranked_cell_lines") return "recommended";
  if (partition === "ranked_beyond_top_n") return "eligible tail";
  if (partition === "low_confidence_lines") return "Low confidence";
  if (partition === "insufficient_evidence_lines") return "insufficient evidence";
  if (partition === "disqualified_lines") return "disqualified";
  return "native partition not reported";
}

export function ComparisonTable({ comparison, partitionByModel }: ComparisonTableProps) {
  return (
    <section className="comparison-table-panel">
      <h3>Side-by-side native evidence</h3>
      <p className="context-boundary">{comparison.boundary}</p>
      <div className="comparison-table-scroll">
        <table className="comparison-table">
          <thead>
            <tr>
              <th>Evidence field</th>
              <th>Source / unit</th>
              {comparison.selected_model_ids.map((modelId) => <th key={modelId}><code>{modelId}</code><span className="partition-badge">{partitionLabel(partitionByModel[modelId])}</span></th>)}
            </tr>
          </thead>
          <tbody>
            {comparison.rows.map((row) => (
              <tr key={row.row_id}>
                <td><strong>{row.label}</strong>{row.gene_id && <span className="comparison-subtext">{row.gene_id}{row.gene_role ? ` (${row.gene_role})` : ""}</span>}</td>
                <td>{row.layer ?? row.section}{row.unit ? `; ${row.unit}` : ""}<span className="comparison-subtext">{row.provenance}</span></td>
                {comparison.selected_model_ids.map((modelId) => <td key={modelId}>{displayValue(row.values[modelId])}</td>)}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      {comparison.limitations.length > 0 && <ul className="context-caveats">{comparison.limitations.map((limitation, index) => <li key={index}>{limitation}</li>)}</ul>}
    </section>
  );
}
