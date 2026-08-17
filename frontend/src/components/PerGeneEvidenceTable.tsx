import { useState } from "react";
import type { PerGeneResult } from "../types/scoring";
import { formatScore } from "../lib/formatScore";
import { stateLabel } from "../lib/layerStateLabel";

interface PerGeneEvidenceTableProps {
  perGene: PerGeneResult[];
}

const LAYER_ORDER: Array<keyof PerGeneResult["layers"]> = [
  "rna",
  "protein",
  "fusion",
  "copy_number",
  "mutation",
  "dependency",
];

// Collapsed by default -- this is the "don't overburden the user" detail level: the raw
// value/unit/desirability/state for every layer of every gene, plus the auto-generated
// narrative sentence (scoring/explain.py::gene_narrative), shown only once the researcher
// asks for it.
export function PerGeneEvidenceTable({ perGene }: PerGeneEvidenceTableProps) {
  const [isExpanded, setIsExpanded] = useState(false);

  return (
    <div className="evidence-table-wrapper">
      <button type="button" className="text-button" onClick={() => setIsExpanded((v) => !v)}>
        {isExpanded ? "Hide per-gene evidence table" : "Show per-gene evidence table"}
      </button>
      {isExpanded && (
        <div className="evidence-table-scroll">
          <table className="evidence-table">
            <thead>
              <tr>
                <th>Gene</th>
                <th>Role</th>
                <th>d_gene</th>
                {LAYER_ORDER.map((layer) => (
                  <th key={layer}>{layer}</th>
                ))}
              </tr>
            </thead>
            <tbody>
              {perGene.map((gene) => (
                <tr key={gene.ensembl_id}>
                  <td>{gene.symbol}</td>
                  <td>{gene.role}</td>
                  <td>{formatScore(gene.d_gene)}</td>
                  {LAYER_ORDER.map((layer) => {
                    const detail = gene.layers[layer];
                    return (
                      <td key={layer} className={`layer-state-${detail.state === "not_assayed" ? "not_assayed" : "known"}`}>
                        {detail.d === null ? "—" : formatScore(detail.d)}
                        <span className="layer-value-note">
                          {stateLabel(detail.state)}{detail.value === null ? "" : `; ${detail.value} ${detail.unit}`}
                        </span>
                      </td>
                    );
                  })}
                </tr>
              ))}
            </tbody>
          </table>
          <ul className="gene-narratives">
            {perGene.map((gene) => (
              <li key={gene.ensembl_id}>{gene.narrative ?? "No native narrative was generated for this gene."}</li>
            ))}
          </ul>
        </div>
      )}
    </div>
  );
}
