import { PIPELINE_STEPS } from "../content/pipelineExplainers";
import type { RankedCellLine } from "../types/scoring";
import { PerGeneEvidenceTable } from "./PerGeneEvidenceTable";
import { StepInfoTooltip } from "./StepInfoTooltip";
import { BiologicalContextPanel } from "../research/BiologicalContextPanel";
import { AlternativesPanel } from "../research/AlternativesPanel";
import { formatScore } from "../lib/formatScore";

interface CellLineDetailProps {
  line: RankedCellLine;
  researchQueryId?: string | null;
  researchStatus?: string;
  researchReason?: string;
  onAddToCompare?: (modelId: string) => void;
}

function formatD(d: number | null): string {
  return formatScore(d);
}

function agreementLabel(value: boolean | null): string {
  if (value === true) return "consistent";
  if (value === false) return "not consistent";
  return "not assessed";
}

// The confirmed 4-step pipeline stepper: per-gene desirability -> within-gene layer
// combination -> across-gene veto combination -> correlation discount. Each step shows this
// result's actual numbers; the (i) icon reveals the general formula/biology, which is the
// same for every cell line and every query.
export function CellLineDetail({ line, researchQueryId, researchStatus, researchReason, onAddToCompare }: CellLineDetailProps) {
  const inclusionGenes = line.per_gene.filter((g) => g.role === "inclusion");
  const exclusionGenes = line.per_gene.filter((g) => g.role === "exclusion");

  return (
    <div className="cell-line-detail">
      <ol className="pipeline-steps">
        <li className="pipeline-step">
          <div className="pipeline-step-header">
            <span className="step-number">{PIPELINE_STEPS[0].step}</span>
            <span className="step-title">{PIPELINE_STEPS[0].title}</span>
            <StepInfoTooltip {...PIPELINE_STEPS[0]} />
          </div>
          <div className="pipeline-step-body">
            {line.per_gene.map((gene) => (
              <div key={gene.ensembl_id} className="gene-layer-summary">
                <strong>
                  {gene.symbol} ({gene.role})
                </strong>
                {": "}
                {Object.entries(gene.layers)
                  .filter(([, detail]) => detail.d !== null)
                  .map(([layer, detail]) => `${layer} d=${formatScore(detail.d)}`)
                  .join(", ") || "no evidence in any layer"}
              </div>
            ))}
          </div>
        </li>

        <li className="pipeline-step">
          <div className="pipeline-step-header">
            <span className="step-number">{PIPELINE_STEPS[1].step}</span>
            <span className="step-title">{PIPELINE_STEPS[1].title}</span>
            <StepInfoTooltip {...PIPELINE_STEPS[1]} />
          </div>
          <div className="pipeline-step-body">
            {line.per_gene.map((gene) => (
              <div key={gene.ensembl_id} className="gene-d">
                {gene.symbol}: d_gene = {formatD(gene.d_gene)}
              </div>
            ))}
          </div>
        </li>

        <li className="pipeline-step">
          <div className="pipeline-step-header">
            <span className="step-number">{PIPELINE_STEPS[2].step}</span>
            <span className="step-title">{PIPELINE_STEPS[2].title}</span>
            <StepInfoTooltip {...PIPELINE_STEPS[2]} />
          </div>
          <div className="pipeline-step-body">
            <div className="headline-score">D = {formatD(line.D)}</div>
            {line.veto ? (
              <div className="veto-notice">
                Vetoed by <strong>{line.veto.symbol}</strong> (
                {line.veto.reason === "exclusion_expressed"
                  ? "exclusion gene is expressed"
                  : "inclusion gene is below its detection floor"}
                )
              </div>
            ) : (
              <div className="veto-clear">No veto triggered.</div>
            )}
          </div>
        </li>

        <li className="pipeline-step">
          <div className="pipeline-step-header">
            <span className="step-number">{PIPELINE_STEPS[3].step}</span>
            <span className="step-title">{PIPELINE_STEPS[3].title}</span>
            <StepInfoTooltip {...PIPELINE_STEPS[3]} />
          </div>
          <div className="pipeline-step-body">
            <div>
              m_eff = {formatScore(line.m_eff)} of {line.m} scored genes (some overlap discounted)
            </div>
            {[...inclusionGenes, ...exclusionGenes].map((gene) => {
              const rho = line.rho_bar[gene.ensembl_id];
              return (
                <div key={gene.ensembl_id} className="rho-bar-entry">
                  {gene.symbol}: rho-bar = {rho === null || rho === undefined ? "n/a" : formatScore(rho)}
                </div>
              );
            })}
          </div>
        </li>
      </ol>

      <section className="native-evidence-summary" aria-label="Native evidence limitations">
        <h3>Native evidence and confidence boundary</h3>
        <dl className="evidence-status-list">
          <div>
            <dt>HPA corroboration</dt>
            <dd>{agreementLabel(line.hpa_agreement)}</dd>
          </div>
          <div>
            <dt>GEO corroboration</dt>
            <dd>{agreementLabel(line.geo_agreement)}</dd>
          </div>
          <div>
            <dt>Missing evidence</dt>
            <dd>
              {line.missing_evidence.length === 0
                ? "No native missing-evidence entries reported."
                : `${line.missing_evidence.length} entry/entries; missing or unassayed evidence is not interpreted as zero.`}
            </dd>
          </div>
        </dl>
        {line.warnings.length > 0 && (
          <ul className="native-warning-list">
            {line.warnings.map((warning, index) => <li key={`${warning.type}-${index}`}>{warning.message}</li>)}
          </ul>
        )}
        {line.missing_evidence.length > 0 && (
          <ul className="missing-evidence-list">
            {line.missing_evidence.map((item, index) => (
              <li key={`${item.ensembl_id}-${item.layer ?? "all"}-${index}`}>
                {item.symbol ?? item.ensembl_id}: {item.severity}{item.layer ? ` (${item.layer})` : ""}
              </li>
            ))}
          </ul>
        )}
        <p className="confidence-boundary-note">
          Confidence describes consistency of available evidence; it is not experimental validation. An unassayed layer remains unknown.
        </p>
      </section>

      <PerGeneEvidenceTable perGene={line.per_gene} />
      <BiologicalContextPanel
        queryId={researchQueryId}
        modelId={line.model_id}
        researchStatus={researchStatus}
        researchReason={researchReason}
      />
      <AlternativesPanel queryId={researchQueryId} modelId={line.model_id} onAddToCompare={onAddToCompare ?? (() => undefined)} />
    </div>
  );
}
