import { useEffect, useRef, useState } from "react";
import { rankCellLinesStream, RankRequestError, type RankStageEvent } from "../api/scoringServiceClient";
import { AiAssistantWidget } from "../components/AiAssistantWidget";
import { BoundaryStatement } from "../components/BoundaryStatement";
import { DiagnosticBanner } from "../components/DiagnosticBanner";
import { InsufficientEvidencePanel } from "../components/InsufficientEvidencePanel";
import { QueryForm, type QuerySubmission } from "../components/QueryForm";
import { RankedResultsList } from "../components/RankedResultsList";
import { ResultCategories } from "../research/ResultCategories";
import { CompareTray } from "../research/CompareTray";
import { ComparisonTable } from "../research/ComparisonTable";
import { compareModels, downloadResearchExport, ResearchRequestError } from "../research/client";
import { RankingProgress } from "../components/RankingProgress";
import { BOUNDARY_STATEMENT_FALLBACK } from "../content/pipelineExplainers";
import type { RankResponse } from "../types/scoring";
import type { ComparisonResponse } from "../research/types";
import "../App.css";
import "../research/research.css";

function downloadJson(data: unknown, filename: string) {
  const blob = new Blob([JSON.stringify(data, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
}

function saveBlob(blob: Blob, filename: string) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  link.click();
  URL.revokeObjectURL(url);
}

function RankPage() {
  const [result, setResult] = useState<RankResponse | null>(null);
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [stages, setStages] = useState<RankStageEvent[]>([]);
  const [selectedModelIds, setSelectedModelIds] = useState<string[]>([]);
  const [comparison, setComparison] = useState<ComparisonResponse | null>(null);
  const [comparisonLoading, setComparisonLoading] = useState(false);
  const [comparisonError, setComparisonError] = useState<string | null>(null);
  const [comparisonNotice, setComparisonNotice] = useState<string | null>(null);
  const comparisonAbortRef = useRef<AbortController | null>(null);
  const activeResearchQueryRef = useRef<string | null>(null);

  const activeResearchQueryId = result?.research_query_id ?? null;
  activeResearchQueryRef.current = activeResearchQueryId;

  useEffect(() => {
    comparisonAbortRef.current?.abort();
    setSelectedModelIds([]);
    setComparison(null);
    setComparisonError(null);
    setComparisonNotice(null);
    setComparisonLoading(false);
  }, [activeResearchQueryId]);

  function toggleComparison(modelId: string) {
    if (!activeResearchQueryId) {
      setComparisonNotice("Comparison needs a trusted research snapshot for this query.");
      return;
    }
    setComparisonError(null);
    setComparisonNotice(null);
    setComparison(null);
    setSelectedModelIds((current) => {
      if (current.includes(modelId)) return current.filter((id) => id !== modelId);
      if (current.length >= 3) {
        setComparisonNotice("Comparison is limited to three models. Remove one before adding another.");
        return current;
      }
      return [...current, modelId];
    });
  }

  async function buildComparison() {
    if (!activeResearchQueryId || selectedModelIds.length < 2 || selectedModelIds.length > 3) return;
    comparisonAbortRef.current?.abort();
    const controller = new AbortController();
    comparisonAbortRef.current = controller;
    const expectedQueryId = activeResearchQueryId;
    const expectedSelection = [...selectedModelIds];
    setComparisonLoading(true);
    setComparisonError(null);
    try {
      const payload = await compareModels(expectedQueryId, expectedSelection, controller.signal);
      if (!controller.signal.aborted && activeResearchQueryRef.current === expectedQueryId) setComparison(payload);
    } catch (requestError) {
      if (!controller.signal.aborted) setComparisonError(requestError instanceof ResearchRequestError ? requestError.message : "Comparison could not be loaded.");
    } finally {
      if (!controller.signal.aborted) setComparisonLoading(false);
    }
  }

  async function downloadResearch(format: "json" | "csv") {
    if (!activeResearchQueryId || selectedModelIds.length < 2 || selectedModelIds.length > 3) return;
    const expectedQueryId = activeResearchQueryId;
    setComparisonError(null);
    try {
      const blob = await downloadResearchExport(expectedQueryId, selectedModelIds, format);
      if (activeResearchQueryRef.current === expectedQueryId) saveBlob(blob, format === "json" ? "cell_line_research_record.json" : "cell_line_comparison.csv");
    } catch (requestError) {
      setComparisonError(requestError instanceof ResearchRequestError ? requestError.message : "Research download failed.");
    }
  }

  const partitionByModel = result ? Object.fromEntries([
    ...result.ranked_cell_lines.map((line) => [line.model_id, "ranked_cell_lines"]),
    ...result.ranked_beyond_top_n.map((line) => [line.model_id, "ranked_beyond_top_n"]),
    ...result.low_confidence_lines.map((line) => [line.model_id, "low_confidence_lines"]),
    ...result.insufficient_evidence_lines.map((line) => [line.model_id, "insufficient_evidence_lines"]),
    ...result.disqualified_lines.map((line) => [line.model_id, "disqualified_lines"]),
  ]) : {};

  async function handleSubmit(query: QuerySubmission) {
    setIsSubmitting(true);
    setError(null);
    setResult(null);
    setStages([]);
    try {
      const response = await rankCellLinesStream(
        {
          // Sent as ensembl_id, not symbol -- the researcher picked one exact gene from the
          // dropdown, and some symbols resolve to more than one Ensembl ID (ALGORITHM_SPEC.md
          // edge case 6); the ID is what they actually selected.
          inclusion_genes: query.inclusionGenes.map((g) => g.ensembl_id),
          exclusion_genes: query.exclusionGenes.map((g) => g.ensembl_id),
          lineage: query.lineage,
          primary_disease: query.primaryDisease,
          exclude_problematic: query.excludeProblematic,
          msi_high: query.msiHigh,
          ploidy_min: query.ploidyMin,
          ploidy_max: query.ploidyMax,
          require_metabolomics: query.requireMetabolomics,
          require_mirna: query.requireMirna,
          top_k: query.topK,
        },
        { onStage: (stage) => setStages((prev) => [...prev, stage]) }
      );
      setResult(response);
    } catch (err) {
      setError(err instanceof RankRequestError ? err.message : "Ranking request failed. Is scoring_service running?");
    } finally {
      setIsSubmitting(false);
    }
  }

  return (
    <div className="app-shell">
      <header>
        <h1>CellLineSelector</h1>
        <p className="subtitle">
          Enter inclusion and exclusion genes to rank candidate cell lines, with a step-by-step
          mathematical and biological explanation for every result.
        </p>
      </header>

      <QueryForm onSubmit={handleSubmit} isSubmitting={isSubmitting} />

      {isSubmitting && <RankingProgress stages={stages} />}

      {error && <p className="error-notice">{error}</p>}

      {result && (
        <section className="results-section">
          <div className="results-header">
            <h2>Results</h2>
            <button type="button" onClick={() => downloadJson(result, "cell_line_selector_native_ranking.json")}>
              Download native ranking JSON
            </button>
          </div>

          <BoundaryStatement text={result.boundary_statement || BOUNDARY_STATEMENT_FALLBACK} />

          <CompareTray
            queryId={result.research_query_id}
            selectedModelIds={selectedModelIds}
            notice={comparisonNotice}
            isLoading={comparisonLoading}
            error={comparisonError}
            onRemove={toggleComparison}
            onCompare={buildComparison}
            onDownload={downloadResearch}
          />

          {comparison && <ComparisonTable comparison={comparison} partitionByModel={partitionByModel} />}

          {result.diagnostic && <DiagnosticBanner diagnostic={result.diagnostic} />}

          {result.ranked_cell_lines.length > 0 ? (
            <RankedResultsList
              lines={result.ranked_cell_lines}
              candidateCounts={result.candidate_counts}
              researchQueryId={result.research_query_id}
              researchStatus={result.research_status}
              researchReason={result.research_reason}
              selectedModelIds={selectedModelIds}
              onToggleCompare={toggleComparison}
              onAddToCompare={toggleComparison}
            />
          ) : (
            <p className="no-recommendations">
              No eligible High- or Moderate-confidence recommendations were returned. Inspect the diagnostic and
              non-recommended partitions below; they are not rank positions.
            </p>
          )}

          <ResultCategories
            rankedBeyondTopN={result.ranked_beyond_top_n}
            lowConfidence={result.low_confidence_lines}
            disqualified={result.disqualified_lines}
            insufficient={<InsufficientEvidencePanel lines={result.insufficient_evidence_lines} researchQueryId={result.research_query_id} researchStatus={result.research_status} researchReason={result.research_reason} selectedModelIds={selectedModelIds} onToggleCompare={toggleComparison} onAddToCompare={toggleComparison} />}
            researchQueryId={result.research_query_id}
            researchStatus={result.research_status}
            researchReason={result.research_reason}
            selectedModelIds={selectedModelIds}
            onToggleCompare={toggleComparison}
            onAddToCompare={toggleComparison}
          />
        </section>
      )}

      <AiAssistantWidget result={result} />
    </div>
  );
}

export default RankPage;
