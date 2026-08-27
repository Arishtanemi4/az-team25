import { useRef, useState } from "react";
import { searchGenes } from "../api/geneServiceClient";
import { expandQuery, RagRequestError } from "../api/ragServiceClient";
import type { ExpansionSuggestion } from "../types/rag";
import type { GeneMatch } from "../types/scoring";

interface QueryExpansionPanelProps {
  inclusionGenes: GeneMatch[];
  onAddGene: (gene: GeneMatch) => void;
}

// Suggestions are proposals only -- clicking "Add" is the researcher's own decision, nothing
// here ever changes the query on its own (rag/query_expansion_agent.py's own R2 rule).
export function QueryExpansionPanel({ inclusionGenes, onAddGene }: QueryExpansionPanelProps) {
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [suggestions, setSuggestions] = useState<ExpansionSuggestion[]>([]);
  const [addingGene, setAddingGene] = useState<string | null>(null);
  // /expand is a multi-turn LLM tool-calling agent, routinely 60-90s+ -- not a hang, but nothing
  // distinguishes "still working" from "frozen" without this, so a researcher can cancel out of
  // a slow call instead of wondering whether the page needs a refresh.
  const abortRef = useRef<AbortController | null>(null);

  async function handleSuggest() {
    setIsLoading(true);
    setError(null);
    const controller = new AbortController();
    abortRef.current = controller;
    try {
      const response = await expandQuery(inclusionGenes.map((g) => g.symbol), controller.signal);
      // Suggestions missing source_db/edge_type failed the structural provenance check --
      // shown de-emphasized below, never offered as addable.
      setSuggestions(response.suggestions);
    } catch (err) {
      if (err instanceof DOMException && err.name === "AbortError") {
        // Cancelled by the researcher -- not a failure, say nothing.
      } else {
        setError(err instanceof RagRequestError ? err.message : "Query expansion failed. Is rag_service running?");
      }
    } finally {
      abortRef.current = null;
      setIsLoading(false);
    }
  }

  function handleCancel() {
    abortRef.current?.abort();
  }

  async function handleAdd(suggestion: ExpansionSuggestion) {
    setAddingGene(suggestion.gene);
    try {
      // gene_service's search isn't exact-match-first (e.g. "GRB2" ranks ADGRB2 above GRB2
      // itself, since ADGRB2's symbol contains "GRB2" as a substring) -- take a few results and
      // prefer the one whose symbol matches the suggestion exactly, so "Add" never silently adds
      // a different gene than the one the knowledge graph actually suggested.
      const matches = await searchGenes(suggestion.gene, 5);
      const exact = matches.find((g) => g.symbol.toUpperCase() === suggestion.gene.toUpperCase());
      const match = exact ?? matches[0];
      if (match) {
        onAddGene(match);
      } else {
        setError(`No gene match found for "${suggestion.gene}".`);
      }
    } catch (err) {
      // A failed resolution just leaves the suggestion in place -- the researcher can retry --
      // but they need to be told it failed rather than seeing "Add" silently reset with no
      // explanation.
      setError(err instanceof Error ? err.message : `Failed to resolve "${suggestion.gene}". Is gene_service running?`);
    } finally {
      setAddingGene(null);
    }
  }

  return (
    <div className="expansion-panel">
      <button type="button" onClick={handleSuggest} disabled={inclusionGenes.length === 0 || isLoading}>
        {isLoading ? "Consulting knowledge graph..." : "Suggest related genes"}
      </button>
      {isLoading && (
        <span className="expansion-loading-note">
          This runs a multi-step reasoning agent and can take 1-2 minutes.{" "}
          <button type="button" onClick={handleCancel}>Cancel</button>
        </span>
      )}
      {error && <p className="error-notice">{error}</p>}
      {suggestions.length > 0 && (
        <ul className="chip-row">
          {suggestions.map((suggestion) => {
            const hasProvenance = Boolean(suggestion.source_db && suggestion.edge_type);
            return (
              <li key={suggestion.gene} className="expansion-suggestion">
                <span className="chip">
                  {suggestion.gene}
                  {hasProvenance && (
                    <span className="provenance-badge">
                      {suggestion.source_db}/{suggestion.edge_type}
                    </span>
                  )}
                  {hasProvenance && (
                    <button
                      type="button"
                      onClick={() => handleAdd(suggestion)}
                      disabled={addingGene === suggestion.gene}
                    >
                      {addingGene === suggestion.gene ? "Adding..." : "Add"}
                    </button>
                  )}
                </span>
                <span className="muted">{suggestion.reason}</span>
              </li>
            );
          })}
        </ul>
      )}
    </div>
  );
}
