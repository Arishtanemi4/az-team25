import { useRef, useState } from "react";
import { RagRequestError, searchLiterature } from "../api/ragServiceClient";
import type { LiteratureResponse } from "../types/rag";

// Kept separate from MethodologyQABox -- different source (live PubMed, not project docs),
// different response shape (quote+PMID findings, not prose), and merging the two behind a mode
// toggle would be exactly the kind of unrequested configurability root _.md SS7 warns against.
export function LiteratureSearchBox() {
  const [question, setQuestion] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<LiteratureResponse | null>(null);
  const [showHint, setShowHint] = useState(false);
  // find_evidence() is a multi-turn agent hitting live PubMed/NCBI, routinely slow -- lets the
  // researcher cancel out instead of wondering whether the search has hung (same convention as
  // QueryExpansionPanel/NarratorPanel). No unmount-abort effect needed: this box stays mounted
  // for the AI assistant widget's whole lifetime (AiAssistantWidget.tsx), only CSS-hidden when
  // its tab isn't active.
  const abortRef = useRef<AbortController | null>(null);

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    if (!question.trim() || isLoading) return;
    const controller = new AbortController();
    abortRef.current = controller;
    setIsLoading(true);
    setError(null);
    try {
      setResult(await searchLiterature(question, controller.signal));
    } catch (err) {
      if (err instanceof DOMException && err.name === "AbortError") return;
      setError(err instanceof RagRequestError ? err.message : "Literature search failed. Is rag_service running?");
    } finally {
      abortRef.current = null;
      setIsLoading(false);
    }
  }

  function handleCancel() {
    abortRef.current?.abort();
  }

  return (
    <section className="literature-search">
      <h2>Search the literature</h2>
      <form onSubmit={handleSubmit}>
        <input
          type="text"
          value={question}
          placeholder="e.g. Is BRAF V600E associated with vemurafenib sensitivity in melanoma?"
          onChange={(e) => setQuestion(e.target.value)}
        />
        <button type="submit" disabled={!question.trim() || isLoading}>
          {isLoading ? "Searching PubMed..." : "Search"}
        </button>
        <span className="info-tooltip-wrapper">
          <button
            type="button"
            className="info-icon"
            aria-label="Example questions"
            aria-expanded={showHint}
            onClick={() => setShowHint((v) => !v)}
            onMouseEnter={() => setShowHint(true)}
            onMouseLeave={() => setShowHint(false)}
          >
            i
          </button>
          {showHint && (
            <div className="info-popover" role="tooltip">
              <p>
                Ask about a specific gene-drug or gene-outcome relationship it can verify against
                PubMed, e.g.:
              </p>
              <ul>
                <li>Is BRAF V600E associated with vemurafenib sensitivity in melanoma?</li>
                <li>Does HER2 amplification predict trastuzumab response in breast cancer?</li>
                <li>Is EGFR exon 19 deletion associated with TKI sensitivity in lung cancer?</li>
              </ul>
              <p className="muted">Broad or off-topic questions are less likely to get a useful answer.</p>
            </div>
          )}
        </span>
        {isLoading && (
          <span className="literature-loading-note">
            This searches live PubMed through a multi-step LLM agent and can take 1-2 minutes or
            more.{" "}
            <button type="button" onClick={handleCancel}>Cancel</button>
          </span>
        )}
      </form>
      {error && <p className="error-notice">{error}</p>}
      {result && (
        <div className="literature-results">
          <p className="muted">
            {result.n_quote_confirmed} of {result.n_reported_by_model} claimed findings were
            quote-verified against a fetched abstract.
          </p>
          {result.findings.length === 0 ? (
            <p className="muted">No quote-verified findings.</p>
          ) : (
            <ul>
              {result.findings.map((finding) => (
                <li key={finding.pmid + finding.quote} className="literature-finding">
                  <blockquote>{finding.quote}</blockquote>
                  <p>
                    <a
                      href={`https://pubmed.ncbi.nlm.nih.gov/${finding.pmid}/`}
                      target="_blank"
                      rel="noreferrer"
                    >
                      PMID {finding.pmid}
                    </a>
                    {" -- "}
                    {finding.relevance}
                  </p>
                </li>
              ))}
            </ul>
          )}
        </div>
      )}
    </section>
  );
}
