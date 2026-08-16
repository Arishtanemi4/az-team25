import { useState } from "react";
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

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    if (!question.trim() || isLoading) return;
    setIsLoading(true);
    setError(null);
    try {
      setResult(await searchLiterature(question));
    } catch (err) {
      setError(err instanceof RagRequestError ? err.message : "Literature search failed. Is rag_service running?");
    } finally {
      setIsLoading(false);
    }
  }

  return (
    <section className="literature-search">
      <h2>Search the literature</h2>
      <form onSubmit={handleSubmit}>
        <input
          type="text"
          value={question}
          placeholder="e.g. does EGFR expression predict drug response?"
          onChange={(e) => setQuestion(e.target.value)}
        />
        <button type="submit" disabled={!question.trim() || isLoading}>
          {isLoading ? "Searching PubMed..." : "Search"}
        </button>
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
