import { useState } from "react";
import { askMethodologyQuestion, RagRequestError } from "../api/ragServiceClient";
import type { MethodologyAnswer } from "../types/rag";

// Standalone -- not tied to a specific ranking result, since a researcher may want to ask about
// the methodology (e.g. "why is protein weighted 0.27?") before or independently of any query.
export function MethodologyQABox() {
  const [question, setQuestion] = useState("");
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [answer, setAnswer] = useState<MethodologyAnswer | null>(null);

  async function handleSubmit(e: React.FormEvent) {
    e.preventDefault();
    if (!question.trim() || isLoading) return;
    setIsLoading(true);
    setError(null);
    try {
      setAnswer(await askMethodologyQuestion(question));
    } catch (err) {
      setError(err instanceof RagRequestError ? err.message : "Methodology question failed. Is rag_service running?");
    } finally {
      setIsLoading(false);
    }
  }

  return (
    <section className="methodology-qa">
      <h2>Ask about the methodology</h2>
      <form onSubmit={handleSubmit}>
        <input
          type="text"
          value={question}
          placeholder="e.g. why is protein weighted 0.27?"
          onChange={(e) => setQuestion(e.target.value)}
        />
        <button type="submit" disabled={!question.trim() || isLoading}>
          {isLoading ? "Asking..." : "Ask"}
        </button>
      </form>
      {error && <p className="error-notice">{error}</p>}
      {answer && (
        <div className="methodology-answer">
          {answer.abstained && (
            <p className="muted">The project documentation does not cover this question.</p>
          )}
          <p>{answer.answer}</p>
        </div>
      )}
    </section>
  );
}
