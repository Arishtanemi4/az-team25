import type {
  ExpansionResponse,
  LiteratureResponse,
  MethodologyAnswer,
  NarrationResponse,
} from "../types/rag";
import type { RankResponse } from "../types/scoring";
import type { GraphNeighborhoodResponse } from "../types/graph";

const BASE_URL = import.meta.env.VITE_RAG_SERVICE_URL ?? "http://localhost:8000";

export class RagRequestError extends Error {}

async function postJson<T>(path: string, body: unknown, signal?: AbortSignal): Promise<T> {
  const response = await fetch(`${BASE_URL}${path}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
    signal,
  });
  if (!response.ok) {
    const problem = await response.json().catch(() => null);
    // FastAPI's error detail is either a string (our own HTTPException) or a list of pydantic
    // validation errors -- handle both so the user sees a readable message either way (same
    // convention as scoringServiceClient.ts).
    const detail = problem?.detail;
    const message = Array.isArray(detail)
      ? detail.map((d: { msg: string }) => d.msg).join("; ")
      : (detail ?? `Request failed (${response.status})`);
    throw new RagRequestError(message);
  }
  return response.json();
}

export function narrateResult(
  evidenceRecord: RankResponse,
  topKContext = 5,
  signal?: AbortSignal,
): Promise<NarrationResponse> {
  return postJson("/narrate", { evidence_record: evidenceRecord, top_k_context: topKContext }, signal);
}

export function askMethodologyQuestion(question: string): Promise<MethodologyAnswer> {
  return postJson("/methodology/ask", { question });
}

// /expand runs a multi-turn tool-calling LLM agent (rag/query_expansion_agent.py) against a
// large model -- routinely 60-90s+ for a handful of turns, not a hang. `signal` lets the panel
// offer a real Cancel button instead of the researcher wondering whether the UI has frozen.
export function expandQuery(inclusionGenes: string[], signal?: AbortSignal): Promise<ExpansionResponse> {
  return postJson("/expand", { inclusion_genes: inclusionGenes }, signal);
}

export function searchLiterature(question: string): Promise<LiteratureResponse> {
  return postJson("/literature/search", { question });
}

// The one GET call in this client -- /graph/neighborhood (PRODUCT_SURFACE.md SS3.5) reads a
// bounded k-hop neighbourhood, it doesn't post a request body. readErrorDetail's FastAPI-detail
// convention is duplicated inline here rather than shared with postJson, since the two throw
// different error shapes (postJson throws RagRequestError; this stays a plain Error to match
// dataServiceClient.ts's own GET-call convention, which this call is closer in shape to).
export async function getGraphNeighborhood(
  gene: string,
  hops: number,
  maxNeighbors: number,
): Promise<GraphNeighborhoodResponse> {
  const url = `${BASE_URL}/graph/neighborhood?gene=${encodeURIComponent(gene)}&hops=${hops}&max_neighbors=${maxNeighbors}`;
  const response = await fetch(url);
  if (!response.ok) {
    const problem = await response.json().catch(() => null);
    const detail = problem?.detail;
    const message = Array.isArray(detail)
      ? detail.map((d: { msg: string }) => d.msg).join("; ")
      : (detail ?? `Request failed (${response.status})`);
    throw new Error(message);
  }
  return response.json();
}
