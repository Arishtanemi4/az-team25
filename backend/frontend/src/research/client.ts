import type { AlternativesResponse, ComparisonResponse, ContextResponse } from "./types";

const BASE_URL = import.meta.env.VITE_RESEARCH_API_URL
  ?? import.meta.env.VITE_SCORING_SERVICE_URL
  ?? "http://localhost:8000";

export class ResearchRequestError extends Error {}

async function researchJson<T>(path: string, options: RequestInit = {}): Promise<T> {
  const response = await fetch(`${BASE_URL}${path}`, options);
  if (!response.ok) {
    const problem = await response.json().catch(() => null);
    throw new ResearchRequestError(problem?.detail ?? `Research request failed (${response.status}).`);
  }
  return response.json() as Promise<T>;
}

// Trusted query/model IDs come only from the current native result. This client never sends
// scores, native evidence, or local paths to the extension API.
export async function getBiologicalContext(queryId: string, modelId: string, signal?: AbortSignal): Promise<ContextResponse> {
  const payload = await researchJson<ContextResponse>(
    `/research/queries/${encodeURIComponent(queryId)}/models/${encodeURIComponent(modelId)}/context`, { signal },
  );
  if (payload.query_id !== queryId || payload.model_id !== modelId) {
    throw new ResearchRequestError("The biological context response did not match the active result.");
  }
  return payload;
}

export async function getAlternatives(queryId: string, modelId: string, signal?: AbortSignal): Promise<AlternativesResponse> {
  const payload = await researchJson<AlternativesResponse>(
    `/research/queries/${encodeURIComponent(queryId)}/models/${encodeURIComponent(modelId)}/alternatives`, { signal },
  );
  if (payload.query_id !== queryId || payload.anchor_id !== modelId) {
    throw new ResearchRequestError("The alternatives response did not match the active result.");
  }
  return payload;
}

export async function compareModels(queryId: string, modelIds: string[], signal?: AbortSignal): Promise<ComparisonResponse> {
  const payload = await researchJson<ComparisonResponse>(
    `/research/queries/${encodeURIComponent(queryId)}/compare`,
    { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ model_ids: modelIds }), signal },
  );
  if (payload.query_id !== queryId || payload.selected_model_ids.join("\u0000") !== modelIds.join("\u0000")) {
    throw new ResearchRequestError("The comparison response did not match the active selection.");
  }
  return payload;
}

export async function downloadResearchExport(queryId: string, modelIds: string[], format: "json" | "csv"): Promise<Blob> {
  const response = await fetch(`${BASE_URL}/research/queries/${encodeURIComponent(queryId)}/export-comparison`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ model_ids: modelIds, format }),
  });
  if (!response.ok) {
    const problem = await response.json().catch(() => null);
    throw new ResearchRequestError(problem?.detail ?? `Research export failed (${response.status}).`);
  }
  return response.blob();
}
