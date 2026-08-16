import type { RankResponse } from "../types/scoring";

const BASE_URL = import.meta.env.VITE_SCORING_SERVICE_URL ?? "http://localhost:8000";

export interface RankRequestBody {
  inclusion_genes: string[];
  exclusion_genes: string[];
  lineage: string | null;
  primary_disease: string | null;
  exclude_problematic: boolean;
  msi_high: boolean | null;
  ploidy_min: number | null;
  ploidy_max: number | null;
  require_metabolomics: boolean;
  require_mirna: boolean;
  top_k: number;
}

export interface RankStageEvent {
  label: string;
  index: number;
  total: number;
}

export interface RankStreamCallbacks {
  onStage?: (stage: RankStageEvent) => void;
}

export class RankRequestError extends Error {}

// /rank streams Server-Sent Events (a `stage` event per pipeline checkpoint, then one `result`
// event) rather than plain JSON, so the response can't just be awaited as .json() -- native
// EventSource can't POST a body either, so this hand-rolls SSE parsing over a fetch() reader.
export async function rankCellLinesStream(
  body: RankRequestBody,
  callbacks: RankStreamCallbacks = {}
): Promise<RankResponse> {
  const response = await fetch(`${BASE_URL}/rank`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });

  if (!response.ok) {
    // Pre-stream validation failure (422) -- same error contract as before, since headers
    // arrived before any SSE framing began.
    const problem = await response.json().catch(() => null);
    const detail = problem?.detail;
    const message = Array.isArray(detail)
      ? detail.map((d: { msg: string }) => d.msg).join("; ")
      : (detail ?? `Ranking request failed (${response.status})`);
    throw new RankRequestError(message);
  }
  if (!response.body) {
    throw new RankRequestError("Streaming is not supported by this browser.");
  }

  const reader = response.body.getReader();
  const decoder = new TextDecoder();
  let buffer = "";
  let result: RankResponse | null = null;

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;
    buffer += decoder.decode(value, { stream: true });

    let boundary = buffer.indexOf("\n\n");
    while (boundary !== -1) {
      const chunk = buffer.slice(0, boundary);
      buffer = buffer.slice(boundary + 2);
      const lines = chunk.split("\n");
      const eventLine = lines.find((l) => l.startsWith("event: "));
      const dataLine = lines.find((l) => l.startsWith("data: "));
      const event = eventLine?.slice("event: ".length) ?? "message";
      const data = dataLine ? JSON.parse(dataLine.slice("data: ".length)) : null;

      if (event === "stage") {
        callbacks.onStage?.(data as RankStageEvent);
      } else if (event === "result") {
        result = data as RankResponse;
      } else if (event === "error") {
        throw new RankRequestError((data as { message: string }).message);
      }

      boundary = buffer.indexOf("\n\n");
    }
  }

  if (!result) {
    throw new RankRequestError("Ranking stream ended before a result was received.");
  }
  return result;
}
