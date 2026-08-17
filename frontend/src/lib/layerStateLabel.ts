// Shared human-readable labels for a layer's native `state` string. Extracted from
// PerGeneEvidenceTable.tsx (V7-8) so it and the result-path flowchart agree by construction
// rather than by copy-paste.
export function stateLabel(state: string): string {
  if (state === "measured") return "measured";
  if (state === "measured_absent") return "measured; no event observed";
  if (state === "not_assayed") return "not assayed (unknown, not zero)";
  if (state === "non_detected") return "not detected";
  if (state === "excluded_pan_essential") return "excluded (pan-essential gate)";
  return `unknown native state: ${state}`;
}
