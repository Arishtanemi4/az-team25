// Shared 4-decimal formatter for every desirability-like number (D, d_gene, per-layer d, m_eff,
// rho_bar, similarity S). Four decimals is arithmetic precision, not measurement accuracy
// (CONSTRAINTS.md F5) -- the note that says so lives beside the ranked list, not in this helper.
export function formatScore(value: number | null | undefined): string {
  if (value === null || value === undefined) return "abstained";
  return value.toFixed(4);
}
