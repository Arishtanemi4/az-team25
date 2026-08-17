// A cell line's human name is display metadata joined from data/processed/cell_lines.csv by
// backend/scoring_service (see types/scoring.ts's cell_line_name comment) -- it is optional, so a
// line with no name renders its bare ModelID rather than "undefined (...)".
export function formatCellLineLabel(modelId: string, cellLineName?: string | null): string {
  return cellLineName ? `${cellLineName} (${modelId})` : modelId;
}
