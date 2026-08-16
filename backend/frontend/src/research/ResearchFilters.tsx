export interface ResearchFilterValues {
  excludeProblematic: boolean;
  msiHigh: boolean | null;
  ploidyMin: number | null;
  ploidyMax: number | null;
  requireMetabolomics: boolean;
  requireMirna: boolean;
}

interface ResearchFiltersProps {
  value: ResearchFilterValues;
  onChange: (value: ResearchFilterValues) => void;
}

function nullableNumber(value: string): number | null {
  if (value.trim() === "") return null;
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : null;
}

// These controls map one-to-one to the already-native S02 request fields. They are optional
// cohort restrictions, not evidence weights or an algorithm selector.
export function ResearchFilters({ value, onChange }: ResearchFiltersProps) {
  function patch(changes: Partial<ResearchFilterValues>) {
    onChange({ ...value, ...changes });
  }

  return (
    <details className="research-filters">
      <summary>Research filters (optional)</summary>
      <p className="research-filter-note">
        These restrict the candidate cohort before native scoring. They do not change the desirability formula.
      </p>
      <div className="research-filter-grid">
        <label className="research-filter-checkbox">
          <input type="checkbox" checked={value.excludeProblematic} onChange={(event) => patch({ excludeProblematic: event.target.checked })} />
          Exclude lines flagged as problematic
        </label>
        <label className="research-filter-checkbox">
          <input type="checkbox" checked={value.requireMetabolomics} onChange={(event) => patch({ requireMetabolomics: event.target.checked })} />
          Require metabolomics availability
        </label>
        <label className="research-filter-checkbox">
          <input type="checkbox" checked={value.requireMirna} onChange={(event) => patch({ requireMirna: event.target.checked })} />
          Require miRNA availability
        </label>
        <label className="research-filter-select">
          MSI status
          <select value={value.msiHigh === null ? "any" : value.msiHigh ? "high" : "not-high"} onChange={(event) => patch({ msiHigh: event.target.value === "any" ? null : event.target.value === "high" })}>
            <option value="any">Any</option>
            <option value="high">MSI-high</option>
            <option value="not-high">Not MSI-high</option>
          </select>
        </label>
        <label className="research-filter-select">
          Minimum ploidy
          <input type="number" inputMode="decimal" step="any" value={value.ploidyMin ?? ""} onChange={(event) => patch({ ploidyMin: nullableNumber(event.target.value) })} />
        </label>
        <label className="research-filter-select">
          Maximum ploidy
          <input type="number" inputMode="decimal" step="any" value={value.ploidyMax ?? ""} onChange={(event) => patch({ ploidyMax: nullableNumber(event.target.value) })} />
        </label>
      </div>
    </details>
  );
}
