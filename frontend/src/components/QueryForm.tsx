import { useEffect, useState } from "react";
import { listLineages, listPrimaryDiseases } from "../api/geneServiceClient";
import type { GeneMatch } from "../types/scoring";
import { FilterSelect } from "./FilterSelect";
import { GeneCombobox } from "./GeneCombobox";
import { QueryExpansionPanel } from "./QueryExpansionPanel";
import { ResearchFilters, type ResearchFilterValues } from "../research/ResearchFilters";

export interface QuerySubmission {
  inclusionGenes: GeneMatch[];
  exclusionGenes: GeneMatch[];
  lineage: string | null;
  primaryDisease: string | null;
  excludeProblematic: boolean;
  msiHigh: boolean | null;
  ploidyMin: number | null;
  ploidyMax: number | null;
  requireMetabolomics: boolean;
  requireMirna: boolean;
  topK: number;
}

interface QueryFormProps {
  onSubmit: (query: QuerySubmission) => void;
  isSubmitting: boolean;
}

export function QueryForm({ onSubmit, isSubmitting }: QueryFormProps) {
  const [inclusionGenes, setInclusionGenes] = useState<GeneMatch[]>([]);
  const [exclusionGenes, setExclusionGenes] = useState<GeneMatch[]>([]);
  const [lineage, setLineage] = useState<string | null>(null);
  const [primaryDisease, setPrimaryDisease] = useState<string | null>(null);
  const [topK, setTopK] = useState(10);
  const [researchFilters, setResearchFilters] = useState<ResearchFilterValues>({
    excludeProblematic: false,
    msiHigh: null,
    ploidyMin: null,
    ploidyMax: null,
    requireMetabolomics: false,
    requireMirna: false,
  });
  const [lineages, setLineages] = useState<string[]>([]);
  const [primaryDiseases, setPrimaryDiseases] = useState<string[]>([]);

  useEffect(() => {
    listLineages().then(setLineages).catch(() => setLineages([]));
  }, []);

  // The disease dropdown cascades from whichever lineage is chosen, so a researcher never sees
  // a disease that can't occur in the selected lineage.
  useEffect(() => {
    listPrimaryDiseases(lineage).then(setPrimaryDiseases).catch(() => setPrimaryDiseases([]));
  }, [lineage]);

  const exclusionOnly = inclusionGenes.length === 0 && exclusionGenes.length > 0;
  const invalidPloidyRange = researchFilters.ploidyMin !== null
    && researchFilters.ploidyMax !== null
    && researchFilters.ploidyMin > researchFilters.ploidyMax;
  const canSubmit = inclusionGenes.length > 0 && !invalidPloidyRange && !isSubmitting;

  function addSuggestedGene(gene: GeneMatch) {
    // Never auto-added -- this only runs when the researcher clicks "Add" on a suggestion
    // (QueryExpansionPanel). Skip if it's already in either list.
    const alreadyPresent = [...inclusionGenes, ...exclusionGenes].some(
      (g) => g.ensembl_id === gene.ensembl_id
    );
    if (!alreadyPresent) setInclusionGenes([...inclusionGenes, gene]);
  }

  return (
    <form
      className="query-form"
      onSubmit={(e) => {
        e.preventDefault();
        if (!canSubmit) return;
        onSubmit({
          inclusionGenes,
          exclusionGenes,
          lineage,
          primaryDisease,
          ...researchFilters,
          topK,
        });
      }}
    >
      <GeneCombobox
        label="Inclusion genes (must be present/active)"
        selected={inclusionGenes}
        onChange={setInclusionGenes}
        excluded={exclusionGenes}
      />
      <QueryExpansionPanel inclusionGenes={inclusionGenes} onAddGene={addSuggestedGene} />
      <GeneCombobox
        label="Exclusion genes (must be absent/inactive)"
        selected={exclusionGenes}
        onChange={setExclusionGenes}
        excluded={inclusionGenes}
      />
      <div className="filter-row">
        <FilterSelect label="Lineage" options={lineages} value={lineage} onChange={setLineage} />
        <FilterSelect
          label="Primary disease"
          options={primaryDiseases}
          value={primaryDisease}
          onChange={setPrimaryDisease}
        />
        <div className="combobox">
          <label className="combobox-label" htmlFor="top-k">
            Ranked lines to show
          </label>
          <input
            id="top-k"
            type="number"
            min={1}
            max={100}
            value={topK}
            onChange={(e) => setTopK(Number(e.target.value))}
          />
        </div>
      </div>
      <ResearchFilters value={researchFilters} onChange={setResearchFilters} />
      <button type="submit" disabled={!canSubmit}>
        {isSubmitting ? "Scoring cell lines..." : "Rank cell lines"}
      </button>
      {exclusionOnly && !isSubmitting && (
        <p className="form-notice" role="alert">
          An exclusion-only query is not supported: add at least one inclusion gene. Exclusions are optional criteria.
        </p>
      )}
      {invalidPloidyRange && !isSubmitting && (
        <p className="form-notice" role="alert">Minimum ploidy cannot exceed maximum ploidy.</p>
      )}
      {!canSubmit && !isSubmitting && !exclusionOnly && !invalidPloidyRange && (
        <p className="hint">Add at least one inclusion gene to submit. Exclusion genes are optional criteria.</p>
      )}
    </form>
  );
}
