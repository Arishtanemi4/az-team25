import type { GeneMatch } from "../types/scoring";

const BASE_URL = import.meta.env.VITE_GENE_SERVICE_URL ?? "http://localhost:8000";

export async function searchGenes(query: string, limit = 20): Promise<GeneMatch[]> {
  const url = `${BASE_URL}/genes/search?q=${encodeURIComponent(query)}&limit=${limit}`;
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Gene search failed (${response.status})`);
  return response.json();
}

export async function listLineages(): Promise<string[]> {
  const response = await fetch(`${BASE_URL}/filters/lineages`);
  if (!response.ok) throw new Error(`Lineage lookup failed (${response.status})`);
  const body = await response.json();
  return body.lineages;
}

export async function listPrimaryDiseases(lineage: string | null): Promise<string[]> {
  const url = lineage
    ? `${BASE_URL}/filters/primary-diseases?lineage=${encodeURIComponent(lineage)}`
    : `${BASE_URL}/filters/primary-diseases`;
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Disease lookup failed (${response.status})`);
  const body = await response.json();
  return body.primary_diseases;
}
