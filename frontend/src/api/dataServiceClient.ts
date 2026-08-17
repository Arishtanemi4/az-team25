import type {
  CellLineDiseaseResponse,
  CoverageResponse,
  EdaConcordanceResponse,
  GeneCellLineResponse,
  GeneDiseaseResponse,
  LineageCoverageResponse,
  PreviewResponse,
  SchemaResponse,
  TablesResponse,
  ValidationStudiesResponse,
  ValidationStudyResponse,
} from "../types/data";

const BASE_URL = import.meta.env.VITE_DATA_SERVICE_URL ?? "http://localhost:8000";

// Reads the FastAPI error body so a caller can show *why* a request was refused (e.g. a
// row limit above the cap gets a real 422 with a validation message), not just a status code.
async function readErrorDetail(response: Response): Promise<string> {
  try {
    const body = await response.json();
    if (typeof body?.detail === "string") return body.detail;
    if (Array.isArray(body?.detail)) {
      return body.detail.map((item: { msg?: string }) => item.msg ?? JSON.stringify(item)).join("; ");
    }
  } catch {
    // response body wasn't JSON -- fall through to the status-only message
  }
  return `request failed (${response.status})`;
}

export async function getSchema(): Promise<SchemaResponse> {
  const response = await fetch(`${BASE_URL}/data/schema`);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}

export async function getTables(): Promise<TablesResponse> {
  const response = await fetch(`${BASE_URL}/data/tables`);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}

export async function getTablePreview(
  table: string,
  limit: number,
  offset: number,
): Promise<PreviewResponse> {
  const url = `${BASE_URL}/data/tables/${encodeURIComponent(table)}/preview?limit=${limit}&offset=${offset}`;
  const response = await fetch(url);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}

export async function getRelationshipsGeneCellLine(
  ensemblId: string,
  layer: string,
  limit: number,
): Promise<GeneCellLineResponse> {
  const url = `${BASE_URL}/data/relationships/gene-cell-line?ensembl_id=${encodeURIComponent(ensemblId)}&layer=${encodeURIComponent(layer)}&limit=${limit}`;
  const response = await fetch(url);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}

export async function getRelationshipsCellLineDisease(
  lineage: string | null,
): Promise<CellLineDiseaseResponse> {
  const url = lineage
    ? `${BASE_URL}/data/relationships/cell-line-disease?lineage=${encodeURIComponent(lineage)}`
    : `${BASE_URL}/data/relationships/cell-line-disease`;
  const response = await fetch(url);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}

export async function getRelationshipsGeneDisease(
  ensemblId: string,
  layer: string,
  minN: number,
): Promise<GeneDiseaseResponse> {
  const url = `${BASE_URL}/data/relationships/gene-disease?ensembl_id=${encodeURIComponent(ensemblId)}&layer=${encodeURIComponent(layer)}&min_n=${minN}`;
  const response = await fetch(url);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}

export async function getEdaCoverage(): Promise<CoverageResponse> {
  const response = await fetch(`${BASE_URL}/data/eda/coverage`);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}

export async function getEdaLineageCoverage(): Promise<LineageCoverageResponse> {
  const response = await fetch(`${BASE_URL}/data/eda/lineage-coverage`);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}

export async function getEdaConcordance(): Promise<EdaConcordanceResponse> {
  const response = await fetch(`${BASE_URL}/data/eda/concordance`);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}

export async function getValidationStudies(): Promise<ValidationStudiesResponse> {
  const response = await fetch(`${BASE_URL}/data/validation/studies`);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}

export async function getValidationStudy(study: string): Promise<ValidationStudyResponse> {
  const response = await fetch(`${BASE_URL}/data/validation/${encodeURIComponent(study)}`);
  if (!response.ok) throw new Error(await readErrorDetail(response));
  return response.json();
}
