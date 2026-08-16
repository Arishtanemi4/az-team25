import { useEffect, useState } from "react";
import { getTablePreview, getTables } from "../api/dataServiceClient";
import type { PreviewResponse, TableSummary } from "../types/data";
import "./dataPages.css";

const DEFAULT_LIMIT = 50;

function DataQueryPage() {
  const [tables, setTables] = useState<TableSummary[] | null>(null);
  const [tablesError, setTablesError] = useState<string | null>(null);
  const [selectedTable, setSelectedTable] = useState<string>("");
  const [limit, setLimit] = useState(DEFAULT_LIMIT);
  const [offset, setOffset] = useState(0);
  const [preview, setPreview] = useState<PreviewResponse | null>(null);
  const [previewError, setPreviewError] = useState<string | null>(null);
  const [previewLoading, setPreviewLoading] = useState(false);

  useEffect(() => {
    getTables()
      .then((data) => {
        setTables(data.tables);
        if (data.tables.length > 0) setSelectedTable(data.tables[0].name);
      })
      .catch((err: Error) => setTablesError(err.message));
  }, []);

  useEffect(() => {
    if (!selectedTable) return;
    let cancelled = false;
    setPreviewLoading(true);
    setPreviewError(null);
    getTablePreview(selectedTable, limit, offset)
      .then((data) => {
        if (!cancelled) setPreview(data);
      })
      .catch((err: Error) => {
        if (!cancelled) {
          setPreview(null);
          setPreviewError(err.message);
        }
      })
      .finally(() => {
        if (!cancelled) setPreviewLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, [selectedTable, limit, offset]);

  function handleTableChange(name: string) {
    setSelectedTable(name);
    setOffset(0);
  }

  function handleLimitChange(value: number) {
    // Deliberately no client-side upper clamp -- typing a value above the API's cap (500) is
    // what actually exercises the real 422 refusal, not a client-side stand-in for it (F3).
    setLimit(Number.isFinite(value) && value > 0 ? Math.floor(value) : 1);
    setOffset(0);
  }

  return (
    <div className="app-shell">
      <header>
        <h1>Data query</h1>
        <p className="subtitle">Pick a table and a bounded row limit to see a sample of rows.</p>
      </header>

      {tablesError && <p className="uml-error">Could not load the table list: {tablesError}</p>}

      {tables && (
        <div className="query-controls">
          <label className="query-field">
            Table
            <select value={selectedTable} onChange={(e) => handleTableChange(e.target.value)}>
              {tables.map((t) => (
                <option key={t.name} value={t.name}>
                  {t.name} ({t.row_count.toLocaleString()} rows)
                </option>
              ))}
            </select>
          </label>
          <label className="query-field">
            Row limit
            <input
              type="number"
              min={1}
              value={limit}
              onChange={(e) => handleLimitChange(e.target.valueAsNumber)}
            />
          </label>
          <span className="query-field-note">The API caps this at 500 — try a larger number to see it refused.</span>
        </div>
      )}

      {previewLoading && <p>Loading preview…</p>}
      {previewError && <p className="uml-error">Request refused: {previewError}</p>}

      {preview && !previewError && (
        <>
          <p className="query-denominator">
            Showing <strong>{preview.returned_rows.toLocaleString()}</strong> of{" "}
            <strong>{preview.total_rows.toLocaleString()}</strong> rows — grain: {preview.grain}
          </p>
          <div className="query-table-wrap">
            <table className="query-table">
              <thead>
                <tr>
                  {preview.columns.map((col) => (
                    <th key={col}>{col}</th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {preview.rows.map((row, i) => (
                  <tr key={i}>
                    {preview.columns.map((col) => (
                      <td key={col}>{row[col] === null || row[col] === undefined ? "" : String(row[col])}</td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div className="query-pagination">
            <button type="button" disabled={offset === 0} onClick={() => setOffset(Math.max(0, offset - limit))}>
              Previous
            </button>
            <span>Offset {offset.toLocaleString()}</span>
            <button
              type="button"
              disabled={offset + preview.returned_rows >= preview.total_rows}
              onClick={() => setOffset(offset + limit)}
            >
              Next
            </button>
          </div>
        </>
      )}
    </div>
  );
}

export default DataQueryPage;
