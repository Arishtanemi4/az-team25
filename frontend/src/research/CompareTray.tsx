interface CompareTrayProps {
  queryId: string | null | undefined;
  selectedModelIds: string[];
  notice: string | null;
  isLoading: boolean;
  error: string | null;
  onRemove: (modelId: string) => void;
  onCompare: () => void;
  onDownload: (format: "json" | "csv") => void;
}

export function CompareTray({ queryId, selectedModelIds, notice, isLoading, error, onRemove, onCompare, onDownload }: CompareTrayProps) {
  const canCompare = Boolean(queryId) && selectedModelIds.length >= 2 && selectedModelIds.length <= 3;
  return (
    <section className="compare-tray" aria-label="Compare selected cell lines">
      <h3>Compare cell lines</h3>
      {!queryId ? (
        <p className="context-unavailable">Research comparison needs a trusted snapshot. Enable the local research extension and rerun this query.</p>
      ) : (
        <>
          <p className="compare-help">Select two or three distinct ModelIDs in the order you want columns shown. Comparison preserves native units and states; it does not choose a winner.</p>
          {selectedModelIds.length === 0 ? <p className="context-muted">No models selected.</p> : (
            <ol className="compare-selection-list">
              {selectedModelIds.map((modelId) => <li key={modelId}><code>{modelId}</code><button type="button" onClick={() => onRemove(modelId)} aria-label={`Remove ${modelId} from comparison`}>Remove</button></li>)}
            </ol>
          )}
          {notice && <p className="compare-notice" role="status">{notice}</p>}
          {error && <p className="context-unavailable" role="alert">{error}</p>}
          <div className="compare-actions">
            <button type="button" onClick={onCompare} disabled={!canCompare || isLoading}>{isLoading ? "Building comparison..." : "Compare selected"}</button>
            <button type="button" onClick={() => onDownload("json")} disabled={!canCompare || isLoading}>Download complete research JSON</button>
            <button type="button" onClick={() => onDownload("csv")} disabled={!canCompare || isLoading}>Download comparison CSV</button>
          </div>
          {selectedModelIds.length < 2 && <p className="context-muted">Select one more model to compare.</p>}
        </>
      )}
    </section>
  );
}
