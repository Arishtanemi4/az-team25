import { useEffect, useState } from "react";
import { getBiologicalContext, ResearchRequestError } from "./client";
import type { BoundedGraphEntries, ContextResponse, EventPacket, GraphEdge } from "./types";
import { formatScore } from "../lib/formatScore";

interface BiologicalContextPanelProps {
  queryId: string | null | undefined;
  modelId: string;
  researchStatus?: string;
  researchReason?: string;
}

function displayValue(value: unknown): string {
  if (value === null || value === undefined || value === "") return "not supplied";
  return String(value);
}

function EventEvidence({ packet }: { packet: EventPacket }) {
  if (packet.status !== "ok") return <p className="context-muted">Events: {packet.status}; {packet.reason ?? "no further detail"}</p>;
  return (
    <div className="context-event-records">
      <p>Observed event rows: {packet.total_count ?? packet.rows?.length ?? 0}{packet.truncation_note ? ` — ${packet.truncation_note}` : ""}</p>
      {(packet.rows ?? []).map((row, index) => (
        <code key={index}>{JSON.stringify(row)}</code>
      ))}
    </div>
  );
}

function GraphEntries({ title, entries }: { title: string; entries: BoundedGraphEntries }) {
  return (
    <section className="graph-entry-section">
      <h5>{title} ({entries.displayed.length} shown of {entries.total_count})</h5>
      {entries.truncated && <p className="context-muted">Display is capped; total source records remain in the research export.</p>}
      {entries.displayed.length === 0 ? <p className="context-muted">No retrieved records are displayed.</p> : (
        <div className="context-table-scroll">
          <table className="context-table">
            <thead><tr><th>Source</th><th>Type</th><th>Neighbour / source ID</th><th>Evidence detail</th></tr></thead>
            <tbody>{entries.displayed.map((entry, index) => <GraphEntry key={`${entry.source_db}-${entry.edge_type}-${entry.neighbour_id}-${index}`} entry={entry} />)}</tbody>
          </table>
        </div>
      )}
    </section>
  );
}

function GraphEntry({ entry }: { entry: GraphEdge }) {
  const sourceId = entry.neighbour_id;
  return (
    <tr>
      <td>{displayValue(entry.source_db)}</td>
      <td>{displayValue(entry.edge_type)}</td>
      <td>{entry.url ? <a href={entry.url} target="_blank" rel="noreferrer">{entry.display_name ?? sourceId}</a> : <><span>{entry.display_name ?? sourceId}</span><code>{sourceId}</code></>}</td>
      <td>{displayValue(entry.evidence_detail)}</td>
    </tr>
  );
}

function ContextContent({ context }: { context: ContextResponse }) {
  const metadata = context.model_metadata ?? {};
  const measured = context.measured_context;
  const graph = context.graph_context;
  return (
    <div className="biological-context-content">
      <section className="context-section">
        <h4>Snapshot-pinned model metadata</h4>
        <dl className="context-metadata">
          <div><dt>Model ID</dt><dd>{context.model_id}</dd></div>
          <div><dt>Cell-line name</dt><dd>{displayValue(metadata.cell_line_name)}</dd></div>
          <div><dt>Lineage</dt><dd>{displayValue(metadata.lineage)}</dd></div>
          <div><dt>Native partition</dt><dd>{measured.partition_bucket}</dd></div>
          <div><dt>Problematic flag</dt><dd>{displayValue(measured.model_context.is_problematic)}</dd></div>
        </dl>
      </section>

      <section className="context-section">
        <h4>Measured evidence (native values, units, and states)</h4>
        <p className="context-boundary">This section reuses the query snapshot. It does not rerun ranking or calculate a second score.</p>
        {measured.genes.map((gene) => (
          <article className="context-gene" key={gene.ensembl_id}>
            <h5>{gene.symbol} <code>{gene.ensembl_id}</code> <span className="context-role">{gene.role}</span></h5>
            <p>Native per-gene desirability: {formatScore(gene.d_gene)}; class: {gene.gene_class}.</p>
            <div className="context-table-scroll">
              <table className="context-table">
                <thead><tr><th>Layer</th><th>State</th><th>Value / unit</th><th>Native d</th><th>Caveats</th></tr></thead>
                <tbody>{Object.entries(gene.layers).map(([layer, detail]) => (
                  <tr key={layer}>
                    <td>{layer}</td><td>{detail.state}</td><td>{detail.value === null ? "not reported" : `${detail.value} ${detail.unit}`}</td><td>{formatScore(detail.d)}</td>
                    <td>{detail.caveats?.map((caveat) => caveat.message).join("; ") || "none reported"}</td>
                  </tr>
                ))}</tbody>
              </table>
            </div>
            {Object.values(gene.layers).flatMap((detail) => detail.events ? [detail.events] : []).map((packet, index) => <EventEvidence key={index} packet={packet} />)}
            {gene.evidence_caveats.length > 0 && <ul className="context-caveats">{gene.evidence_caveats.map((caveat, index) => <li key={index}>{caveat.message}</li>)}</ul>}
          </article>
        ))}
      </section>

      <section className="context-section graph-context-section">
        <h4>Generic pathway and interaction reference context</h4>
        <p className="context-boundary">Membership and graph edges are generic reference associations, not activation, enrichment, causality, or cell-line sensitivity.</p>
        <dl className="context-metadata">
          <div><dt>Graph source mode</dt><dd>{graph.graph_source.mode}</dd></div>
          <div><dt>Scope</dt><dd>{displayValue(graph.graph_source.scope)}</dd></div>
        </dl>
        {graph.graph_source.reason && <p className="context-unavailable">{graph.graph_source.reason}</p>}
        {Object.entries(graph.genes).map(([geneId, gene]) => (
          <article className="context-gene" key={geneId}>
            <h5><code>{geneId}</code> <span className="context-role">{gene.role}</span></h5>
            <p>{gene.claim_type}{gene.reason ? `: ${gene.reason}` : ""}</p>
            <GraphEntries title="Direct reference edges" entries={gene.direct_edges} />
            <GraphEntries title="Pathway membership" entries={gene.pathway_membership} />
          </article>
        ))}
        {graph.limitations.length > 0 && <ul className="context-caveats">{graph.limitations.map((limitation, index) => <li key={index}>{limitation.message}</li>)}</ul>}
      </section>
    </div>
  );
}

export function BiologicalContextPanel({ queryId, modelId, researchStatus, researchReason }: BiologicalContextPanelProps) {
  const [isOpen, setIsOpen] = useState(false);
  const [context, setContext] = useState<ContextResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(false);

  // A new native query immediately closes and clears this panel. The request's AbortSignal and
  // query/model equality checks prevent a late response from attaching to a newer result.
  useEffect(() => {
    setIsOpen(false);
    setContext(null);
    setError(null);
  }, [queryId, modelId]);

  useEffect(() => {
    if (!isOpen || !queryId) return;
    const controller = new AbortController();
    setIsLoading(true);
    setError(null);
    getBiologicalContext(queryId, modelId, controller.signal)
      .then((response) => {
        if (response.query_id === queryId && response.model_id === modelId) setContext(response);
      })
      .catch((requestError: unknown) => {
        if (controller.signal.aborted) return;
        setContext(null);
        setError(requestError instanceof ResearchRequestError ? requestError.message : "Biological context could not be loaded.");
      })
      .finally(() => {
        if (!controller.signal.aborted) setIsLoading(false);
      });
    return () => controller.abort();
  }, [isOpen, queryId, modelId]);

  if (!queryId) {
    return <section className="biological-context-panel"><h3>Biological context</h3><p className="context-unavailable">{researchReason ?? (researchStatus === "unavailable" ? "The research snapshot was unavailable; native evidence remains available above." : "This result has no trusted research snapshot. Enable the local research extension and run the query again.")}</p></section>;
  }

  return (
    <section className="biological-context-panel">
      <button type="button" className="context-toggle" onClick={() => setIsOpen((value) => !value)} aria-expanded={isOpen}>
        {isOpen ? "Hide biological context" : "Show biological context"}
      </button>
      {isOpen && <>
        {isLoading && <p className="context-muted" role="status">Loading deterministic biological context; native result details remain above.</p>}
        {error && <p className="context-unavailable" role="alert">{error} Native result details remain available.</p>}
        {context && !isLoading && <ContextContent context={context} />}
      </>}
    </section>
  );
}
