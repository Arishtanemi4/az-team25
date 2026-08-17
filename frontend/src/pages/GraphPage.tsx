import { useEffect, useRef, useState } from "react";
import cytoscape from "cytoscape";
import { GeneSearchSelect } from "../components/GeneSearchSelect";
import { getGraphNeighborhood } from "../api/ragServiceClient";
import type { GraphNeighborhoodResponse } from "../types/graph";
import type { GeneMatch } from "../types/scoring";
import "./graphPage.css";

const MAX_NEIGHBORS_CAP = 200; // the server's own le=200 cap -- shown, never silently clamped (F3)
const DEFAULT_MAX_NEIGHBORS = 50;

// Exact edge_type/source_db pairs rag/build_knowledge_graph.py writes -- the exhaustive legend,
// not a guess (confirmed by reading that script directly).
const EDGE_LEGEND: { edge_type: string; label: string; color: string }[] = [
  { edge_type: "gene_in_pathway", label: "Gene in pathway (Reactome)", color: "#1c7ed6" },
  { edge_type: "pathway_parent_of", label: "Pathway hierarchy (Reactome)", color: "#74b9ff" },
  { edge_type: "protein_interaction", label: "Protein interaction (STRING v12)", color: "#2f9e44" },
  { edge_type: "curated_interaction", label: "Curated interaction (BioGRID)", color: "#f59f00" },
  { edge_type: "fusion_partner", label: "Fusion partner (fusions.csv)", color: "#e64980" },
  { edge_type: "codependency", label: "Co-dependency (dependency.csv)", color: "#7048e8" },
];
const EDGE_COLOR_BY_TYPE = new Map(EDGE_LEGEND.map((e) => [e.edge_type, e.color]));

const NODE_LEGEND = [
  { node_type: "gene", label: "Gene", color: "#aa3bff" },
  { node_type: "pathway", label: "Pathway (Reactome)", color: "#12b886" },
];
const NODE_COLOR_BY_TYPE = new Map(NODE_LEGEND.map((n) => [n.node_type, n.color]));

// Builds/destroys the Cytoscape instance in a useEffect keyed on `data`, mirroring
// ChartCanvas.tsx's build-then-destroy-on-change lifecycle for Chart.js. Cytoscape needs a
// graph layout (Chart.js has none), which is the one reason this page pulls in a second chart
// library rather than reusing ChartCanvas (PRODUCT_SURFACE.md SS1 decision 1).
function GraphCanvas({ data }: { data: GraphNeighborhoodResponse }) {
  const containerRef = useRef<HTMLDivElement>(null);
  const cyRef = useRef<cytoscape.Core | null>(null);

  useEffect(() => {
    if (!containerRef.current) return;
    cyRef.current?.destroy();

    const elements = [
      ...data.nodes.map((n) => ({
        data: { id: n.node_id, label: n.display_name, node_type: n.node_type },
      })),
      ...data.edges.map((e, i) => ({
        data: { id: `edge-${i}`, source: e.source, target: e.target, edge_type: e.edge_type },
      })),
    ];

    cyRef.current = cytoscape({
      container: containerRef.current,
      elements,
      // animate: false -- cose's default animated layout runs its own requestAnimationFrame
      // loop; switching genes quickly destroys the Cytoscape instance while that loop is still
      // ticking, which throws "Cannot read properties of null (reading 'notify')" from inside
      // cytoscape.js itself once it tries to notify the now-destroyed core. A static explanatory
      // neighbourhood view doesn't need the animation, so removing it removes the race entirely.
      layout: { name: "cose", animate: false },
      style: [
        {
          selector: "node",
          style: {
            "background-color": (ele) => NODE_COLOR_BY_TYPE.get(ele.data("node_type")) ?? "#868e96",
            label: "data(label)",
            "font-size": 8,
            width: 14,
            height: 14,
          },
        },
        {
          selector: "edge",
          style: {
            "line-color": (ele) => EDGE_COLOR_BY_TYPE.get(ele.data("edge_type")) ?? "#ced4da",
            width: 1,
            "curve-style": "haystack",
          },
        },
      ],
    });

    return () => {
      cyRef.current?.destroy();
      cyRef.current = null;
    };
  }, [data]);

  return <div ref={containerRef} className="graph-canvas" />;
}

function GraphPage() {
  const [gene, setGene] = useState<GeneMatch | null>(null);
  const [hops, setHops] = useState<1 | 2>(1);
  const [maxNeighbors, setMaxNeighbors] = useState(DEFAULT_MAX_NEIGHBORS);
  const [data, setData] = useState<GraphNeighborhoodResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    if (!gene) {
      setData(null);
      return;
    }
    let cancelled = false;
    setLoading(true);
    setError(null);
    // Pass the symbol, not the Ensembl ID -- the backend echoes this same token back as the
    // seed node's display_name (graph_service.get_neighborhood), so passing the ID would label
    // the seed node "ENSG00000146648" on the canvas instead of "EGFR". resolve_symbol accepts
    // either form; the symbol is what a person actually reads.
    getGraphNeighborhood(gene.symbol, hops, maxNeighbors)
      .then((res) => {
        if (!cancelled) setData(res);
      })
      .catch((err: Error) => {
        if (!cancelled) {
          setData(null);
          setError(err.message);
        }
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, [gene, hops, maxNeighbors]);

  return (
    <div className="app-shell">
      <header>
        <h1>Knowledge graph</h1>
        <p className="subtitle">
          A seed gene's k-hop neighbourhood — Reactome pathways, STRING and BioGRID interactions, plus
          this project's own fusion-partner and co-dependency edges. This page explains relationships; it
          never scores, ranks, or feeds a suggestion into a query on its own.
        </p>
      </header>

      <div className="explore-controls">
        <GeneSearchSelect label="Seed gene" selected={gene} onChange={setGene} />
        <label className="query-field">
          Hops
          <select value={hops} onChange={(e) => setHops(Number(e.target.value) as 1 | 2)}>
            <option value={1}>1</option>
            <option value={2}>2</option>
          </select>
        </label>
        <label className="query-field">
          Max neighbours
          <input
            type="number"
            min={1}
            value={maxNeighbors}
            onChange={(e) =>
              setMaxNeighbors(Number.isFinite(e.target.valueAsNumber) && e.target.valueAsNumber > 0
                ? Math.floor(e.target.valueAsNumber)
                : 1)
            }
          />
        </label>
        <span className="query-field-note">
          The API caps this at {MAX_NEIGHBORS_CAP} — try a larger number to see it refused. Hops is capped
          at 2: EGFR alone has 3,270 one-hop and 24,257 two-hop neighbours, so a third hop would be a
          browser-killing response, not a useful view.
        </span>
      </div>

      {!gene && <p className="explore-placeholder">Pick a seed gene to see its knowledge-graph neighbourhood.</p>}
      {loading && <p>Loading…</p>}
      {error && <p className="uml-error">{error}</p>}

      {data && !data.available && (
        <div className="diagnostic-banner">
          <p>{data.reason}</p>
          {data.build_command && (
            <p>
              Build it with: <code>{data.build_command}</code>
            </p>
          )}
        </div>
      )}

      {data && data.available && data.nodes.length === 0 && (
        <p className="explore-placeholder">{data.reason}</p>
      )}

      {data && data.available && data.nodes.length > 0 && (
        <>
          <p className="query-denominator">
            <strong>{data.nodes.length.toLocaleString()}</strong> nodes,{" "}
            <strong>{data.edges.length.toLocaleString()}</strong> edges
            {data.truncated && <> — capped at {maxNeighbors} neighbours; more exist but are not shown</>}
          </p>
          <GraphCanvas data={data} />
          <div className="graph-legend">
            <div>
              <h3>Node types</h3>
              {NODE_LEGEND.map((n) => (
                <span key={n.node_type} className="graph-legend-item">
                  <span className="graph-legend-swatch" style={{ background: n.color }} /> {n.label}
                </span>
              ))}
            </div>
            <div>
              <h3>Edge types</h3>
              {EDGE_LEGEND.map((e) => (
                <span key={e.edge_type} className="graph-legend-item">
                  <span className="graph-legend-swatch" style={{ background: e.color }} /> {e.label}
                </span>
              ))}
            </div>
          </div>
        </>
      )}
    </div>
  );
}

export default GraphPage;
