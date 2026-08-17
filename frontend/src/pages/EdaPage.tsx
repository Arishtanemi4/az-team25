import { useEffect, useMemo, useState } from "react";
import type { ChartConfiguration } from "chart.js/auto";
import { ChartCanvas } from "../components/charts/ChartCanvas";
import { getEdaConcordance, getEdaCoverage, getEdaLineageCoverage } from "../api/dataServiceClient";
import type { ConcordanceSummary, CoverageResponse, EdaConcordanceResponse, LineageCoverageResponse } from "../types/data";
import "./explorePages.css";

const STATE_COLORS: Record<string, string> = {
  measured: "#2f9e44",
  not_assayed: "#adb5bd",
  measured_absent: "#e64980",
  true: "#2f9e44",
  false: "#adb5bd",
};

function stripSuffix(key: string): string {
  return key.replace(/_state$/, "").replace(/_available$/, "").replace(/_/g, " ");
}

function CoverageSection() {
  const [data, setData] = useState<CoverageResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    let cancelled = false;
    getEdaCoverage()
      .then((res) => {
        if (!cancelled) setData(res);
      })
      .catch((err: Error) => {
        if (!cancelled) setError(err.message);
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, []);

  const { stateChart, availableChart, rows } = useMemo(() => {
    if (!data) return { stateChart: null, availableChart: null, rows: [] as { layer: string; n: number; states: Record<string, number> }[] };
    const entries = Object.entries(data.coverage_state_counts);
    const stateEntries = entries.filter(([key]) => key.endsWith("_state"));
    const availableEntries = entries.filter(([key]) => key.endsWith("_available"));

    function buildStacked(pairs: [string, Record<string, number>][], stateKeys: string[]) {
      if (pairs.length === 0) return null;
      return {
        type: "bar" as const,
        data: {
          labels: pairs.map(([key]) => stripSuffix(key)),
          datasets: stateKeys.map((state) => ({
            label: state,
            data: pairs.map(([, counts]) => counts[state] ?? 0),
            backgroundColor: STATE_COLORS[state] ?? "#495057",
          })),
        },
        options: {
          responsive: true,
          scales: { x: { stacked: true }, y: { stacked: true, title: { display: true, text: "n cell lines" } } },
        },
      } as ChartConfiguration;
    }

    const rows = entries.map(([layer, counts]) => ({
      layer: stripSuffix(layer),
      n: Object.values(counts).reduce((sum, v) => sum + v, 0),
      states: counts,
    }));

    return {
      stateChart: buildStacked(stateEntries, ["measured", "not_assayed", "measured_absent"]),
      availableChart: buildStacked(availableEntries, ["true", "false"]),
      rows,
    };
  }, [data]);

  return (
    <section className="explore-section">
      <h2>Per-layer coverage</h2>
      <p className="subtitle">Straight from build_manifest.json's coverage_state_counts, no recomputation.</p>
      {loading && <p>Loading…</p>}
      {error && <p className="uml-error">{error}</p>}
      {data && (
        <>
          {stateChart && (
            <div className="explore-chart-wrap explore-chart-scroll">
              <ChartCanvas config={stateChart} />
            </div>
          )}
          {availableChart && (
            <div className="explore-chart-wrap explore-chart-scroll">
              <ChartCanvas config={availableChart} />
            </div>
          )}
          <div className="explore-table-wrap">
            <table className="query-table">
              <thead>
                <tr>
                  <th>Layer</th>
                  <th>n (all states)</th>
                  <th>Breakdown</th>
                </tr>
              </thead>
              <tbody>
                {rows.map((row) => (
                  <tr key={row.layer}>
                    <td>{row.layer}</td>
                    <td>{row.n.toLocaleString()}</td>
                    <td>
                      {Object.entries(row.states)
                        .map(([state, count]) => state + ": " + count.toLocaleString())
                        .join(" · ")}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </>
      )}
    </section>
  );
}

function LineageCoverageSection() {
  const [data, setData] = useState<LineageCoverageResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    let cancelled = false;
    getEdaLineageCoverage()
      .then((res) => {
        if (!cancelled) setData(res);
      })
      .catch((err: Error) => {
        if (!cancelled) setError(err.message);
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, []);

  const matrix = useMemo(() => {
    if (!data) return null;
    const lineages = Array.from(new Set(data.rows.map((r) => r.lineage))).sort();
    const layers = Array.from(new Set(data.rows.map((r) => r.layer))).sort();
    const lookup = new Map(data.rows.map((r) => [r.lineage + "|" + r.layer, r]));
    const nByLineage = new Map<string, number>();
    for (const r of data.rows) nByLineage.set(r.lineage, r.n);
    return { lineages, layers, lookup, nByLineage };
  }, [data]);

  return (
    <section className="explore-section">
      <h2>Lineage × layer coverage bias</h2>
      <p className="subtitle">
        Measured fraction per (lineage, layer). A lineage below the minimum group size is flagged, not
        dropped.
      </p>
      {loading && <p>Loading…</p>}
      {error && <p className="uml-error">{error}</p>}
      {data && matrix && (
        <>
          <p className="explore-legend">
            <span className="explore-flag-swatch" /> flagged — fewer than {data.min_lineage_n} modelled
            lines in that lineage
          </p>
          <div className="explore-table-wrap explore-matrix-wrap">
            <table className="query-table explore-matrix">
              <thead>
                <tr>
                  <th>Lineage</th>
                  <th>n</th>
                  {matrix.layers.map((layer) => (
                    <th key={layer}>{stripSuffix(layer)}</th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {matrix.lineages.map((lineage) => (
                  <tr key={lineage}>
                    <td>{lineage}</td>
                    <td>{matrix.nByLineage.get(lineage)?.toLocaleString() ?? ""}</td>
                    {matrix.layers.map((layer) => {
                      const cell = matrix.lookup.get(lineage + "|" + layer);
                      if (!cell) return <td key={layer}>—</td>;
                      const alpha = Math.max(0.06, cell.measured_fraction);
                      return (
                        <td
                          key={layer}
                          className={cell.flagged ? "explore-matrix-cell explore-matrix-cell-flagged" : "explore-matrix-cell"}
                          style={{ background: "rgba(170, 59, 255, " + alpha + ")" }}
                          title={cell.n_measured + " of " + cell.n + " measured"}
                        >
                          {(cell.measured_fraction * 100).toFixed(0)}%{cell.flagged ? "*" : ""}
                        </td>
                      );
                    })}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </>
      )}
    </section>
  );
}

const PAIR_LABELS: Record<string, string> = {
  depmap_vs_hpa: "DepMap vs HPA",
  depmap_vs_geo: "DepMap vs GEO",
  depmap_rna_vs_protein: "RNA vs Protein",
};

function buildConcordanceChart(summaries: Record<string, ConcordanceSummary>) {
  const entries = Object.entries(summaries).filter(([, s]) => s.median_rho !== null);
  if (entries.length === 0) return null;
  const labels = entries.map(([key]) => PAIR_LABELS[key] ?? key);
  return {
    type: "bar" as const,
    data: {
      labels,
      datasets: [
        {
          type: "bar" as const,
          label: "IQR",
          data: entries.map(([, s]) => [s.iqr_low ?? s.median_rho, s.iqr_high ?? s.median_rho]),
          backgroundColor: "rgba(170, 59, 255, 0.35)",
        },
        {
          type: "scatter" as const,
          label: "median ρ",
          data: labels.map((label, i) => ({ x: label, y: entries[i][1].median_rho })),
          backgroundColor: "#aa3bff",
          pointRadius: 5,
        },
      ],
    },
    options: {
      responsive: true,
      scales: { y: { title: { display: true, text: "Spearman ρ" } } },
    },
  } as ChartConfiguration;
}

function ConcordanceGroup({ title, summaries }: { title: string; summaries: Record<string, ConcordanceSummary> }) {
  const chart = useMemo(() => buildConcordanceChart(summaries), [summaries]);
  return (
    <div className="explore-concordance-group">
      <h3>{title}</h3>
      {chart ? (
        <div className="explore-chart-wrap">
          <ChartCanvas config={chart} height={260} />
        </div>
      ) : (
        <p className="explore-placeholder">No pair in this group had an overlapping model set.</p>
      )}
      <div className="explore-table-wrap">
        <table className="query-table">
          <thead>
            <tr>
              <th>Pair</th>
              <th>Models scored</th>
              <th>Models in intersection</th>
              <th>Genes used (median)</th>
              <th>Genes floor</th>
            </tr>
          </thead>
          <tbody>
            {Object.entries(summaries).map(([key, s]) => (
              <tr key={key}>
                <td>{PAIR_LABELS[key] ?? key}</td>
                <td>{s.n_models_scored.toLocaleString()}</td>
                <td>{s.n_models_in_intersection.toLocaleString()}</td>
                <td>{s.median_n_genes_used ?? "—"}</td>
                <td>{s.min_genes_floor}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function ConcordanceSection() {
  const [data, setData] = useState<EdaConcordanceResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    let cancelled = false;
    getEdaConcordance()
      .then((res) => {
        if (!cancelled) setData(res);
      })
      .catch((err: Error) => {
        if (!cancelled) setError(err.message);
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, []);

  return (
    <section className="explore-section">
      <h2>RNA cross-source and RNA ↔ protein concordance</h2>
      <p className="subtitle">
        Precomputed by <code>backend/data_service/build_eda_aggregates.py</code>, never scanned at
        request time.
      </p>
      {loading && <p>Loading…</p>}
      {error && <p className="uml-error">{error}</p>}
      {data && !data.available && (
        <div className="diagnostic-banner">
          <p>{data.reason}</p>
          <p>
            Build it with: <code>{data.build_command}</code>
          </p>
        </div>
      )}
      {data && data.available && (
        <>
          <ConcordanceGroup title="RNA cross-source" summaries={data.rna_cross_source} />
          <ConcordanceGroup title="RNA ↔ protein" summaries={data.rna_protein} />
        </>
      )}
    </section>
  );
}

function EdaPage() {
  return (
    <div className="app-shell">
      <header>
        <h1>EDA</h1>
        <p className="subtitle">
          Per-layer coverage, lineage × layer coverage bias, and RNA concordance — recomputed from{" "}
          <code>data/processed/</code>, never a shipped image.
        </p>
      </header>
      <CoverageSection />
      <LineageCoverageSection />
      <ConcordanceSection />
    </div>
  );
}

export default EdaPage;
