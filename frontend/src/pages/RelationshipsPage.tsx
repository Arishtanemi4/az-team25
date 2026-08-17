import { useEffect, useMemo, useState } from "react";
import type { ChartConfiguration } from "chart.js/auto";
import { ChartCanvas } from "../components/charts/ChartCanvas";
import { FilterSelect } from "../components/FilterSelect";
import { GeneSearchSelect } from "../components/GeneSearchSelect";
import { listLineages } from "../api/geneServiceClient";
import {
  getRelationshipsCellLineDisease,
  getRelationshipsGeneCellLine,
  getRelationshipsGeneDisease,
} from "../api/dataServiceClient";
import type {
  CellLineDiseaseResponse,
  GeneCellLineResponse,
  GeneDiseaseResponse,
  RelationshipLayer,
} from "../types/data";
import type { GeneMatch } from "../types/scoring";
import "./explorePages.css";

const LAYER_OPTIONS: { value: RelationshipLayer; label: string }[] = [
  { value: "rna", label: "RNA (DepMap)" },
  { value: "rna_hpa", label: "RNA (HPA)" },
  { value: "rna_geo", label: "RNA (GEO)" },
  { value: "protein", label: "Protein" },
  { value: "dependency", label: "Dependency (CRISPR)" },
  { value: "copy_number", label: "Copy number" },
];

// A fixed, deterministic palette so lineage colors never depend on fetch order. Cycled with
// modulo if a query returns more distinct lineages than colors -- readability degrades
// gracefully rather than the app crashing on a long tail of rare lineages.
const PALETTE = [
  "#aa3bff", "#ff6b3b", "#2f9e44", "#1c7ed6", "#e64980", "#f59f00",
  "#12b886", "#7048e8", "#e03131", "#0ca678", "#495057", "#a04000",
];

// Small deterministic string hash -> [0, 1), used only to jitter scatter points horizontally
// within their lineage's category band for readability. Never used to compute or alter `value`.
function hashUnit(text: string): number {
  let h = 0;
  for (let i = 0; i < text.length; i++) h = (h * 31 + text.charCodeAt(i)) >>> 0;
  return (h % 1000) / 1000;
}

function GeneCellLineSection() {
  const [gene, setGene] = useState<GeneMatch | null>(null);
  const [layer, setLayer] = useState<RelationshipLayer>("rna");
  const [data, setData] = useState<GeneCellLineResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!gene) {
      setData(null);
      return;
    }
    let cancelled = false;
    setLoading(true);
    setError(null);
    getRelationshipsGeneCellLine(gene.ensembl_id, layer, 2000)
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
  }, [gene, layer]);

  const chart = useMemo(() => {
    if (!data || data.points.length === 0) return null;
    const categories = Array.from(
      new Set(data.points.map((p) => p.lineage ?? "unknown")),
    ).sort();
    const indexOf = new Map(categories.map((c, i) => [c, i]));
    const byLineage = new Map<string, { x: number; y: number; label: string }[]>();
    for (const p of data.points) {
      const lineage = p.lineage ?? "unknown";
      const center = indexOf.get(lineage) ?? 0;
      const jitter = (hashUnit(p.model_id) - 0.5) * 0.7;
      const list = byLineage.get(lineage) ?? [];
      list.push({ x: center + jitter, y: p.value, label: p.cell_line_name ?? p.model_id });
      byLineage.set(lineage, list);
    }
    const datasets = Array.from(byLineage.entries()).map(([lineage, points], i) => ({
      label: lineage,
      data: points,
      backgroundColor: PALETTE[i % PALETTE.length],
      pointRadius: 3.5,
    }));
    return {
      type: "scatter" as const,
      data: { datasets },
      options: {
        responsive: true,
        scales: {
          x: {
            min: -0.5,
            max: categories.length - 0.5,
            ticks: {
              stepSize: 1,
              callback: (val: number | string) => categories[Number(val)] ?? "",
            },
            title: { display: true, text: "lineage" },
          },
          y: { title: { display: true, text: data.unit } },
        },
        plugins: {
          tooltip: {
            callbacks: {
              // Chart.js's TooltipItem<T> is a large generic union; a plain `any` here is the
              // pragmatic choice over fighting the tooltip callback's own strict typing.
              label: (ctx: any) => (ctx.raw.label ?? "") + ": " + ctx.raw.y,
            },
          },
        },
      },
    } as ChartConfiguration;
  }, [data]);

  return (
    <section className="explore-section">
      <h2>Gene ↔ cell line</h2>
      <p className="subtitle">One gene's measured value across cell lines, coloured by lineage.</p>
      <div className="explore-controls">
        <GeneSearchSelect label="Gene" selected={gene} onChange={setGene} />
        <label className="query-field">
          Layer
          <select value={layer} onChange={(e) => setLayer(e.target.value as RelationshipLayer)}>
            {LAYER_OPTIONS.map((opt) => (
              <option key={opt.value} value={opt.value}>
                {opt.label}
              </option>
            ))}
          </select>
        </label>
      </div>

      {!gene && <p className="explore-placeholder">Pick a gene to plot its values across cell lines.</p>}
      {loading && <p>Loading…</p>}
      {error && <p className="uml-error">{error}</p>}
      {data && (
        <>
          <p className="query-denominator">
            Showing <strong>{data.points.length.toLocaleString()}</strong> of{" "}
            <strong>{data.n_models_total.toLocaleString()}</strong> cell lines — layer unit:{" "}
            {data.unit}
          </p>
          {chart ? (
            <div className="explore-chart-wrap">
              <ChartCanvas config={chart} />
            </div>
          ) : (
            <p className="explore-placeholder">No measured values for this gene on this layer.</p>
          )}
        </>
      )}
    </section>
  );
}

function CellLineDiseaseSection() {
  const [lineages, setLineages] = useState<string[]>([]);
  const [lineage, setLineage] = useState<string | null>(null);
  const [data, setData] = useState<CellLineDiseaseResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    listLineages()
      .then(setLineages)
      .catch(() => setLineages([]));
  }, []);

  useEffect(() => {
    let cancelled = false;
    setLoading(true);
    setError(null);
    getRelationshipsCellLineDisease(lineage)
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
  }, [lineage]);

  const { chart, shownTotal } = useMemo(() => {
    if (!data || data.cells.length === 0) return { chart: null, shownTotal: 0 };
    const byDisease = new Map<string, number>();
    for (const cell of data.cells) {
      byDisease.set(cell.primary_disease, (byDisease.get(cell.primary_disease) ?? 0) + cell.n_models);
    }
    const entries = Array.from(byDisease.entries()).sort((a, b) => b[1] - a[1]);
    const total = entries.reduce((sum, [, n]) => sum + n, 0);
    return {
      shownTotal: total,
      chart: {
        type: "bar" as const,
        data: {
          labels: entries.map(([disease]) => disease),
          datasets: [{ label: "cell lines", data: entries.map(([, n]) => n), backgroundColor: PALETTE[0] }],
        },
        options: {
          responsive: true,
          scales: { y: { title: { display: true, text: "n cell lines" } } },
          plugins: { legend: { display: false } },
        },
      } as ChartConfiguration,
    };
  }, [data]);

  return (
    <section className="explore-section">
      <h2>Cell line ↔ disease</h2>
      <p className="subtitle">How many modelled cell lines exist per disease, optionally filtered by lineage.</p>
      <div className="explore-controls">
        <FilterSelect label="Lineage" options={lineages} value={lineage} onChange={setLineage} />
      </div>
      {loading && <p>Loading…</p>}
      {error && <p className="uml-error">{error}</p>}
      {data && (
        <>
          <p className="query-denominator">
            Showing <strong>{shownTotal.toLocaleString()}</strong> of{" "}
            <strong>{data.n_models_total.toLocaleString()}</strong> cell lines
          </p>
          {chart && (
            <div className="explore-chart-wrap explore-chart-scroll">
              <ChartCanvas config={chart} />
            </div>
          )}
        </>
      )}
    </section>
  );
}

function GeneDiseaseSection() {
  const [gene, setGene] = useState<GeneMatch | null>(null);
  const [layer, setLayer] = useState<RelationshipLayer>("rna");
  const [minN, setMinN] = useState(5);
  const [data, setData] = useState<GeneDiseaseResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!gene) {
      setData(null);
      return;
    }
    let cancelled = false;
    setLoading(true);
    setError(null);
    getRelationshipsGeneDisease(gene.ensembl_id, layer, minN)
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
  }, [gene, layer, minN]);

  const chart = useMemo(() => {
    if (!data) return null;
    const plottable = data.groups.filter((g) => g.median !== null);
    if (plottable.length === 0) return null;
    return {
      type: "bar" as const,
      data: {
        labels: plottable.map((g) => g.primary_disease),
        datasets: [
          {
            label: "median",
            data: plottable.map((g) => g.median),
            backgroundColor: plottable.map((g) => (g.flagged ? "#f59f0080" : PALETTE[0])),
          },
        ],
      },
      options: {
        responsive: true,
        scales: { y: { title: { display: true, text: "median (" + data.layer + ")" } } },
        plugins: { legend: { display: false } },
      },
    } as ChartConfiguration;
  }, [data]);

  return (
    <section className="explore-section">
      <h2>Gene ↔ disease</h2>
      <p className="subtitle">
        Median measured value per disease. A disease with fewer than the minimum measured lines is
        flagged, never dropped.
      </p>
      <div className="explore-controls">
        <GeneSearchSelect label="Gene" selected={gene} onChange={setGene} />
        <label className="query-field">
          Layer
          <select value={layer} onChange={(e) => setLayer(e.target.value as RelationshipLayer)}>
            {LAYER_OPTIONS.map((opt) => (
              <option key={opt.value} value={opt.value}>
                {opt.label}
              </option>
            ))}
          </select>
        </label>
        <label className="query-field">
          Minimum measured lines
          <input
            type="number"
            min={1}
            value={minN}
            onChange={(e) => setMinN(Number.isFinite(e.target.valueAsNumber) ? Math.max(1, Math.floor(e.target.valueAsNumber)) : 1)}
          />
        </label>
      </div>

      {!gene && <p className="explore-placeholder">Pick a gene to compare it across diseases.</p>}
      {loading && <p>Loading…</p>}
      {error && <p className="uml-error">{error}</p>}
      {data && (
        <>
          <p className="explore-legend">
            <span className="explore-flag-swatch" /> flagged — fewer than {data.min_n} measured lines
          </p>
          {chart ? (
            <div className="explore-chart-wrap explore-chart-scroll">
              <ChartCanvas config={chart} />
            </div>
          ) : (
            <p className="explore-placeholder">No disease group has a measured value on this layer.</p>
          )}
          <div className="explore-table-wrap">
            <table className="query-table">
              <thead>
                <tr>
                  <th>Disease</th>
                  <th>Measured</th>
                  <th>Total lines</th>
                  <th>Median</th>
                  <th>IQR</th>
                  <th>Flagged</th>
                </tr>
              </thead>
              <tbody>
                {data.groups.map((g) => (
                  <tr key={g.primary_disease} className={g.flagged ? "explore-row-flagged" : ""}>
                    <td>{g.primary_disease}</td>
                    <td>{g.n_measured}</td>
                    <td>{g.n_models}</td>
                    <td>{g.median === null ? "—" : g.median.toFixed(3)}</td>
                    <td>{g.iqr_low === null || g.iqr_high === null ? "—" : g.iqr_low.toFixed(3) + "–" + g.iqr_high.toFixed(3)}</td>
                    <td>{g.flagged ? "yes" : "no"}</td>
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

function RelationshipsPage() {
  return (
    <div className="app-shell">
      <header>
        <h1>Relationships</h1>
        <p className="subtitle">
          Gene ↔ cell line, cell line ↔ disease, and gene ↔ disease — every panel states how many
          lines were measured against how many exist, and a small group is flagged, never dropped.
        </p>
      </header>
      <GeneCellLineSection />
      <CellLineDiseaseSection />
      <GeneDiseaseSection />
    </div>
  );
}

export default RelationshipsPage;
