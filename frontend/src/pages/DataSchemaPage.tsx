import { useEffect, useState, type ReactNode } from "react";
import { getSchema } from "../api/dataServiceClient";
import type { ColumnInfo, SchemaResponse, TableInfo } from "../types/data";
import "./dataPages.css";

// Fixed hand-authored layout (PRODUCT_SURFACE.md SS1 decision 4): `cell_lines` and
// `gene_reference` are the two hubs every other table joins to via ModelID / ensembl_id
// (table_registry.py). Positions are hardcoded; every box's *content* (grain, row count,
// columns) is read live from GET /data/schema, so the diagram can never silently drift from
// the backend registry the way a static image would.
const BOX_W = 210;
const BOX_H = 210;
const VIEW_W = 1660;
const VIEW_H = 950;

const LAYOUT: Record<string, { x: number; y: number }> = {
  coverage: { x: 70, y: 40 },
  genome_signatures: { x: 70, y: 260 },
  mirna: { x: 70, y: 480 },
  metabolomics: { x: 70, y: 700 },
  cell_lines: { x: 400, y: 370 },
  expression_rna: { x: 740, y: 10 },
  expression_rna_hpa: { x: 740, y: 230 },
  expression_rna_geo: { x: 740, y: 450 },
  protein: { x: 740, y: 670 },
  dependency: { x: 1070, y: 10 },
  copy_number: { x: 1070, y: 230 },
  mutations: { x: 1070, y: 450 },
  fusions: { x: 1070, y: 670 },
  gene_reference: { x: 1400, y: 340 },
};

function boxCenter(name: string): { x: number; y: number } | null {
  const pos = LAYOUT[name];
  if (!pos) return null;
  return { x: pos.x + BOX_W / 2, y: pos.y + BOX_H / 2 };
}

// ModelID and ensembl_id are the only two automatic join keys (PROJECT_ARCHITECTURE.md SS4).
// Coloring both the column name and the edge it participates in by the same key name is what
// makes a join key readable without the legend.
function keyColor(columnName: string): string | undefined {
  if (columnName === "ModelID") return "var(--accent)";
  if (columnName === "ensembl_id") return "#a04000";
  return undefined;
}

const ROLE_ORDER: Record<string, number> = { join_key: 0, value: 1, state: 2 };
const MAX_COLUMNS_SHOWN = 5;

function criticalColumns(table: TableInfo): ColumnInfo[] {
  return table.columns
    .filter((c) => c.role !== "context")
    .sort((a, b) => (ROLE_ORDER[a.role] ?? 3) - (ROLE_ORDER[b.role] ?? 3));
}

function TableBox({ table }: { table: TableInfo }) {
  const pos = LAYOUT[table.name];
  if (!pos) return null;
  const critical = criticalColumns(table);
  const contextCount = table.columns.length - critical.length;
  const shown = critical.slice(0, MAX_COLUMNS_SHOWN);
  const hiddenCritical = critical.length - shown.length;

  return (
    <foreignObject x={pos.x} y={pos.y} width={BOX_W} height={BOX_H}>
      <div className="uml-box">
        <div className="uml-box-title">{table.name}</div>
        <div className="uml-box-grain">{table.grain}</div>
        <div className="uml-box-rows">{table.row_count.toLocaleString()} rows</div>
        <ul className="uml-col-list">
          {shown.map((col) => {
            const color = keyColor(col.name);
            return (
              <li
                key={col.name}
                className={color ? "uml-col uml-col-key" : `uml-col uml-col-${col.role}`}
                style={color ? { color } : undefined}
              >
                {col.name}
                {col.role === "value" && <span className="uml-col-dtype"> · {col.dtype}</span>}
              </li>
            );
          })}
        </ul>
        {hiddenCritical > 0 && (
          <div className="uml-col-more">+{hiddenCritical} more critical column{hiddenCritical === 1 ? "" : "s"}</div>
        )}
        {contextCount > 0 && (
          <div className="uml-col-more">{contextCount} context column{contextCount === 1 ? "" : "s"} not shown</div>
        )}
      </div>
    </foreignObject>
  );
}

function JoinEdges({ tables }: { tables: TableInfo[] }) {
  const edges: ReactNode[] = [];
  for (const table of tables) {
    const from = boxCenter(table.name);
    if (!from) continue;
    for (const join of table.joins) {
      const to = boxCenter(join.to);
      if (!to) continue;
      const color = keyColor(join.on) ?? "var(--border)";
      edges.push(
        <line
          key={`${table.name}-${join.to}-${join.on}`}
          x1={from.x}
          y1={from.y}
          x2={to.x}
          y2={to.y}
          stroke={color}
          strokeWidth={2}
          opacity={0.55}
        />,
      );
    }
  }
  return <g>{edges}</g>;
}

function DataSchemaPage() {
  const [schema, setSchema] = useState<SchemaResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    let cancelled = false;
    getSchema()
      .then((data) => {
        if (!cancelled) setSchema(data);
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
    <div className="app-shell">
      <header>
        <h1>Data schema</h1>
        <p className="subtitle">
          The 14 `data/processed/` tables, critical columns only. `ModelID` and `ensembl_id` are the
          only two automatic join keys.
        </p>
      </header>

      {loading && <p>Loading schema…</p>}
      {error && <p className="uml-error">Could not load the data schema: {error}</p>}

      {schema && (
        <>
          <div className="uml-legend">
            <span className="uml-legend-item">
              <span className="uml-legend-swatch" style={{ background: "var(--accent)" }} />
              ModelID join
            </span>
            <span className="uml-legend-item">
              <span className="uml-legend-swatch" style={{ background: "#a04000" }} />
              ensembl_id join
            </span>
            <span className="uml-legend-note">
              Only these two columns are ever automatic join keys — every other identifier
              (names, CCLE labels, RRIDs) is a cross-check, never a join key.
            </span>
          </div>
          <div className="uml-svg-wrap">
            <svg viewBox={`0 0 ${VIEW_W} ${VIEW_H}`} width="100%" role="img" aria-label="Data schema diagram">
              <JoinEdges tables={schema.tables} />
              {schema.tables.map((table) => (
                <TableBox key={table.name} table={table} />
              ))}
            </svg>
          </div>
        </>
      )}
    </div>
  );
}

export default DataSchemaPage;
