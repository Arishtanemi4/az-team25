import { useEffect, useState, type ReactNode } from "react";
import { getValidationStudies, getValidationStudy } from "../api/dataServiceClient";
import type { ValidationStudySummary, ValidationStudyResponse } from "../types/data";
import "./validationPage.css";

function ComparableTag({ value }: { value: string }) {
  const slug = value.toLowerCase().startsWith("partial") ? "partial" : value.toLowerCase() === "yes" ? "yes" : "no";
  return <span className={`validation-tag validation-tag--${slug}`}>{value}</span>;
}

function ManifestTable({ studies }: { studies: ValidationStudySummary[] }) {
  return (
    <div className="query-table-wrap">
      <table className="query-table validation-manifest-table">
        <thead>
          <tr>
            <th>Study</th>
            <th>Mode</th>
            <th>Comparable to our ranking?</th>
            <th>Reason</th>
          </tr>
        </thead>
        <tbody>
          {studies.map((s) => (
            <tr key={s.study}>
              <td>{s.study}</td>
              <td>{s.mode}</td>
              <td>
                <ComparableTag value={s.comparable_to_ranking} />
              </td>
              <td className="validation-reason-cell">{s.reason}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function StudyTable({ study }: { study: ValidationStudyResponse }) {
  return (
    <>
      <p className="query-denominator">
        Showing <strong>{study.returned_rows.toLocaleString()}</strong> of{" "}
        <strong>{study.total_rows.toLocaleString()}</strong> rows
      </p>
      <div className="query-table-wrap">
        <table className="query-table">
          <thead>
            <tr>
              {study.columns.map((col) => (
                <th key={col}>{col}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {study.rows.map((row, i) => (
              <tr key={i}>
                {study.columns.map((col) => (
                  <td key={col}>{row[col] === null || row[col] === undefined ? "" : String(row[col])}</td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </>
  );
}

// Fetches and renders one allow-listed study table (getValidationStudy's own row/column cap).
function StudyTableSection({ studyKey, label }: { studyKey: string; label: string }) {
  const [data, setData] = useState<ValidationStudyResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    let cancelled = false;
    getValidationStudy(studyKey)
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
  }, [studyKey]);

  return (
    <div className="validation-subtable">
      <h3>{label}</h3>
      {loading && <p>Loading…</p>}
      {error && <p className="uml-error">Could not load this table: {error}</p>}
      {data && <StudyTable study={data} />}
    </div>
  );
}

// One section per manifest row -- each of its tables fetches independently so a slow/failing
// study (e.g. jin2023's 500-row preview) never blocks the others from rendering.
function StudySection({
  tables,
  title,
  verdict,
  children,
}: {
  tables: { key: string; label: string }[];
  title: string;
  verdict: string;
  children: ReactNode;
}) {
  return (
    <section className="validation-section">
      <h2>{title}</h2>
      {children}
      <p className="validation-verdict">{verdict}</p>
      {tables.map((t) => (
        <StudyTableSection key={t.key} studyKey={t.key} label={t.label} />
      ))}
    </section>
  );
}

function ValidationPage() {
  const [studies, setStudies] = useState<ValidationStudySummary[] | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;
    getValidationStudies()
      .then((res) => {
        if (!cancelled) setStudies(res.studies);
      })
      .catch((err: Error) => {
        if (!cancelled) setError(err.message);
      });
    return () => {
      cancelled = true;
    };
  }, []);

  return (
    <div className="app-shell">
      <header>
        <h1>Validation</h1>
        <p className="subtitle">
          Four published external studies, read as corroboration only — never used to fit or calibrate the
          scorer. This page reproduces <code>validation/validation_report.ipynb</code>'s own recorded
          verdicts; it never recomputes them, and it shows the negative verdicts as prominently as the
          positive one.
        </p>
      </header>

      {error && <p className="uml-error">Could not load the study manifest: {error}</p>}

      {studies && (
        <section className="validation-section">
          <h2>Overview — what's comparable, at a glance</h2>
          <ManifestTable studies={studies} />
        </section>
      )}

      <StudySection
        tables={[{ key: "jin2023", label: "Per-(cohort, cell line) comparison" }]}
        title="Jin et al. (2023) — whole-transcriptome tumour representativeness"
        verdict="Verdict: PARTIAL / descriptive corroboration only. Jin et al.'s top-ranked LIHC lines score High confidence and land in the upper half of our own ALB/liver ranking — encouraging, but a single gene-grounded case out of 26 cohorts, and a moderate correlation is not proof either method is right."
      >
        <p>
          Jin et al. ranked ~1,000 cell lines against 26 TCGA cancer cohorts by whole-transcriptome
          correlation and GSEA — a different method from this project's per-gene desirability scorer. Their
          26 worksheets are keyed by TCGA cohort, not by gene.
        </p>
        <p className="validation-note">
          <strong>26 of 26 cohorts abandoned in <code>lineage_only</code> mode</strong> — a real bug, not
          fixed here: <code>score.score_panel([], [], ...)</code> (no genes at all) raises inside{" "}
          <code>scoring/score.py::_load_filtered</code>'s Parquet predicate-pushdown path (
          <code>pyarrow</code> cannot type an empty <code>ensembl_id</code> filter). The 1,249{" "}
          <code>abandoned</code> rows below are shown, not hidden.
        </p>
        <p className="validation-note">
          The one gene-grounded case: LIHC / albumin (<code>ALB</code>), Jin et al.'s own worked example.
          24 of 25 Jin-ranked LIHC lines were scored by our ALB/liver query. Spearman rho = 0.69 (p = 0.000,
          n = 24) against Jin's overall rank — a rough summary over a small sample, not a validation metric;
          the two methods measure different things.
        </p>
      </StudySection>

      <StudySection
        tables={[
          { key: "tumorcomparer", label: "Tier 1 — our own algorithm, RTK-RAS/WNT panel, SKCM/LIHC" },
          { key: "tumorcomparer_tier2", label: "Tier 2 — TumorComparer's own precomputed .rds files, per-file status" },
        ]}
        title="TumorComparer (Sinha et al. 2021) — the one new gene-level dataset"
        verdict="Verdict: reported honestly, not forced. Tier 2's numeric cross-check against TumorComparer's own published numbers could not be attempted at all (data-format blocker, not a unit mismatch)."
      >
        <p>
          TumorComparer's Figure 4 restricts SKCM and LIHC cell lines to two pathway gene panels — RTK-RAS
          and WNT — claiming SKCM scores higher on RTK-RAS and LIHC higher on WNT. The exact gene lists used
          in Figure 4 are never published; the 152-gene panel used here (Sanchez-Vega et al. 2018, pulled
          from TumorComparer's own GitHub repo) is a legitimately-cited proxy, not a certified match.
        </p>
        <p className="validation-note">
          <strong>Tier 1 (our own algorithm):</strong> SKCM — RTK-RAS mean D=0.33 vs. WNT mean D=0.29 —
          matches the paper's claim. LIHC — WNT mean D=0.17 vs. RTK-RAS mean D=0.25 — does <strong>NOT</strong>{" "}
          match the paper's claim.
        </p>
        <p className="validation-note">
          <strong>Tier 2 (TumorComparer's own precomputed output):</strong> blocked for all 6 precomputed
          files — two parse cleanly but hold only genome-wide, not pathway-restricted, scores; the other
          four parse into zero extractable objects (nested R lists <code>pyreadr</code> cannot flatten
          without R itself).
        </p>
      </StudySection>

      <StudySection
        tables={[{ key: "celligner", label: "Coverage / lineage-agreement, per cell line" }]}
        title="Celligner (Warren et al. 2021) — not comparable as a ranking validation"
        verdict="Verdict: NOT COMPARABLE as a ranking validation — coverage check only. Celligner's downloaded file has no gene axis, so it cannot validate a gene-in/gene-out ranking algorithm."
      >
        <p>
          Celligner's downloaded file carries only per-sample lineage/subtype/cluster metadata — the paper's
          gene-level aligned-expression matrix was deliberately not downloaded. With no gene axis,{" "}
          <code>score.score_panel</code> is never called here; this is a data-coverage check instead.
        </p>
        <p className="validation-note">
          1243/1249 Celligner cell-line rows resolve to a <code>ModelID</code> we have; of those, 1228/1243
          agree on lineage label. Most disagreements are <code>engineered_*</code> lineage labels resolving
          to the underlying tissue (e.g. <code>engineered_kidney</code> vs. <code>kidney</code>) — a real
          annotation difference, not noise.
        </p>
      </StudySection>

      <StudySection
        tables={[{ key: "netcellmatch", label: "Coverage / lineage-agreement, per cell line" }]}
        title="NetCellMatch (Desai et al. 2022) — not comparable as a ranking validation"
        verdict="Verdict: NOT COMPARABLE as a ranking validation — coverage check only. Protein/phospho epitope level, not gene level, and only 1 of the paper's 3 cancer cohorts is even on disk."
      >
        <p>
          NetCellMatch matches cell lines to patients using a 233-antibody RPPA (protein/phospho) panel —
          not gene expression, and breast-cohort only on disk. Same treatment as Celligner: no{" "}
          <code>score.score_panel</code> call, a coverage-only check instead.
        </p>
        <p className="validation-note">
          52/58 NetCellMatch breast-cohort cell-line rows resolve to a name we have; of those, 52/52 agree
          on lineage (expected — the file is breast-only).
        </p>
      </StudySection>

      <section className="validation-section">
        <h2>Summary</h2>
        <p>
          None of these four studies feeds back into <code>scoring/</code>. All four stay external
          corroboration until a separate label and evaluation design is approved.
        </p>
      </section>
    </div>
  );
}

export default ValidationPage;
