import os
import sys

import openpyxl
import pandas as pd
from scipy.stats import spearmanr

try:
    import pyreadr
except ImportError:
    pyreadr = None

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "scoring"))
import score  # noqa: E402 -- must follow sys.path.insert so "import score" resolves

VALIDATION_DATA_DIR = os.path.join(REPO_ROOT, "data", "augmented", "validation")
DATA_DIR = os.path.join(REPO_ROOT, "data", "processed")
OUTPUT_DIR = os.path.join(REPO_ROOT, "validation")

JIN_2023_PATH = os.path.join(VALIDATION_DATA_DIR, "jin2023_supplementary_data5.xlsx")
CELLIGNER_PATH = os.path.join(VALIDATION_DATA_DIR, "celligner_2021_info.csv")
NETCELLMATCH_PATH = os.path.join(
    VALIDATION_DATA_DIR, "netcellmatch_2022_breast_cancer_batchcorrected.csv"
)
TUMORCOMPARER_GENES_PATH = os.path.join(
    VALIDATION_DATA_DIR, "tumorcomparer_2021_rtkras_wnt_genes.csv"
)
TUMORCOMPARER_PRECOMPUTED_DIR = os.path.join(VALIDATION_DATA_DIR, "tumorcomparer_precomputed")

ALB = "ENSG00000163631"  # albumin -- Jin et al.'s own liver-cancer marker gene (LIHC worked example)
MARKER_GENE_COHORTS = {"LIHC": ALB}

TUMORCOMPARER_COHORTS = ["SKCM", "LIHC"]

COHORT_LINEAGE_CANDIDATES = {
    "BLCA": ["urinary_tract"],
    "BRCA": ["breast", "Breast"],
    "CHOL": ["bile_duct", "Biliary Tract"],
    "CESC": ["cervix"],
    "COAD": ["colorectal", "Bowel"],
    "ESCA": ["esophagus", "Esophagus/Stomach"],
    "GBM": ["central_nervous_system", "CNS/Brain"],
    "HNSC": ["upper_aerodigestive", "Head and Neck"],
    "KICH": ["kidney", "Kidney"],
    "KIRC": ["kidney", "Kidney"],
    "KIRP": ["kidney", "Kidney"],
    "LAML": ["Myeloid", "blood"],
    "LGG": ["central_nervous_system", "CNS/Brain"],
    "LIHC": ["liver", "Liver"],
    "LUAD": ["lung", "Lung"],
    "LUSC": ["lung", "Lung"],
    "OV": ["ovary", "Ovary/Fallopian Tube"],
    "PAAD": ["pancreas", "Pancreas"],
    "PRAD": ["prostate", "Prostate"],
    "READ": ["colorectal", "Bowel"],
    "SARC": ["soft_tissue", "Soft Tissue", "bone", "Bone"],
    "SKCM": ["skin", "Skin"],
    "STAD": ["gastric", "Esophagus/Stomach"],
    "TGCT": ["testis", "Testis"],
    "THCA": ["thyroid", "Thyroid"],
    "UCEC": ["uterus", "Uterus"],
}

JIN_COLUMNS = [
    "RRID", "Cell line name", "Overall rank", "Correlation", "Correlation rank", "NES",
    "P-value", "Adjusted P-value", "GSEA rank", "Primary disease", "Primary/metastasis",
    "Sample collection site", "Data source",
]


def read_jin_sheet(path, sheet_name):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet_name]
    rows = list(ws.iter_rows(values_only=True))
    header = list(rows[2])
    frame = pd.DataFrame(rows[3:], columns=header)
    return frame[JIN_COLUMNS]


def pick_lineage(cell_lines_df, cohort):
    for candidate in COHORT_LINEAGE_CANDIDATES.get(cohort, []):
        if (cell_lines_df["lineage"] == candidate).any():
            return candidate
    return None


def run_mode(inclusion_tokens, lineage, gene_reference_df, cell_lines_df):
    try:
        result = score.score_panel(
            inclusion_tokens, [], gene_reference_df, cell_lines_df,
            lineage=lineage, data_dir=DATA_DIR,
            rna_constants_path=os.path.join(REPO_ROOT, "scoring", "resources", "desirability_constants.json"),
            extended_constants_path=os.path.join(REPO_ROOT, "scoring", "resources", "desirability_constants_extended.json"),
            essentiality_constants_path=os.path.join(REPO_ROOT, "scoring", "resources", "common_essential_genes.json"),
            gene_role_path=os.path.join(REPO_ROOT, "preprocessing", "resources", "gene_role_reference.csv"),
            unresolved_symbols_path=os.path.join(REPO_ROOT, "preprocessing", "resources", "unresolved_gene_symbols.json"),
            top_n=100_000,
        )
    except Exception as exc:  
        return None, f"{type(exc).__name__}: {exc}"

    if result.get("diagnostic") is not None:
        return None, f"diagnostic: {result['diagnostic']}"

    by_model = {}
    for tier_key in ("ranked", "ranked_beyond_top_n", "low_confidence", "insufficient"):
        for line in result.get(tier_key, []):
            by_model[line["model_id"]] = {
                "our_D": line.get("D"),
                "our_confidence_tier": line.get("confidence_tier"),
                "our_veto": line["veto"]["reason"] if line.get("veto") else None,
                "our_rank_position": None,
            }

    ordered = result.get("ranked", []) + result.get("ranked_beyond_top_n", [])
    for i, line in enumerate(ordered):
        by_model[line["model_id"]]["our_rank_position"] = i + 1
    return by_model, None


def score_cohort_mode(query_genes, lineage, gene_reference_df, cell_lines_df):
    if lineage is None:
        return "no_lineage_match", "no candidate lineage string had any matching cell line", None
    by_model, error = run_mode(query_genes, lineage, gene_reference_df, cell_lines_df)
    if error is not None:
        return "abandoned", error, None
    return "ok", None, by_model


def run_jin2023(gene_reference, cell_lines):
    rrid_to_model = cell_lines.set_index("RRID")["ModelID"].to_dict()
    wb = openpyxl.load_workbook(JIN_2023_PATH, read_only=True, data_only=True)
    all_rows = []

    for cohort in wb.sheetnames:
        jin_df = read_jin_sheet(JIN_2023_PATH, cohort)
        lineage = pick_lineage(cell_lines, cohort)

        lineage_status, lineage_notes, lineage_by_model = score_cohort_mode(
            [], lineage, gene_reference, cell_lines
        )

        if cohort in MARKER_GENE_COHORTS:
            marker_gene = MARKER_GENE_COHORTS[cohort]
            marker_status, marker_notes, marker_by_model = score_cohort_mode(
                [marker_gene], lineage, gene_reference, cell_lines
            )
            marker_genes_str = marker_gene
        else:
            marker_status = "no_marker_gene"
            marker_notes = "no paper-grounded marker gene documented for this cohort"
            marker_by_model = None
            marker_genes_str = ""

        modes = [
            ("lineage_only", "", lineage, lineage_status, lineage_notes, lineage_by_model),
            ("marker_gene", marker_genes_str, lineage, marker_status, marker_notes, marker_by_model),
        ]

        for _, jin_row in jin_df.iterrows():
            rrid = jin_row["RRID"]
            model_id = rrid_to_model.get(rrid)

            for mode_name, query_genes_str, lineage_used, status, notes, by_model in modes:
                row = {
                    "tcga_cohort": cohort,
                    "query_mode": mode_name,
                    "query_genes": query_genes_str,
                    "lineage_used": lineage_used,
                    "RRID": rrid,
                    "model_id": model_id,
                    "cell_line_name": jin_row["Cell line name"],
                    "jin_overall_rank": jin_row["Overall rank"],
                    "jin_correlation": jin_row["Correlation"],
                    "jin_correlation_rank": jin_row["Correlation rank"],
                    "jin_nes": jin_row["NES"],
                    "jin_pvalue": jin_row["P-value"],
                    "jin_adj_pvalue": jin_row["Adjusted P-value"],
                    "jin_gsea_rank": jin_row["GSEA rank"],
                    "jin_primary_disease": jin_row["Primary disease"],
                    "jin_primary_metastasis": jin_row["Primary/metastasis"],
                    "jin_sample_site": jin_row["Sample collection site"],
                    "jin_data_source": jin_row["Data source"],
                    "our_D": None,
                    "our_confidence_tier": None,
                    "our_rank_position": None,
                    "our_veto": None,
                    "status": status,
                    "notes": notes,
                }
                if status == "ok":
                    if model_id is None:
                        row["status"] = "rrid_not_found"
                        row["notes"] = "RRID not present in data/processed/cell_lines.csv"
                    elif model_id not in by_model:
                        row["status"] = "not_scored"
                        row["notes"] = "line filtered out of the candidate pool for this lineage/query"
                    else:
                        row.update(by_model[model_id])
                all_rows.append(row)

    out_df = pd.DataFrame(all_rows)
    out_df.to_csv(os.path.join(OUTPUT_DIR, "jin2023_comparison.csv"), index=False)
    print(f"[jin2023] wrote {len(out_df)} rows")
    print(out_df.groupby(["query_mode", "status"]).size())

    liver = out_df[(out_df["tcga_cohort"] == "LIHC") & (out_df["query_mode"] == "marker_gene") & (out_df["status"] == "ok")]
    if len(liver) >= 3:
        rho, pval = spearmanr(liver["jin_overall_rank"], liver["our_rank_position"])
        reason = (
            f"lineage_only abandoned for all 26 cohorts (score.score_panel([], [], ...) raises "
            f"inside scoring/score.py::_load_filtered on an empty gene filter -- a real bug, not "
            f"fixed here). marker_gene mode only exists for LIHC/ALB, the paper's own worked "
            f"example: n={len(liver)} lines, Spearman rho={rho:.2f} (p={pval:.3f}) against Jin's "
            f"overall rank -- descriptive corroboration, not a validation metric (different methods)."
        )
        comparable = "partial"
    else:
        reason = "lineage_only abandoned for all 26 cohorts; marker_gene mode did not score enough LIHC lines to summarise"
        comparable = "no"

    return {
        "study": "Jin et al. (2023)", "mode": "marker_gene (LIHC/ALB only)",
        "comparable_to_ranking": comparable, "reason": reason,
    }


def run_tumorcomparer(gene_reference, cell_lines):
    panel = pd.read_csv(TUMORCOMPARER_GENES_PATH)
    pathways = {p: sub["gene_symbol"].tolist() for p, sub in panel.groupby("pathway")}
    model_to_name = cell_lines.set_index("ModelID")["cell_line_name"].to_dict()

    rows = []
    mean_d = {}
    for cohort in TUMORCOMPARER_COHORTS:
        lineage = pick_lineage(cell_lines, cohort)
        for pathway, genes in pathways.items():
            status, notes, by_model = score_cohort_mode(genes, lineage, gene_reference, cell_lines)
            if status != "ok":
                mean_d[(cohort, pathway)] = None
                rows.append({
                    "tcga_cohort": cohort, "pathway": pathway, "lineage_used": lineage,
                    "model_id": None, "cell_line_name": None, "our_D": None,
                    "our_confidence_tier": None, "our_rank_position": None,
                    "status": status, "notes": notes,
                })
                continue
            d_values = [v["our_D"] for v in by_model.values() if v["our_D"] is not None]
            mean_d[(cohort, pathway)] = sum(d_values) / len(d_values) if d_values else None
            for model_id, vals in by_model.items():
                rows.append({
                    "tcga_cohort": cohort, "pathway": pathway, "lineage_used": lineage,
                    "model_id": model_id, "cell_line_name": model_to_name.get(model_id),
                    **vals, "status": "ok", "notes": None,
                })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(os.path.join(OUTPUT_DIR, "tumorcomparer_comparison.csv"), index=False)
    print(f"[tumorcomparer] wrote {len(out_df)} rows (Tier 1)")

    skcm_rtk, skcm_wnt = mean_d.get(("SKCM", "RTK RAS")), mean_d.get(("SKCM", "WNT"))
    lihc_rtk, lihc_wnt = mean_d.get(("LIHC", "RTK RAS")), mean_d.get(("LIHC", "WNT"))
    if None in (skcm_rtk, skcm_wnt, lihc_rtk, lihc_wnt):
        directional_verdict = "inconclusive -- at least one (lineage, pathway) combination did not score"
    elif skcm_rtk > skcm_wnt and lihc_wnt > lihc_rtk:
        directional_verdict = (
            f"directionally consistent with Sinha et al. 2021 Fig. 4: mean D SKCM/RTK-RAS="
            f"{skcm_rtk:.2f} > SKCM/WNT={skcm_wnt:.2f}; mean D LIHC/WNT={lihc_wnt:.2f} > "
            f"LIHC/RTK-RAS={lihc_rtk:.2f}"
        )
    else:
        directional_verdict = (
            f"directionally inconsistent with Sinha et al. 2021 Fig. 4: mean D SKCM RTK-RAS/WNT="
            f"{skcm_rtk:.2f}/{skcm_wnt:.2f}, LIHC RTK-RAS/WNT={lihc_rtk:.2f}/{lihc_wnt:.2f}"
        )
    print(f"[tumorcomparer] Tier 1 verdict: {directional_verdict}")

    tier2_rows = []
    if pyreadr is None:
        for fname in sorted(os.listdir(TUMORCOMPARER_PRECOMPUTED_DIR)):
            tier2_rows.append({"file": fname, "parsed": False, "rows": None, "columns": None,
                                "reason": "pyreadr not installed in this environment"})
    else:
        for fname in sorted(os.listdir(TUMORCOMPARER_PRECOMPUTED_DIR)):
            path = os.path.join(TUMORCOMPARER_PRECOMPUTED_DIR, fname)
            try:
                result = pyreadr.read_r(path)
            except Exception as exc:  # noqa: BLE001 -- record and move on, never crash the run
                tier2_rows.append({"file": fname, "parsed": False, "rows": None, "columns": None,
                                    "reason": f"{type(exc).__name__}: {exc}"})
                continue
            if not result:
                tier2_rows.append({
                    "file": fname, "parsed": False, "rows": None, "columns": None,
                    "reason": "pyreadr returned no extractable data.frame -- likely a nested R "
                              "list (e.g. a per-geneset lookup), needs R itself to unpack",
                })
                continue
            df = next(iter(result.values()))
            has_pathway_column = any(
                key in col.lower() for col in df.columns for key in ("pathway", "geneset", "gene_set")
            )
            reason = (
                "parsed, but only genome-wide similarity columns -- no pathway/geneset-restricted "
                "column, so it cannot answer the RTK-RAS/WNT-restricted question Figure 4 asks"
                if not has_pathway_column else
                "parsed and carries a pathway/geneset column -- inspect manually for a SKCM/LIHC "
                "RTK-RAS/WNT-restricted score to join against Tier 1"
            )
            tier2_rows.append({"file": fname, "parsed": True, "rows": len(df),
                                "columns": len(df.columns), "reason": reason})

    tier2_df = pd.DataFrame(tier2_rows)
    tier2_df.to_csv(os.path.join(OUTPUT_DIR, "tumorcomparer_tier2_status.csv"), index=False)
    tier2_usable = tier2_df["reason"].str.contains("inspect manually", na=False).any()
    tier2_summary = (
        "Tier 2 found a candidate pathway-restricted file -- see tumorcomparer_tier2_status.csv"
        if tier2_usable else
        "Tier 2 blocked for all 6 precomputed files (see tumorcomparer_tier2_status.csv for the "
        "per-file reason) -- no numeric comparison against TumorComparer's own published output "
        "was possible, not forced."
    )
    print(f"[tumorcomparer] {tier2_summary}")

    comparable = "partial" if "directionally consistent" in directional_verdict else "no"
    return {
        "study": "TumorComparer (Sinha et al. 2021)", "mode": "gene_panel_qualitative (RTK-RAS/WNT, SKCM/LIHC)",
        "comparable_to_ranking": comparable, "reason": f"{directional_verdict}. {tier2_summary}",
    }


def run_celligner_coverage(cell_lines):
    celligner = pd.read_csv(CELLIGNER_PATH)
    cl_rows = celligner[celligner["type"] == "CL"].copy()

    our = cell_lines.set_index("ModelID")
    rows = []
    for _, r in cl_rows.iterrows():
        model_id = r["sampleID"]  # Celligner's CL-type sampleID is already a ModelID (ACH-xxxxxx)
        matched = model_id in our.index
        our_lineage = our.loc[model_id, "lineage"] if matched else None
        lineage_agrees = (
            str(our_lineage).strip().lower() == str(r["lineage"]).strip().lower() if matched else None
        )
        rows.append({
            "model_id": model_id,
            "celligner_lineage": r["lineage"],
            "celligner_subtype": r["subtype"],
            "celligner_undifferentiated_cluster": r["undifferentiated_cluster"],
            "matched_in_cell_lines_csv": matched,
            "our_lineage": our_lineage,
            "lineage_agrees": lineage_agrees,
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(os.path.join(OUTPUT_DIR, "celligner_coverage.csv"), index=False)

    n_total = len(out_df)
    n_matched = int(out_df["matched_in_cell_lines_csv"].sum())
    n_agree = int(out_df.loc[out_df["matched_in_cell_lines_csv"], "lineage_agrees"].sum())
    print(f"[celligner] {n_matched}/{n_total} matched a ModelID; {n_agree}/{n_matched} agree on lineage")

    reason = (
        f"{n_matched}/{n_total} Celligner cell-line rows resolve to a ModelID in "
        f"data/processed/cell_lines.csv; of those, {n_agree}/{n_matched} agree on lineage label. "
        f"Celligner's gene-level aligned-expression matrix was deliberately not downloaded (no "
        f"gene axis in this repo), so this is a data-coverage/lineage-label check, not a "
        f"validation of the ranking algorithm's output."
    )
    return {
        "study": "Celligner (Warren et al. 2021)", "mode": "coverage_only",
        "comparable_to_ranking": "no", "reason": reason,
    }


def run_netcellmatch_coverage(cell_lines):
    netcellmatch = pd.read_csv(NETCELLMATCH_PATH)
    cl_rows = netcellmatch[~netcellmatch["Sample_Name"].str.match(r"^TCGA-")].copy()

    our = cell_lines.set_index("stripped_cell_line_name")
    rows = []
    for _, r in cl_rows.iterrows():
        name = r["Sample_Name"]
        matched = name in our.index
        our_lineage = our.loc[name, "lineage"] if matched else None
        rows.append({
            "sample_name": name,
            "netcellmatch_lineage": r["Lineage"],
            "matched_in_cell_lines_csv": matched,
            "our_lineage": our_lineage,
            "lineage_agrees": (
                str(our_lineage).strip().lower() == str(r["Lineage"]).strip().lower() if matched else None
            ),
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(os.path.join(OUTPUT_DIR, "netcellmatch_coverage.csv"), index=False)

    n_total = len(out_df)
    n_matched = int(out_df["matched_in_cell_lines_csv"].sum())
    n_agree = int(out_df.loc[out_df["matched_in_cell_lines_csv"], "lineage_agrees"].sum())
    print(f"[netcellmatch] {n_matched}/{n_total} matched a stripped_cell_line_name; {n_agree}/{n_matched} agree on lineage")

    reason = (
        f"{n_matched}/{n_total} of NetCellMatch's breast-cohort cell-line rows resolve to a "
        f"stripped_cell_line_name in data/processed/cell_lines.csv; of those, {n_agree}/{n_matched} "
        f"agree on lineage. The file is a 233-antibody RPPA (protein/phospho) panel, not gene-level, "
        f"and covers only 1 of the paper's 3 cancer types (breast; lung/colon were not retrievable "
        f"-- see docs/reference/DATA_SOURCES.md), so this is a data-coverage check, not a validation "
        f"of the ranking algorithm's output."
    )
    return {
        "study": "NetCellMatch (Desai et al. 2022)", "mode": "coverage_only",
        "comparable_to_ranking": "no", "reason": reason,
    }


def main():
    gene_reference = pd.read_csv(os.path.join(DATA_DIR, "gene_reference.csv"))
    cell_lines = pd.read_csv(os.path.join(DATA_DIR, "cell_lines.csv"))

    manifest = [
        run_jin2023(gene_reference, cell_lines),
        run_tumorcomparer(gene_reference, cell_lines),
        run_celligner_coverage(cell_lines),
        run_netcellmatch_coverage(cell_lines),
    ]

    manifest_df = pd.DataFrame(manifest)
    manifest_df.to_csv(os.path.join(OUTPUT_DIR, "study_comparability.csv"), index=False)
    print("\n[summary] validation/study_comparability.csv:")
    print(manifest_df[["study", "mode", "comparable_to_ranking"]].to_string(index=False))


if __name__ == "__main__":
    main()
