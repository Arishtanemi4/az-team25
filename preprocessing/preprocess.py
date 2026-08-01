import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import clean
import impute
import joins
import normalize
import scale

GENE_HEADER_PATTERN = re.compile(r"^(.*) \((ENSG\d+)\)$")

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "data" / "raw"
AUGMENTED = REPO_ROOT / "data" / "augmented"
PROCESSED = REPO_ROOT / "data" / "processed"
RESOURCES = REPO_ROOT / "preprocessing" / "resources"


PATHS = {
    "df1": RAW / "gene_expression" / "1_4_hpa_rna_celline.tsv",
    "df2": RAW / "gene_expression" / "2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.csv",
    "df3": RAW / "gene_expression" / "3_GEOexpression.txt",
    "df4": RAW / "gene_expression" / "4_Harmonized_MS_CCLE_Gygi_subsetted.csv",
    "df5": RAW / "gene_properties" / "5_OmicsFusionFilteredSupplementary.csv",
    "df6": RAW / "gene_properties" / "6_OmicsSomaticMutationsProfile.csv",
    "df7": RAW / "nomenclature" / "7_cellosaurus.csv",
    "df8": RAW / "nomenclature" / "8_DepMap_OmicsProfiles.csv",
    "df9": RAW / "nomenclature" / "9_DepMap_sample_info.csv",
    "df10": RAW / "nomenclature" / "10_GEOInfo.txt",
    "df11": RAW / "nomenclature" / "11_hpa_rna_celline_description.tsv",
    "df12": RAW / "non_gene_expression" / "12_CCLE_metabolomics_20190502.csv",
    "df13": RAW / "non_gene_expression" / "13_CCLE_miRNA_20181103.gct",
    "df14": RAW / "non_gene_expression" / "14_OmicsGlobalSignatures.csv",
    "df15": AUGMENTED / "gene_properties" / "15_CRISPRGeneEffect.csv",
    "df16": AUGMENTED / "gene_properties" / "16_OmicsCNGeneWGS.csv",
    "df17": AUGMENTED / "nomenclature" / "17_Model.csv",
}
GENE_BIOTYPE_REFERENCE = RESOURCES / "gene_biotype_reference.csv"


def resolve_profile_ids(long_df, profile_bridge, datatypes):
    bridge = profile_bridge[profile_bridge["Datatype"].isin(datatypes)][["ProfileID", "ModelID"]]

    return long_df.merge(bridge, on="ProfileID", how="inner", validate="many_to_one")


def scan_model_ids(path, **read_kwargs):
    return set(pd.read_csv(path, **read_kwargs).iloc[:, 0].astype(str))


def two_state_coverage(model_ids, assayed_ids):
    return model_ids.isin(assayed_ids).map({True: "measured", False: "not_assayed"})


def three_state_coverage(model_ids, assayed_ids, has_event_ids):
    state = pd.Series("not_assayed", index=model_ids.index)
    assayed_mask = model_ids.isin(assayed_ids)
    state[assayed_mask] = "measured_absent"
    state[assayed_mask & model_ids.isin(has_event_ids)] = "measured"
    return state


def write_table(df, filename):
    df.to_csv(PROCESSED / filename, index=False)
    df.to_parquet(PROCESSED / filename.replace(".csv", ".parquet"), index=False)


def main():
    files_read = set()
    report = {}
    PROCESSED.mkdir(parents=True, exist_ok=True)

    print("Loading nomenclature files...")
    df7_raw = pd.read_csv(PATHS["df7"])
    df7 = joins.filter_human(df7_raw)
    df8 = pd.read_csv(PATHS["df8"])
    df9 = pd.read_csv(PATHS["df9"])
    df10 = pd.read_csv(PATHS["df10"], sep="\t")
    df11 = pd.read_csv(PATHS["df11"], sep="\t")
    df17 = pd.read_csv(PATHS["df17"])
    files_read |= {"df7", "df8", "df9", "df10", "df11", "df17"}
    report["species_gate"] = {"cellosaurus_rows_in": len(df7_raw), "human_rows": len(df7)}

    profile_bridge = joins.build_profile_bridge(df8)
    ccle_bridge = joins.build_ccle_bridge(df9)
    hpa_bridge, hpa_report = joins.build_hpa_bridge(df11, df7, df9)
    report["hpa_bridge"] = hpa_report

    geo_gsm_columns = pd.read_csv(PATHS["df3"], sep="\t", nrows=0).columns[1:]
    geo_bridge, geo_report = joins.build_geo_bridge(geo_gsm_columns, df10, df7, df9)
    report["geo_bridge"] = geo_report

    print("Building gene_reference (table 2)...")
    gene_reference, gene_ref_report = joins.build_gene_reference(PATHS["df1"], PATHS["df2"])
    files_read |= {"df1", "df2"}
    biotype_reference = pd.read_csv(GENE_BIOTYPE_REFERENCE)

    gene_reference = gene_reference.merge(biotype_reference, on="ensembl_id", how="left", validate="one_to_one")
    gene_ref_report["n_biotype_unresolved"] = int(gene_reference["gene_biotype"].isna().sum())
    gene_reference["gene_biotype"] = gene_reference["gene_biotype"].fillna("Unknown")
    report["gene_reference"] = gene_ref_report

    print("Scanning measurement files for the identity spine, before any full clean/melt...")
    rna_assayed_ids = set(df8.loc[df8["Datatype"] == "rna", "ModelID"])
    dna_assayed_ids = set(df8.loc[df8["Datatype"].isin(["wes", "wgs"]), "ModelID"])
    # DepMap fusion calls come from RNA-seq, not WES/WGS (data/raw/gene_properties/_.md item 5,
    # "SequencingID = PR-... (the RNA profile)") -- the RNA-assayed pool is the right gate.
    fusion_assayed_ids = rna_assayed_ids
    protein_assayed_ids = scan_model_ids(PATHS["df4"], usecols=[0])
    dependency_assayed_ids = scan_model_ids(PATHS["df15"], usecols=[0])
    copy_number_assayed_ids = scan_model_ids(PATHS["df16"], usecols=["ModelID"])
    signatures_assayed_ids = scan_model_ids(PATHS["df14"], usecols=["ModelID"])
    files_read |= {"df4", "df15", "df16", "df14"}

    omics_model_ids = (
        rna_assayed_ids
        | dna_assayed_ids
        | protein_assayed_ids
        | fusion_assayed_ids
        | dependency_assayed_ids
        | copy_number_assayed_ids
        | signatures_assayed_ids
    )

    print("Building cell_lines (table 1) -- the spine, before any measurement table is cleaned...")
    cell_lines, cell_lines_report = joins.build_cell_lines(df9, df7, df11, df17, omics_model_ids)
    cell_lines = clean.clean_cell_lines_meta(cell_lines)
    report["cell_lines"] = cell_lines_report
    valid_model_ids = set(cell_lines["ModelID"])

    print("Cleaning expression_rna (table 3) -- this is the big one, ~80M rows...")
    df2 = pd.read_csv(PATHS["df2"])
    expr_rna = clean.clean_expression_rna(df2)
    del df2
    ensembl_from_header = {}
    for col in expr_rna["gene_col"].unique():
        match = GENE_HEADER_PATTERN.match(col)
        if match:
            ensembl_from_header[col] = match.group(2)
    expr_rna["ensembl_id"] = expr_rna["gene_col"].map(ensembl_from_header)
    expr_rna = expr_rna.dropna(subset=["ensembl_id"]).drop(columns=["gene_col"])
    expr_rna = resolve_profile_ids(expr_rna, profile_bridge, ["rna"])

    expr_rna, n_expr_rna_dupes = joins.average_duplicate_rna_profiles(
        expr_rna, ["ModelID", "ensembl_id"], "log2tpm1"
    )
    report["expression_rna_duplicate_keys_resolved"] = n_expr_rna_dupes
    expr_rna = scale.scale_expression_rna(expr_rna)
    expr_rna = expr_rna[expr_rna["ModelID"].isin(valid_model_ids)]
    expr_rna = expr_rna[["ModelID", "ensembl_id", "log2tpm1"]]

    print("Cleaning expression_rna_hpa (table 4)...")
    df1 = pd.read_csv(PATHS["df1"], sep="\t")
    expr_hpa = clean.clean_expression_hpa(df1)
    del df1

    expr_hpa = expr_hpa.merge(hpa_bridge, on="Cell line", how="inner", validate="many_to_one").drop(columns=["Cell line"])
    expr_hpa = normalize.normalize_hpa(expr_hpa)
    expr_hpa = expr_hpa[expr_hpa["ModelID"].isin(valid_model_ids)]
    expr_hpa = expr_hpa[["ModelID", "ensembl_id", "TPM", "pTPM", "nTPM", "log2ntpm1"]].rename(
        columns={"TPM": "tpm", "pTPM": "ptpm", "nTPM": "ntpm"}
    )

    print("Cleaning expression_rna_geo (new table, corroboration/breadth only)...")
    df3 = pd.read_csv(PATHS["df3"], sep="\t")
    files_read.add("df3")
    expr_geo = clean.clean_expression_geo(df3)
    del df3

    expr_geo = expr_geo.merge(geo_bridge, on="Geo_accession", how="inner", validate="many_to_one").drop(columns=["Geo_accession"])
    expr_geo = scale.scale_expression_geo(expr_geo)
    expr_geo = expr_geo[expr_geo["ModelID"].isin(valid_model_ids)]
    expr_geo = expr_geo[["ModelID", "ensembl_id", "log1p_expr"]]

    print("Cleaning protein (table 5)...")
    df4 = pd.read_csv(PATHS["df4"])
    protein_cols = [c for c in df4.columns if c != df4.columns[0]]
    protein_bridge, protein_bridge_report = joins.build_symbol_bridge(
        protein_cols, gene_reference, symbol_position="inside"
    )
    report["protein_symbol_bridge"] = protein_bridge_report
    protein = clean.clean_protein(df4, protein_bridge)
    del df4
    protein, n_protein_dupes = joins.resolve_duplicate_keys(
        protein, ["ModelID", "ensembl_id"], "protein_col"
    )
    report["protein_duplicate_keys_resolved"] = n_protein_dupes
    protein, impute_report = impute.impute_protein_detection_floor(protein)
    report["protein_imputation"] = impute_report
    protein = scale.scale_protein(protein)
    protein = protein[protein["ModelID"].isin(valid_model_ids)]
    protein = protein[["ModelID", "ensembl_id", "zscore", "detected"]]

    print("Cleaning mutations (table 6)...")
    df6 = pd.read_csv(PATHS["df6"], low_memory=False)
    files_read.add("df6")
    mutations = clean.clean_mutations(df6)
    del df6
    mutations = resolve_profile_ids(mutations, profile_bridge, ["wes", "wgs"])
    valid_genes = set(gene_reference["ensembl_id"])
    n_before = len(mutations)
    mutations = mutations[mutations["ensembl_id"].isin(valid_genes)]
    report["mutations_gene_filter"] = {"dropped_unresolved_ensembl_id": n_before - len(mutations)}
    mutations = mutations[mutations["ModelID"].isin(valid_model_ids)]

    print("Cleaning fusions (table 7)...")
    df5 = pd.read_csv(PATHS["df5"])
    files_read.add("df5")
    fusions = clean.clean_fusions(df5)
    del df5
    n_before = len(fusions)
    fusions = fusions[fusions["ensembl_id"].isin(valid_genes)]
    fusions["partner_ensembl_id"] = fusions["partner_ensembl_id"].where(
        fusions["partner_ensembl_id"].isin(valid_genes)
    )
    report["fusions_gene_filter"] = {"dropped_unresolved_ensembl_id": n_before - len(fusions)}
    fusions = fusions[fusions["ModelID"].isin(valid_model_ids)]

    print("Cleaning dependency (table 8)...")
    df15_header = pd.read_csv(PATHS["df15"], nrows=0)
    dep_cols = [c for c in df15_header.columns if c != df15_header.columns[0]]
    dep_bridge, dep_bridge_report = joins.build_symbol_bridge(
        dep_cols, gene_reference, symbol_position="before"
    )
    report["dependency_symbol_bridge"] = dep_bridge_report
    df15 = pd.read_csv(PATHS["df15"])
    dependency = clean.clean_dependency(df15, dep_bridge)
    del df15
    dependency = dependency[dependency["ModelID"].isin(valid_model_ids)]

    print("Cleaning copy_number (table 9)...")
    df16_header = pd.read_csv(PATHS["df16"], nrows=0)
    cn_cols = [c for c in df16_header.columns if c not in
               ("SequencingID", "ModelConditionID", "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel")]
    cn_bridge, cn_bridge_report = joins.build_symbol_bridge(
        cn_cols, gene_reference, symbol_position="before"
    )
    report["copy_number_symbol_bridge"] = cn_bridge_report
    df16 = pd.read_csv(PATHS["df16"])
    copy_number = clean.clean_copy_number(df16, cn_bridge)
    del df16
    copy_number = copy_number[copy_number["ModelID"].isin(valid_model_ids)]

    print("Cleaning genome signatures (context, feeds coverage.csv)...")
    df14 = pd.read_csv(PATHS["df14"])
    signatures = clean.clean_genome_signatures(df14)
    del df14
    signatures = signatures[signatures["ModelID"].isin(valid_model_ids)]

    print("Cleaning metabolomics (new full-value table, context only -- fails A1)...")
    df12 = pd.read_csv(PATHS["df12"])
    metabolomics = clean.clean_metabolomics(df12, ccle_bridge)
    metabolomics = scale.scale_metabolomics(metabolomics)
    metabolomics = metabolomics[metabolomics["ModelID"].isin(valid_model_ids)]
    del df12
    files_read.add("df12")

    print("Cleaning miRNA (new table, context only -- fails A1)...")
    df13 = pd.read_csv(PATHS["df13"], sep="\t", skiprows=2)
    files_read.add("df13")
    mirna = clean.clean_mirna(df13, ccle_bridge)
    del df13
    mirna = scale.scale_mirna(mirna)
    mirna = mirna[mirna["ModelID"].isin(valid_model_ids)]
    mirna = mirna[["ModelID", "mirna_name", "raw_value", "log1p_value"]]

    print("Building coverage (three-state) ...")
    coverage = cell_lines[["ModelID"]].copy()
    coverage["expression_rna_state"] = two_state_coverage(coverage["ModelID"], rna_assayed_ids)
    coverage["expression_rna_hpa_state"] = two_state_coverage(coverage["ModelID"], set(hpa_bridge["ModelID"]))
    coverage["expression_rna_geo_state"] = two_state_coverage(coverage["ModelID"], set(geo_bridge["ModelID"]))
    coverage["protein_state"] = two_state_coverage(coverage["ModelID"], protein_assayed_ids)
    coverage["dependency_state"] = two_state_coverage(coverage["ModelID"], dependency_assayed_ids)
    coverage["copy_number_state"] = two_state_coverage(coverage["ModelID"], copy_number_assayed_ids)
    coverage["genome_signatures_state"] = two_state_coverage(coverage["ModelID"], signatures_assayed_ids)
    coverage["mutations_state"] = three_state_coverage(
        coverage["ModelID"], dna_assayed_ids, set(mutations["ModelID"])
    )
    coverage["fusions_state"] = three_state_coverage(
        coverage["ModelID"], fusion_assayed_ids, set(fusions["ModelID"])
    )
    coverage["metabolomics_available"] = coverage["ModelID"].isin(metabolomics["ModelID"])
    coverage["mirna_available"] = coverage["ModelID"].isin(mirna["ModelID"])

    coverage = coverage.merge(
        signatures[["ModelID", "msi_score", "msi_high", "loh_fraction", "wgd", "cin", "ploidy", "wgs_available"]],
        on="ModelID",
        how="left",
        validate="one_to_one",
    )

    print("Writing processed tables + parquet mirrors to data/processed/...")
    tables = {
        "cell_lines.csv": cell_lines,
        "gene_reference.csv": gene_reference,
        "expression_rna.csv": expr_rna,
        "expression_rna_hpa.csv": expr_hpa,
        "expression_rna_geo.csv": expr_geo,
        "protein.csv": protein,
        "mutations.csv": mutations,
        "fusions.csv": fusions,
        "dependency.csv": dependency,
        "copy_number.csv": copy_number,
        "genome_signatures.csv": signatures,
        "mirna.csv": mirna,
        "metabolomics.csv": metabolomics,
        "coverage.csv": coverage,
    }
    for filename, table in tables.items():
        write_table(table, filename)

    print("Checking every primary/augmented file was read...")
    unread = set(PATHS) - files_read
    assert not unread, f"preprocess.py never read: {sorted(unread)}"

    dupe_check = {
        "expression_rna": int(expr_rna.duplicated(subset=["ModelID", "ensembl_id"]).sum()),
        "protein": int(protein.duplicated(subset=["ModelID", "ensembl_id"]).sum()),
    }
    assert all(v == 0 for v in dupe_check.values()), f"(ModelID, ensembl_id) still duplicated: {dupe_check}"

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "primary_files_read": sorted(files_read),
        "primary_files_expected": sorted(PATHS),
        "counts": {name: len(table) for name, table in tables.items()},
        "coverage_state_counts": {
            col: coverage[col].value_counts().to_dict()
            for col in coverage.columns
            if col.endswith("_state") or col.endswith("_available")
        },
        "join_bridge_report": report,
        "duplicate_key_check": dupe_check,
    }
    (PROCESSED / "build_manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")

    print("\n=== Consolidated report ===")
    for name, table in tables.items():
        print(f"{name}: {len(table):,} rows")
    print("\nJoin/bridge reports:")
    for key, value in report.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
