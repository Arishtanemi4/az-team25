import pandas as pd

UNKNOWN_FILL_COLUMNS = ["primary_or_metastasis", "lineage_sub_subtype", "Subtype"]


def clean_cell_lines_meta(cell_lines):
    df = cell_lines.copy()
    for col in UNKNOWN_FILL_COLUMNS:
        if col in df.columns:
            df[col] = df[col].fillna("Unknown")
    return df


def clean_expression_rna(df2):
    df = df2.rename(columns={df2.columns[0]: "ProfileID"})
    long = df.melt(id_vars="ProfileID", var_name="gene_col", value_name="log2tpm1")
    return long.dropna(subset=["log2tpm1"])


def clean_expression_hpa(df1):
    return df1[["Gene", "Cell line", "TPM", "pTPM", "nTPM"]].rename(columns={"Gene": "ensembl_id"})


def clean_protein(df4, col_to_ensembl):
    df = df4.rename(columns={df4.columns[0]: "ModelID"})
    protein_cols = [c for c in df.columns if c != "ModelID"]
    keep_cols = [c for c in protein_cols if c in col_to_ensembl]

    long = df[["ModelID"] + keep_cols].melt(
        id_vars="ModelID", var_name="protein_col", value_name="zscore"
    )
    long["ensembl_id"] = long["protein_col"].map(col_to_ensembl)
    long["detected"] = long["zscore"].notna()
    return long[["ModelID", "ensembl_id", "protein_col", "zscore", "detected"]]


def clean_mutations(df6):
    return df6.rename(
        columns={
            "EnsemblGeneID": "ensembl_id",
            "ProteinChange": "protein_change",
            "VepImpact": "vep_impact",
            "HessDriver": "is_driver",
            "Hotspot": "is_hotspot",
        }
    )


def clean_fusions(df5):
    df = df5[df5["IsDefaultEntryForModel"] == "Yes"].copy()
    df = df.drop(columns=["Unnamed: 0"], errors="ignore")
    df["gene1_ensembl"] = df["gene1(ENS ID)"].str.extract(r"\((ENSG\d+)")
    df["gene2_ensembl"] = df["gene2(ENS ID)"].str.extract(r"\((ENSG\d+)")
    df["confidence_high"] = df["confidence"] == "high"
    df["in_frame"] = df["reading_frame"] == "in-frame"

    shared_cols = [c for c in df.columns if c not in ("gene1_ensembl", "gene2_ensembl")]
    forward = df.rename(
        columns={"gene1_ensembl": "ensembl_id", "gene2_ensembl": "partner_ensembl_id"}
    )
    reverse = df.rename(
        columns={"gene2_ensembl": "ensembl_id", "gene1_ensembl": "partner_ensembl_id"}
    )
    cols = shared_cols + ["ensembl_id", "partner_ensembl_id"]
    long = pd.concat([forward[cols], reverse[cols]], ignore_index=True)
    return long.dropna(subset=["ensembl_id"])


def clean_dependency(df15, col_to_ensembl):
    df = df15.rename(columns={df15.columns[0]: "ModelID"})
    gene_cols = [c for c in col_to_ensembl if c in df.columns]
    long = df[["ModelID"] + gene_cols].melt(
        id_vars="ModelID", var_name="gene_col", value_name="dependency_score"
    )
    long["ensembl_id"] = long["gene_col"].map(col_to_ensembl)
    long = long.dropna(subset=["dependency_score"])
    return long[["ModelID", "ensembl_id", "dependency_score"]]


def clean_copy_number(df16, col_to_ensembl):
    df = df16[df16["IsDefaultEntryForModel"] == "Yes"].copy()
    id_vars = ["ModelID", "SequencingID", "ModelConditionID", "IsDefaultEntryForMC"]
    gene_cols = [c for c in col_to_ensembl if c in df.columns]
    long = df[id_vars + gene_cols].melt(
        id_vars=id_vars, var_name="gene_col", value_name="copy_number"
    )
    long["ensembl_id"] = long["gene_col"].map(col_to_ensembl)
    long = long.dropna(subset=["copy_number"])
    long["amp"] = long["copy_number"] > 1.5
    long["del"] = long["copy_number"] < 0.5
    return long[id_vars + ["ensembl_id", "copy_number", "amp", "del"]]


def clean_genome_signatures(df14):
    df = df14[df14["IsDefaultEntryForModel"] == "Yes"].copy()
    df["msi_high"] = df["MSIScore"] >= 20
    df["wgs_available"] = df[["LoHFraction", "WGD", "CIN", "Ploidy"]].notna().all(axis=1)
    keep = [
        "ModelID", "SequencingID", "ModelConditionID",
        "MSIScore", "msi_high", "LoHFraction", "WGD", "CIN", "Ploidy", "Aneuploidy",
        "wgs_available",
    ]
    return df[keep].rename(
        columns={
            "MSIScore": "msi_score",
            "LoHFraction": "loh_fraction",
            "WGD": "wgd",
            "CIN": "cin",
            "Ploidy": "ploidy",
            "Aneuploidy": "aneuploidy",
        }
    )


def clean_metabolomics(df12, ccle_bridge):
    ccle_to_model = ccle_bridge.set_index("CCLE_ID")["ModelID"]
    df = df12.copy()
    df["ModelID"] = df["DepMap_ID"].fillna(df["CCLE_ID"].map(ccle_to_model))
    metabolite_cols = [c for c in df.columns if c not in ("DepMap_ID", "CCLE_ID", "ModelID")]
    long = df[["ModelID"] + metabolite_cols].melt(
        id_vars="ModelID", var_name="metabolite_name", value_name="value"
    )
    return long.dropna(subset=["ModelID", "value"])


def clean_expression_geo(df3):
    df = df3.rename(columns={"Gene": "ensembl_id"})
    long = df.melt(id_vars="ensembl_id", var_name="Geo_accession", value_name="raw_expr")
    return long.dropna(subset=["raw_expr"])


def clean_mirna(df13, ccle_bridge):
    df = df13.drop(columns=["Name"]).set_index("Description")
    long = df.T.reset_index().rename(columns={"index": "CCLE_ID"})
    long = long.melt(id_vars="CCLE_ID", var_name="mirna_name", value_name="raw_value")
    long = long.dropna(subset=["raw_value"])
    ccle_to_model = ccle_bridge.set_index("CCLE_ID")["ModelID"]
    long["ModelID"] = long["CCLE_ID"].map(ccle_to_model)
    return long.dropna(subset=["ModelID"]).drop(columns=["CCLE_ID"])
