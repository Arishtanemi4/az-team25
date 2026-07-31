from pathlib import Path

import pandas as pd

import joins

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "data" / "raw"

PATHS = {
    "df7": RAW / "nomenclature" / "7_cellosaurus.csv",
    "df8": RAW / "nomenclature" / "8_DepMap_OmicsProfiles.csv",
    "df9": RAW / "nomenclature" / "9_DepMap_sample_info.csv",
    "df11": RAW / "nomenclature" / "11_hpa_rna_celline_description.tsv",
    "df3": RAW / "gene_expression" / "3_GEOexpression.txt",
    "df10": RAW / "nomenclature" / "10_GEOInfo.txt",
    "df1": RAW / "gene_expression" / "1_4_hpa_rna_celline.tsv",
    "df2": RAW / "gene_expression" / "2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.csv",
    "df4": RAW / "gene_expression" / "4_Harmonized_MS_CCLE_Gygi_subsetted.csv",
}


def main():
    df7_raw = pd.read_csv(PATHS["df7"])
    df7 = joins.filter_human(df7_raw)

    df9 = pd.read_csv(PATHS["df9"])
    df9, n_df9_dupes = joins.resolve_duplicate_keys(df9, ["DepMap_ID"], "RRID")

    df9_age_avg, n_age_dupes = joins.average_duplicate_rna_profiles(df9, ["patient_id"], "age")

    df8 = pd.read_csv(PATHS["df8"])
    profile_bridge = joins.build_profile_bridge(df8)

    ccle_bridge = joins.build_ccle_bridge(df9)

    df11 = pd.read_csv(PATHS["df11"], sep="\t")
    hpa_bridge, hpa_report = joins.build_hpa_bridge(df11, df7, df9)

    df10 = pd.read_csv(PATHS["df10"], sep="\t")
    geo_gsm_columns = pd.read_csv(PATHS["df3"], sep="\t", nrows=0).columns[1:]
    geo_bridge, geo_report = joins.build_geo_bridge(geo_gsm_columns, df10, df7, df9)

    gene_reference, gene_ref_report = joins.build_gene_reference(PATHS["df1"], PATHS["df2"])

    df4 = pd.read_csv(PATHS["df4"])
    protein_cols = [c for c in df4.columns if c != df4.columns[0]]
    protein_bridge, protein_bridge_report = joins.build_symbol_bridge(
        protein_cols, gene_reference, symbol_position="inside"
    )


if __name__ == "__main__":
    main()
