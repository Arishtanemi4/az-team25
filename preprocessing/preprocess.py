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


if __name__ == "__main__":
    main()
