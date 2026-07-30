from pathlib import Path

import pandas as pd

import joins

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "data" / "raw"

PATHS = {
    "df7": RAW / "nomenclature" / "7_cellosaurus.csv",
    "df8": RAW / "nomenclature" / "8_DepMap_OmicsProfiles.csv",
    "df9": RAW / "nomenclature" / "9_DepMap_sample_info.csv",
}


def main():
    df7_raw = pd.read_csv(PATHS["df7"])
    df7 = joins.filter_human(df7_raw)

    df9 = pd.read_csv(PATHS["df9"])
    df9, n_df9_dupes = joins.resolve_duplicate_keys(df9, ["DepMap_ID"], "RRID")

    df8 = pd.read_csv(PATHS["df8"])
    profile_bridge = joins.build_profile_bridge(df8)

    ccle_bridge = joins.build_ccle_bridge(df9)


if __name__ == "__main__":
    main()
