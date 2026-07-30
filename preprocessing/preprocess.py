from pathlib import Path

import pandas as pd

import joins

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "data" / "raw"

PATHS = {
    "df7": RAW / "nomenclature" / "7_cellosaurus.csv",
}


def main():
    df7_raw = pd.read_csv(PATHS["df7"])
    df7 = joins.filter_human(df7_raw)


if __name__ == "__main__":
    main()
