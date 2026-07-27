import pandas as pd


def build_profile_bridge(df8):
    return df8[["ProfileID", "ModelID", "Datatype"]].copy()


def build_hpa_bridge(df11, df7, df9):
    n_in = df11["Cell line"].nunique()

    step1 = df11.dropna(subset=["Cellosaurus ID"])[["Cell line", "Cellosaurus ID"]]
    step1 = step1.rename(columns={"Cellosaurus ID": "RRID"})

    known_accessions = set(df7["Accession (CVCL_xxxx)"].dropna())
    step2 = step1[step1["RRID"].isin(known_accessions)]

    df9_rrid = df9[["DepMap_ID", "RRID"]].dropna(subset=["RRID"])

    bridge = step2.merge(df9_rrid, on="RRID", how="inner", validate="one_to_many")
    bridge = bridge.rename(columns={"DepMap_ID": "ModelID"})
    bridge = bridge.drop_duplicates(subset="Cell line")[["Cell line", "ModelID"]]

    n_out = bridge["Cell line"].nunique()
    report = {
        "hpa_cell_lines_in": n_in,
        "hpa_cell_lines_resolved": n_out,
        "hpa_survival_pct": round(n_out / n_in * 100, 1),
    }
    return bridge, report
