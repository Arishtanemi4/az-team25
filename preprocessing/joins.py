import re

import pandas as pd

PROBLEM_PATTERN = re.compile(r"problematic|contaminated|misidentified", re.IGNORECASE)


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


def build_geo_bridge(geo_gsm_columns, df10, df7, df9):
    n_in = len(geo_gsm_columns)

    step1 = df10.dropna(subset=["Cellosaurus_ID"])[["Geo_accession", "Cellosaurus_ID", "Matching_Type"]]
    step1 = step1.rename(columns={"Cellosaurus_ID": "RRID"})
    step1 = step1[step1["Geo_accession"].isin(geo_gsm_columns)]

    known_accessions = set(df7["Accession (CVCL_xxxx)"].dropna())
    step2 = step1[step1["RRID"].isin(known_accessions)]

    df9_rrid = df9[["DepMap_ID", "RRID"]].dropna(subset=["RRID"])

    max_right_multiplicity = int(df9_rrid["RRID"].value_counts().max())
    bridge = step2.merge(df9_rrid, on="RRID", how="inner", validate="many_to_many")
    assert len(bridge) <= len(step2) * max_right_multiplicity, (
        f"GEO/RRID merge fanned out beyond the per-key multiplicity bound: "
        f"{len(bridge)} > {len(step2)} x {max_right_multiplicity}"
    )
    bridge = bridge.rename(columns={"DepMap_ID": "ModelID"})
    bridge = bridge.drop_duplicates(subset="Geo_accession")

    n_out = bridge["Geo_accession"].nunique()
    report = {
        "geo_gsm_in": n_in,
        "geo_gsm_resolved": n_out,
        "geo_survival_pct": round(n_out / n_in * 100, 1),
        "geo_matching_type_counts_resolved": bridge["Matching_Type"].value_counts().to_dict(),
    }
    return bridge[["Geo_accession", "ModelID"]], report


def filter_human(df7):
    return df7[df7["Species of origin"].astype(str).str.contains("Homo sapiens", na=False)].copy()


def build_ccle_bridge(df9):
    return (
        df9[["DepMap_ID", "CCLE_Name"]]
        .dropna(subset=["CCLE_Name"])
        .rename(columns={"DepMap_ID": "ModelID", "CCLE_Name": "CCLE_ID"})
        .drop_duplicates(subset="CCLE_ID")
    )


def resolve_duplicate_keys(df, key_cols, tie_break_col):
    n_before = len(df)
    deduped = df.sort_values(tie_break_col).drop_duplicates(subset=key_cols, keep="first")
    deduped = deduped.reset_index(drop=True)
    return deduped, n_before - len(deduped)


def average_duplicate_rna_profiles(df, key_cols, value_col):
    n_before = len(df)
    averaged = df.groupby(key_cols, as_index=False)[value_col].mean()
    return averaged, n_before - len(averaged)


def build_gene_reference(df1_path, df2_path):
    df1_genes = pd.read_csv(df1_path, sep="\t", usecols=["Gene", "Gene name"])
    df1_genes = df1_genes.drop_duplicates().rename(
        columns={"Gene": "ensembl_id", "Gene name": "symbol"}
    )

    header = pd.read_csv(df2_path, nrows=0)
    pairs = []
    for col in header.columns[1:]:
        match = re.match(r"^(.*) \((ENSG\d+)\)$", col)
        if match:
            pairs.append((match.group(2), match.group(1)))
    df2_genes = pd.DataFrame(pairs, columns=["ensembl_id", "symbol"]).drop_duplicates()

    all_genes = pd.concat([df1_genes, df2_genes], ignore_index=True).drop_duplicates()

    symbol_counts = all_genes.groupby("symbol")["ensembl_id"].nunique()
    ambiguous_symbols = symbol_counts[symbol_counts > 1].index.tolist()

    gene_reference = all_genes.drop_duplicates(subset="ensembl_id", keep="first")
    gene_reference = gene_reference.reset_index(drop=True)

    report = {
        "n_genes": len(gene_reference),
        "n_ambiguous_symbols": len(ambiguous_symbols),
        "ambiguous_symbols_sample": ambiguous_symbols[:15],
    }
    return gene_reference, report


def build_symbol_bridge(header_cols, gene_reference, symbol_position="before"):
    symbol_to_ensembl = gene_reference.groupby("symbol")["ensembl_id"].apply(list).to_dict()

    col_to_ensembl = {}
    no_match = 0
    collisions = 0
    collision_candidates = {}
    for col in header_cols:
        if symbol_position == "before":
            symbol = col.split(" (")[0]
        else:
            match = re.search(r"\(([^()]+)\)\s*$", col)
            symbol = match.group(1) if match else None

        candidates = symbol_to_ensembl.get(symbol) if symbol else None
        if not candidates:
            no_match += 1
        elif len(candidates) > 1:
            collisions += 1
            collision_candidates[symbol] = candidates
        else:
            col_to_ensembl[col] = candidates[0]

    report = {
        "n_columns": len(header_cols),
        "n_resolved": len(col_to_ensembl),
        "n_no_match": no_match,
        "n_collisions": collisions,
        "collision_candidates": collision_candidates,
    }
    return col_to_ensembl, report


def build_cell_lines(df9, df7, df11, df17, omics_model_ids):
    df9 = df9.rename(columns={"DepMap_ID": "ModelID"}).copy()

    df7_clean = df7.dropna(subset=["Identifier (cell line name)"]).copy()
    df7_clean["is_problematic"] = (
        df7_clean["Comments"].fillna("").str.contains(PROBLEM_PATTERN)
    )
    problem_accessions = set(
        df7_clean.loc[df7_clean["is_problematic"], "Accession (CVCL_xxxx)"]
    )
    df9["is_problematic"] = df9["RRID"].isin(problem_accessions)

    df9_ids = set(df9["ModelID"])
    orphans = set(omics_model_ids) - df9_ids
    df17_ids = set(df17["ModelID"])
    resolved_orphans = orphans & df17_ids
    still_excluded = orphans - df17_ids

    rename_map = {
        "OncotreeLineage": "lineage",
        "CellLineName": "cell_line_name",
        "StrippedCellLineName": "stripped_cell_line_name",
        "OncotreePrimaryDisease": "primary_disease",
        "PatientID": "patient_id",
        "Sex": "sex",
        "PrimaryOrMetastasis": "primary_or_metastasis",
    }
    backfill = df17[df17["ModelID"].isin(resolved_orphans)].copy()
    backfill = backfill.rename(columns=rename_map)
    if "sex" in backfill.columns:
        backfill["sex"] = backfill["sex"].str.lower()
    keep_cols = ["ModelID", "RRID"] + [c for c in rename_map.values() if c in backfill.columns]
    backfill = backfill[keep_cols]
    backfill["is_problematic"] = backfill["RRID"].isin(problem_accessions)

    cell_lines = pd.concat([df9, backfill], ignore_index=True, sort=False)

    report = {
        "n_orphans_total": len(orphans),
        "n_orphans_resolved_by_df17": len(resolved_orphans),
        "n_orphans_still_excluded": len(still_excluded),
    }
    return cell_lines, report
