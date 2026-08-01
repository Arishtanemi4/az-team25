from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
VOGELSTEIN_PATH = REPO_ROOT / "preprocessing" / "resources" / "vogelstein_2013_table_s2a.csv"
BAILEY_PATH = REPO_ROOT / "preprocessing" / "resources" / "bailey_2018_table_s1.csv"
GENE_REFERENCE_PATH = REPO_ROOT / "data" / "processed" / "gene_reference.csv"
OUT_PATH = REPO_ROOT / "preprocessing" / "resources" / "gene_role_reference.csv"


BAILEY_FUNCTION_TO_CLASS = {
    "oncogene": "oncogene",
    "possible oncogene": "oncogene",
    "tsg": "tumour_suppressor",
    "possible tsg": "tumour_suppressor",
}

VOGELSTEIN_CLASS_TO_CLASS = {"Oncogene": "oncogene", "TSG": "tumour_suppressor"}


def resolve_bailey(bailey_df):
    resolved, unresolved_split = {}, []
    for symbol, rows in bailey_df.groupby("Symbol"):
        official = rows[rows["Decision"] == "official"]
        pool = official if len(official) > 0 else rows  # fall back to rescued only if no official row exists

        pancan = pool[pool["Cancer"] == "PANCAN"]
        pancan_classes = {BAILEY_FUNCTION_TO_CLASS[f] for f in pancan["Function"] if f in BAILEY_FUNCTION_TO_CLASS}
        if len(pancan_classes) == 1:
            resolved[symbol] = pancan_classes.pop()
            continue

        tissue_classes = [BAILEY_FUNCTION_TO_CLASS[f] for f in pool["Function"] if f in BAILEY_FUNCTION_TO_CLASS]
        distinct = set(tissue_classes)
        if len(distinct) == 1:
            resolved[symbol] = distinct.pop()
        elif len(distinct) > 1:
            unresolved_split.append(symbol)
        # len(distinct) == 0: no usable Function anywhere for this gene -- simply not resolved,
        # not logged as a conflict (there is no disagreement, just no signal).
    return resolved, unresolved_split


def resolve_symbols(symbols, gene_reference_df):
    counts = gene_reference_df.groupby("symbol")["ensembl_id"].nunique()
    unique_symbols = set(counts[counts == 1].index)
    symbol_to_id = gene_reference_df[gene_reference_df["symbol"].isin(unique_symbols)].set_index("symbol")["ensembl_id"].to_dict()

    resolved, unresolved, ambiguous = {}, [], []
    for symbol in symbols:
        if symbol in symbol_to_id:
            resolved[symbol] = symbol_to_id[symbol]
        elif symbol in counts.index:  # present, but with >1 ensembl_id
            ambiguous.append(symbol)
        else:
            unresolved.append(symbol)
    return resolved, unresolved, ambiguous


def build_gene_role_reference(vogelstein_df, bailey_df, gene_reference_df):
    vogelstein_by_symbol = dict(zip(vogelstein_df["symbol"], vogelstein_df["classification"].map(VOGELSTEIN_CLASS_TO_CLASS)))
    bailey_by_symbol, bailey_unresolved_split = resolve_bailey(bailey_df)

    disagreements = [
        (symbol, vogelstein_by_symbol[symbol], bailey_by_symbol[symbol])
        for symbol in set(vogelstein_by_symbol) & set(bailey_by_symbol)
        if vogelstein_by_symbol[symbol] != bailey_by_symbol[symbol]
    ]

    gene_class_by_symbol = dict(vogelstein_by_symbol)
    gene_class_by_symbol.update(bailey_by_symbol)

    symbol_to_id, unresolved, ambiguous = resolve_symbols(gene_class_by_symbol.keys(), gene_reference_df)
    rows = [
        {"ensembl_id": symbol_to_id[symbol], "gene_class": gene_class}
        for symbol, gene_class in gene_class_by_symbol.items()
        if symbol in symbol_to_id
    ]
    out_df = pd.DataFrame(rows).drop_duplicates(subset="ensembl_id").sort_values("ensembl_id").reset_index(drop=True)

    report = {
        "n_vogelstein_genes": len(vogelstein_by_symbol),
        "n_bailey_resolved_genes": len(bailey_by_symbol),
        "n_bailey_unresolved_cross_tissue_split": len(bailey_unresolved_split),
        "bailey_unresolved_cross_tissue_split_sample": sorted(bailey_unresolved_split)[:10],
        "n_source_disagreements_bailey_wins": len(disagreements),
        "source_disagreements_sample": disagreements[:10],
        "n_union_genes_before_ensembl_resolution": len(gene_class_by_symbol),
        "n_symbol_unresolved": len(unresolved),
        "n_symbol_ambiguous": len(ambiguous),
        "n_final_rows": len(out_df),
    }
    return out_df, report


def main():
    vogelstein_df = pd.read_csv(VOGELSTEIN_PATH)
    bailey_df = pd.read_csv(BAILEY_PATH)
    gene_reference_df = pd.read_csv(GENE_REFERENCE_PATH)

    out_df, report = build_gene_role_reference(vogelstein_df, bailey_df, gene_reference_df)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_PATH, index=False)
    print(f"Wrote {OUT_PATH} ({len(out_df):,} rows)")
    for key, value in report.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
