import json
from pathlib import Path

import pandas as pd

import joins

REPO_ROOT = Path(__file__).resolve().parent.parent
RAW = REPO_ROOT / "data" / "raw"
AUGMENTED = REPO_ROOT / "data" / "augmented"
PROCESSED = REPO_ROOT / "data" / "processed"
OUTPUT_PATH = REPO_ROOT / "preprocessing" / "resources" / "unresolved_gene_symbols.json"

PATHS = {
    "protein": RAW / "gene_expression" / "4_Harmonized_MS_CCLE_Gygi_subsetted.csv",
    "dependency": AUGMENTED / "gene_properties" / "15_CRISPRGeneEffect.csv",
    "copy_number": AUGMENTED / "gene_properties" / "16_OmicsCNGeneWGS.csv",
}


def _collision_ensembl_ids(header_cols, gene_reference, symbol_position):
    _, report = joins.build_symbol_bridge(header_cols, gene_reference, symbol_position=symbol_position)
    ensembl_ids = set()
    for candidates in report["collision_candidates"].values():
        ensembl_ids.update(candidates)
    return sorted(ensembl_ids), report


def main():
    gene_reference = pd.read_csv(PROCESSED / "gene_reference.csv")

    protein_header = pd.read_csv(PATHS["protein"], nrows=0)
    protein_cols = [c for c in protein_header.columns if c != protein_header.columns[0]]
    protein_ids, protein_report = _collision_ensembl_ids(protein_cols, gene_reference, "inside")

    dep_header = pd.read_csv(PATHS["dependency"], nrows=0)
    dep_cols = [c for c in dep_header.columns if c != dep_header.columns[0]]
    dependency_ids, dependency_report = _collision_ensembl_ids(dep_cols, gene_reference, "before")

    cn_header = pd.read_csv(PATHS["copy_number"], nrows=0)
    cn_cols = [c for c in cn_header.columns if c not in
               ("SequencingID", "ModelConditionID", "ModelID", "IsDefaultEntryForMC", "IsDefaultEntryForModel")]
    copy_number_ids, copy_number_report = _collision_ensembl_ids(cn_cols, gene_reference, "before")

    output = {
        "protein": protein_ids,
        "dependency": dependency_ids,
        "copy_number": copy_number_ids,
        "_counts": {
            "protein": {"n_collisions": protein_report["n_collisions"], "n_no_match_not_included": protein_report["n_no_match"]},
            "dependency": {"n_collisions": dependency_report["n_collisions"], "n_no_match_not_included": dependency_report["n_no_match"]},
            "copy_number": {"n_collisions": copy_number_report["n_collisions"], "n_no_match_not_included": copy_number_report["n_no_match"]},
        },
    }
    OUTPUT_PATH.write_text(json.dumps(output, indent=2) + "\n")
    print(f"Wrote {OUTPUT_PATH} -- "
          f"protein {len(protein_ids)}, dependency {len(dependency_ids)}, copy_number {len(copy_number_ids)} "
          "collision-attributable genes (no-match genes intentionally excluded, see module docstring).")


if __name__ == "__main__":
    main()
