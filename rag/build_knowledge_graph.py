import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
REACTOME_DIR = REPO_ROOT / "data" / "external" / "reactome" / "raw"
STRING_DIR = REPO_ROOT / "data" / "external" / "string_v12" / "raw"
BIOGRID_ZIP = REPO_ROOT / "data" / "external" / "biogrid" / "raw" / "BIOGRID-ORGANISM-5.0.259.tab3.zip"
BIOGRID_HUMAN_MEMBER = "BIOGRID-ORGANISM-Homo_sapiens-5.0.259.tab3.txt"
GENE_REFERENCE_PATH = REPO_ROOT / "data" / "processed" / "gene_reference.csv"
FUSIONS_PATH = REPO_ROOT / "data" / "processed" / "fusions.csv"
DEPENDENCY_PATH = REPO_ROOT / "data" / "processed" / "dependency.csv"
OUT_DIR = REPO_ROOT / "data" / "external" / "knowledge_graph" / "derived"
NODES_OUT_PATH = OUT_DIR / "knowledge_graph_nodes.parquet"
EDGES_OUT_PATH = OUT_DIR / "knowledge_graph_edges.parquet"
MANIFEST_OUT_PATH = OUT_DIR / "knowledge_graph_build_manifest.json"
LEGACY_EDGES_PATH = OUT_DIR / "knowledge_graph_edges.csv"
LEGACY_NODES_PATH = OUT_DIR / "knowledge_graph_nodes.csv"

# STRING v12 documents 0.4/0.7/0.9 as its own medium/high/highest confidence bands -- 700 is
# STRING's "high confidence" cutoff, not an invented threshold (docs/plan/PARAMETERS.md SS13).
STRING_COMBINED_SCORE_THRESHOLD = 700
STRING_LINKS_CHUNKSIZE = 1_000_000  # matches the project's existing convention, PARAMETERS.md row 39

# Co-dependency sparsification -- never materialise the full ~18,398x18,398 correlation matrix
# (that class of memory pressure is what forced V6-6's batched RNA reshape). Chunk over genes,
# sparsify immediately, keep only the strongest partners.
CODEPENDENCY_TOP_K = 20
CODEPENDENCY_MIN_ABS_R = 0.2
CODEPENDENCY_MIN_N_SHARED = 30  # reuses PARAMETERS.md row 9's existing min_calibration_n floor
CODEPENDENCY_CHUNK_SIZE = 500


def resolve_symbols(symbols, gene_reference_df):
    counts = gene_reference_df.groupby("symbol")["ensembl_id"].nunique()
    unique_symbols = set(counts[counts == 1].index)
    symbol_to_id = (
        gene_reference_df[gene_reference_df["symbol"].isin(unique_symbols)]
        .set_index("symbol")["ensembl_id"]
        .to_dict()
    )

    resolved, unresolved, ambiguous = {}, [], []
    for symbol in symbols:
        if symbol in symbol_to_id:
            resolved[symbol] = symbol_to_id[symbol]
        elif symbol in counts.index:
            ambiguous.append(symbol)
        else:
            unresolved.append(symbol)
    return resolved, unresolved, ambiguous


def validate_reactome_ensembl_ids(ids, gene_reference_df):
    known = set(gene_reference_df["ensembl_id"])
    ids = set(ids)
    valid = ids & known
    not_in_gene_reference = sorted(ids - known)
    return valid, not_in_gene_reference


def load_reactome(reactome_dir=REACTOME_DIR):
    e2r = pd.read_csv(
        reactome_dir / "Ensembl2Reactome.txt",
        sep="\t",
        header=None,
        names=["gene_id", "reactome_id", "url", "pathway_name", "evidence_code", "species"],
        usecols=["gene_id", "reactome_id", "species"],
    )
    human_e2r = e2r[(e2r["species"] == "Homo sapiens") & e2r["gene_id"].str.startswith("ENSG")]
    gene_pathway_df = human_e2r[["gene_id", "reactome_id"]].drop_duplicates().reset_index(drop=True)

    pathways = pd.read_csv(
        reactome_dir / "ReactomePathways.txt",
        sep="\t",
        header=None,
        names=["reactome_id", "pathway_name", "species"],
    )
    pathway_meta_df = pathways[pathways["species"] == "Homo sapiens"][["reactome_id", "pathway_name"]].reset_index(drop=True)
    human_pathway_ids = set(pathway_meta_df["reactome_id"])

    relation = pd.read_csv(
        reactome_dir / "ReactomePathwaysRelation.txt",
        sep="\t",
        header=None,
        names=["parent_id", "child_id"],
    )
    hierarchy_df = relation[
        relation["parent_id"].isin(human_pathway_ids) & relation["child_id"].isin(human_pathway_ids)
    ].reset_index(drop=True)

    return gene_pathway_df, pathway_meta_df, hierarchy_df


def load_string(string_dir=STRING_DIR, min_combined_score=STRING_COMBINED_SCORE_THRESHOLD):
    info = pd.read_csv(string_dir / "9606.protein.info.v12.0.txt.gz", sep="\t", usecols=["#string_protein_id", "preferred_name"])
    protein_to_symbol = dict(zip(info["#string_protein_id"], info["preferred_name"]))

    kept_chunks = []
    n_rows_seen = 0
    for chunk in pd.read_csv(
        string_dir / "9606.protein.links.v12.0.txt.gz", sep=" ", chunksize=STRING_LINKS_CHUNKSIZE
    ):
        n_rows_seen += len(chunk)
        kept_chunks.append(chunk[chunk["combined_score"] >= min_combined_score])
    high_confidence = pd.concat(kept_chunks, ignore_index=True)

    high_confidence["symbol_a"] = high_confidence["protein1"].map(protein_to_symbol)
    high_confidence["symbol_b"] = high_confidence["protein2"].map(protein_to_symbol)
    # STRING lists both (A,B) and (B,A) -- collapse to one unordered row per pair. Vectorised
    # (np.where, not a row-wise .apply) -- a Python-level apply over this many rows is slow
    # enough to look like a hang.
    a, b = high_confidence["symbol_a"].to_numpy(), high_confidence["symbol_b"].to_numpy()
    lo, hi = np.where(a < b, a, b), np.where(a < b, b, a)
    high_confidence = high_confidence.assign(symbol_a=lo, symbol_b=hi).drop_duplicates(subset=["symbol_a", "symbol_b"])

    string_edges_df = high_confidence[["symbol_a", "symbol_b", "combined_score"]].reset_index(drop=True)
    return string_edges_df, n_rows_seen


def load_biogrid(biogrid_zip=BIOGRID_ZIP, member=BIOGRID_HUMAN_MEMBER):
    with zipfile.ZipFile(biogrid_zip) as zf:
        with zf.open(member) as fh:
            df = pd.read_csv(
                fh,
                sep="\t",
                usecols=[
                    "Official Symbol Interactor A",
                    "Official Symbol Interactor B",
                    "Experimental System Type",
                    "Organism ID Interactor A",
                    "Organism ID Interactor B",
                ],
            )
    both_human = (df["Organism ID Interactor A"] == 9606) & (df["Organism ID Interactor B"] == 9606)
    df = df[both_human].rename(
        columns={"Official Symbol Interactor A": "symbol_a", "Official Symbol Interactor B": "symbol_b"}
    )
    # Vectorised unordered-pair key (np.where, not a row-wise .apply -- BioGRID is ~1.4M rows,
    # where a Python-level apply is slow enough to look like a hang).
    a, b = df["symbol_a"].to_numpy(), df["symbol_b"].to_numpy()
    df = df.assign(symbol_a=np.where(a < b, a, b), symbol_b=np.where(a < b, b, a))
    grouped = (
        df.groupby(["symbol_a", "symbol_b"])["Experimental System Type"]
        .agg(experimental_system_types=lambda s: "|".join(sorted(set(s))), n_observations="size")
        .reset_index()
    )
    return grouped, len(df)


def build_fusion_partner_edges(fusions_path=FUSIONS_PATH):
    fusions_df = pd.read_csv(
        fusions_path, usecols=["ModelID", "ensembl_id", "partner_ensembl_id", "confidence_high", "in_frame"]
    )
    fusions_df = fusions_df.dropna(subset=["ensembl_id", "partner_ensembl_id"])
    # Vectorised unordered-pair key -- gene_a/gene_b are already deterministically sorted per row,
    # so the groupby below needs no separate re-sort step afterwards.
    a, b = fusions_df["ensembl_id"].to_numpy(), fusions_df["partner_ensembl_id"].to_numpy()
    fusions_df = fusions_df.assign(gene_a=np.where(a < b, a, b), gene_b=np.where(a < b, b, a))
    grouped = (
        fusions_df.groupby(["gene_a", "gene_b"])
        .agg(
            n_models_observed=("ModelID", "nunique"),
            any_high_confidence=("confidence_high", "any"),
            any_in_frame=("in_frame", "any"),
        )
        .reset_index()
    )
    return grouped


def compute_codependency_edges(
    dependency_path=DEPENDENCY_PATH,
    top_k=CODEPENDENCY_TOP_K,
    min_abs_r=CODEPENDENCY_MIN_ABS_R,
    min_n_shared=CODEPENDENCY_MIN_N_SHARED,
    chunk_size=CODEPENDENCY_CHUNK_SIZE,
):
    long_df = pd.read_csv(dependency_path, usecols=["ModelID", "ensembl_id", "dependency_score"])
    wide = long_df.pivot(index="ModelID", columns="ensembl_id", values="dependency_score")
    genes = wide.columns.to_numpy()
    values = wide.to_numpy(dtype=np.float32)  # NaN where a (model, gene) pair was never assayed
    mask = ~np.isnan(values)

    mean = np.nanmean(values, axis=0)
    std = np.nanstd(values, axis=0)
    standardized = (values - mean) / std
    standardized = np.nan_to_num(standardized, nan=0.0)  # missing entries contribute 0 to every dot product
    mask_f = mask.astype(np.float32)

    n_genes = len(genes)
    edge_rows = []
    for start in range(0, n_genes, chunk_size):
        end = min(start + chunk_size, n_genes)
        z_chunk = standardized[:, start:end]  # (n_models, chunk_width)
        mask_chunk = mask_f[:, start:end]

        numerator = z_chunk.T @ standardized  # (chunk_width, n_genes)
        n_shared = mask_chunk.T @ mask_f  # (chunk_width, n_genes), pairwise-present model count
        with np.errstate(divide="ignore", invalid="ignore"):
            r = numerator / np.maximum(n_shared - 1, 1)

        for local_i, global_i in enumerate(range(start, end)):
            row_r = r[local_i]
            row_n = n_shared[local_i]
            candidate = np.where((np.abs(row_r) >= min_abs_r) & (row_n >= min_n_shared))[0]
            candidate = candidate[candidate != global_i]  # drop self-pair
            if candidate.size == 0:
                continue
            top = candidate[np.argsort(-np.abs(row_r[candidate]))[:top_k]]
            for j in top:
                edge_rows.append((genes[global_i], genes[j], float(row_r[j]), int(row_n[j])))

    edges_df = pd.DataFrame(edge_rows, columns=["gene_a", "gene_b", "r", "n_shared"])
    # Both directions of a strong pair can each independently qualify for their own top-k list --
    # collapse to one undirected row per pair, keeping the higher |r| if the two computed values
    # differ. Vectorised (np.where, not a row-wise .apply).
    a, b = edges_df["gene_a"].to_numpy(), edges_df["gene_b"].to_numpy()
    edges_df = edges_df.assign(gene_a=np.where(a < b, a, b), gene_b=np.where(a < b, b, a), _abs_r=edges_df["r"].abs())
    edges_df = edges_df.sort_values("_abs_r", ascending=False).drop_duplicates(subset=["gene_a", "gene_b"])
    return edges_df[["gene_a", "gene_b", "r", "n_shared"]].reset_index(drop=True)


def build_nodes_table(gene_ensembl_ids, pathway_meta_df, gene_reference_df):
    symbol_by_ensembl = gene_reference_df.set_index("ensembl_id")["symbol"].to_dict()
    gene_rows = [
        {
            "node_id": f"gene:{eid}",
            "node_type": "gene",
            "source_db": "gene_reference",
            "display_name": symbol_by_ensembl.get(eid, eid),
            "ensembl_id": eid,
            "reactome_id": None,
            "species": "Homo sapiens",
        }
        for eid in sorted(set(gene_ensembl_ids))
    ]
    pathway_rows = [
        {
            "node_id": f"pathway:{row.reactome_id}",
            "node_type": "pathway",
            "source_db": "reactome",
            "display_name": row.pathway_name,
            "ensembl_id": None,
            "reactome_id": row.reactome_id,
            "species": "Homo sapiens",
        }
        for row in pathway_meta_df.itertuples()
    ]
    nodes_df = pd.DataFrame(gene_rows + pathway_rows)
    return nodes_df


def _gene_node(ensembl_id):
    return f"gene:{ensembl_id}"


def _pathway_node(reactome_id):
    return f"pathway:{reactome_id}"


def build_edges_table(
    reactome_gene_pathway_resolved,
    reactome_hierarchy_df,
    string_edges_resolved,
    biogrid_edges_resolved,
    fusion_edges_df,
    codependency_edges_df,
):
    frames = []

    gp = reactome_gene_pathway_resolved
    frames.append(
        pd.DataFrame(
            {
                "source_node_id": gp["gene_id"].map(_gene_node),
                "target_node_id": gp["reactome_id"].map(_pathway_node),
                "edge_type": "gene_in_pathway",
                "source_db": "reactome",
                "weight": np.nan,
                "directed": False,
                "evidence_detail": None,
                "n_source_observations": 1,
            }
        )
    )

    rh = reactome_hierarchy_df
    frames.append(
        pd.DataFrame(
            {
                "source_node_id": rh["parent_id"].map(_pathway_node),
                "target_node_id": rh["child_id"].map(_pathway_node),
                "edge_type": "pathway_parent_of",
                "source_db": "reactome",
                "weight": np.nan,
                "directed": True,
                "evidence_detail": None,
                "n_source_observations": 1,
            }
        )
    )

    se = string_edges_resolved
    frames.append(
        pd.DataFrame(
            {
                "source_node_id": se["ensembl_a"].map(_gene_node),
                "target_node_id": se["ensembl_b"].map(_gene_node),
                "edge_type": "protein_interaction",
                "source_db": "string_v12",
                "weight": se["combined_score"] / 1000.0,
                "directed": False,
                "evidence_detail": None,
                "n_source_observations": 1,
            }
        )
    )

    be = biogrid_edges_resolved
    frames.append(
        pd.DataFrame(
            {
                "source_node_id": be["ensembl_a"].map(_gene_node),
                "target_node_id": be["ensembl_b"].map(_gene_node),
                "edge_type": "curated_interaction",
                "source_db": "biogrid",
                "weight": np.nan,
                "directed": False,
                "evidence_detail": be["experimental_system_types"],
                "n_source_observations": be["n_observations"],
            }
        )
    )

    fe = fusion_edges_df
    frames.append(
        pd.DataFrame(
            {
                "source_node_id": fe["gene_a"].map(_gene_node),
                "target_node_id": fe["gene_b"].map(_gene_node),
                "edge_type": "fusion_partner",
                "source_db": "project_fusions",
                "weight": np.nan,
                "directed": False,
                "evidence_detail": fe.apply(
                    lambda r: f"high_confidence={r['any_high_confidence']};in_frame={r['any_in_frame']}", axis=1
                ),
                "n_source_observations": fe["n_models_observed"],
            }
        )
    )

    ce = codependency_edges_df
    frames.append(
        pd.DataFrame(
            {
                "source_node_id": ce["gene_a"].map(_gene_node),
                "target_node_id": ce["gene_b"].map(_gene_node),
                "edge_type": "codependency",
                "source_db": "project_dependency",
                "weight": ce["r"],
                "directed": False,
                "evidence_detail": None,
                "n_source_observations": ce["n_shared"],
            }
        )
    )

    edges_df = pd.concat(frames, ignore_index=True)
    return edges_df


def cross_check_against_legacy(edges_df, nodes_df):
    report = {"legacy_files_found": False}
    if LEGACY_EDGES_PATH.exists() and LEGACY_NODES_PATH.exists():
        try:
            legacy_edges = pd.read_csv(LEGACY_EDGES_PATH, low_memory=False)
            legacy_nodes = pd.read_csv(LEGACY_NODES_PATH, low_memory=False)
            report["legacy_files_found"] = True
            report["legacy_edge_count"] = len(legacy_edges)
            report["legacy_node_count"] = len(legacy_nodes)
            report["fresh_edge_count"] = len(edges_df)
            report["fresh_node_count"] = len(nodes_df)
        except Exception as exc:  # legacy file schema is not this script's responsibility
            report["legacy_read_error"] = str(exc)
    return report


def write_outputs(nodes_df, edges_df, manifest):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    nodes_df.to_parquet(NODES_OUT_PATH, index=False)
    edges_df.to_parquet(EDGES_OUT_PATH, index=False)
    with open(MANIFEST_OUT_PATH, "w") as fh:
        json.dump(manifest, fh, indent=2, default=str)


def main():
    gene_reference_df = pd.read_csv(GENE_REFERENCE_PATH)
    report = {"built_utc": datetime.now(timezone.utc).isoformat()}

    print("Loading Reactome (human-filtered)...")
    gene_pathway_df, pathway_meta_df, hierarchy_df = load_reactome()
    valid_ids, not_in_gene_reference = validate_reactome_ensembl_ids(gene_pathway_df["gene_id"], gene_reference_df)
    reactome_gene_pathway_resolved = gene_pathway_df[gene_pathway_df["gene_id"].isin(valid_ids)].reset_index(drop=True)
    report["reactome"] = {
        "human_gene_pathway_rows_raw": len(gene_pathway_df),
        "resolved_against_gene_reference": len(reactome_gene_pathway_resolved),
        "n_ensembl_ids_not_in_gene_reference": len(not_in_gene_reference),
        "n_pathways": len(pathway_meta_df),
        "n_hierarchy_edges": len(hierarchy_df),
    }

    print("Loading STRING v12 (combined_score >= %d)..." % STRING_COMBINED_SCORE_THRESHOLD)
    string_edges_df, string_rows_seen = load_string()
    string_symbols = set(string_edges_df["symbol_a"]) | set(string_edges_df["symbol_b"])
    string_resolved, string_unresolved, string_ambiguous = resolve_symbols(string_symbols, gene_reference_df)
    string_edges_df["ensembl_a"] = string_edges_df["symbol_a"].map(string_resolved)
    string_edges_df["ensembl_b"] = string_edges_df["symbol_b"].map(string_resolved)
    string_edges_resolved = string_edges_df.dropna(subset=["ensembl_a", "ensembl_b"]).reset_index(drop=True)
    report["string_v12"] = {
        "raw_link_rows_scanned": string_rows_seen,
        "high_confidence_pairs": len(string_edges_df),
        "resolved_edges": len(string_edges_resolved),
        "n_symbols_unresolved": len(string_unresolved),
        "n_symbols_ambiguous": len(string_ambiguous),
    }

    print("Loading BioGRID (Homo sapiens member)...")
    biogrid_edges_df, biogrid_rows_seen = load_biogrid()
    biogrid_symbols = set(biogrid_edges_df["symbol_a"]) | set(biogrid_edges_df["symbol_b"])
    biogrid_resolved, biogrid_unresolved, biogrid_ambiguous = resolve_symbols(biogrid_symbols, gene_reference_df)
    biogrid_edges_df["ensembl_a"] = biogrid_edges_df["symbol_a"].map(biogrid_resolved)
    biogrid_edges_df["ensembl_b"] = biogrid_edges_df["symbol_b"].map(biogrid_resolved)
    biogrid_edges_resolved = biogrid_edges_df.dropna(subset=["ensembl_a", "ensembl_b"]).reset_index(drop=True)
    report["biogrid"] = {
        "raw_interaction_rows": biogrid_rows_seen,
        "unique_pairs": len(biogrid_edges_df),
        "resolved_edges": len(biogrid_edges_resolved),
        "n_symbols_unresolved": len(biogrid_unresolved),
        "n_symbols_ambiguous": len(biogrid_ambiguous),
    }

    print("Building fusion-partner edges from data/processed/fusions.csv...")
    fusion_edges_df = build_fusion_partner_edges()
    report["fusion_partner"] = {"n_edges": len(fusion_edges_df)}

    print("Computing co-dependency edges from data/processed/dependency.csv (chunked, sparsified)...")
    codependency_edges_df = compute_codependency_edges()
    report["codependency"] = {
        "n_edges": len(codependency_edges_df),
        "top_k": CODEPENDENCY_TOP_K,
        "min_abs_r": CODEPENDENCY_MIN_ABS_R,
        "min_n_shared": CODEPENDENCY_MIN_N_SHARED,
    }

    print("Assembling combined edge table...")
    edges_df = build_edges_table(
        reactome_gene_pathway_resolved,
        hierarchy_df,
        string_edges_resolved,
        biogrid_edges_resolved,
        fusion_edges_df,
        codependency_edges_df,
    )

    all_gene_ensembl_ids = set()
    for series in [
        reactome_gene_pathway_resolved["gene_id"],
        string_edges_resolved["ensembl_a"],
        string_edges_resolved["ensembl_b"],
        biogrid_edges_resolved["ensembl_a"],
        biogrid_edges_resolved["ensembl_b"],
        fusion_edges_df["gene_a"],
        fusion_edges_df["gene_b"],
        codependency_edges_df["gene_a"],
        codependency_edges_df["gene_b"],
    ]:
        all_gene_ensembl_ids.update(series.tolist())

    print("Building nodes table...")
    nodes_df = build_nodes_table(all_gene_ensembl_ids, pathway_meta_df, gene_reference_df)

    report["combined"] = {
        "n_nodes": len(nodes_df),
        "n_gene_nodes": int((nodes_df["node_type"] == "gene").sum()),
        "n_pathway_nodes": int((nodes_df["node_type"] == "pathway").sum()),
        "n_edges": len(edges_df),
        "n_edges_by_type": edges_df["edge_type"].value_counts().to_dict(),
        "n_edges_by_source_db": edges_df["source_db"].value_counts().to_dict(),
    }

    print("Cross-checking against the legacy (non-authoritative) CSVs, if present...")
    report["legacy_cross_check"] = cross_check_against_legacy(edges_df, nodes_df)

    write_outputs(nodes_df, edges_df, report)

    print(f"\nWrote {NODES_OUT_PATH} ({len(nodes_df):,} rows)")
    print(f"Wrote {EDGES_OUT_PATH} ({len(edges_df):,} rows)")
    print(f"Wrote {MANIFEST_OUT_PATH}")
    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
