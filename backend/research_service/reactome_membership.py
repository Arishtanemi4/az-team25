"""Membership-only Reactome fallback (extensions/plans/CONTRACTS.md C4) -- used only when the
full built graph (rag/knowledge_graph.py's Parquet outputs) is absent but the raw Reactome
source files are present.

Parses data/external/reactome/raw/{Ensembl2Reactome.txt, ReactomePathways.txt} directly, using
the exact same column schema rag/build_knowledge_graph.py::load_reactome uses -- reimplemented
here, not imported, since that module is the frozen builder and C4 is explicit: "imports must
not invoke builder" and "do not run frozen builders that write frozen outputs." The schema
itself (Reactome's own published tab-separated, headerless format) is not proprietary code; only
the module that writes frozen outputs is off-limits to import.

Explicitly partial: no STRING, no BioGRID, no fusion/co-dependency edges -- pathway membership
only. Any derived cache this module writes lands under backend/research_service/runtime/, never inside
data/external/ -- extending the baseline's own outputs is not this package's job.
"""

import os
from pathlib import Path

import pandas as pd

_DEFAULT_CACHE_DIR = Path(__file__).resolve().parent / "runtime" / "reactome_membership_cache"


_LFS_POINTER_SIGNATURE = "version https://git-lfs.github.com/spec/v1"


def _is_lfs_pointer(path):
    """A Git LFS-tracked file that was never `git lfs pull`-ed still exists on disk as a small
    text stub (this worktree's own data/external/reactome/raw/ files are exactly this -- V6-7,
    docs/plan/STATUS.md item 60/61, materialised only a scoped subset elsewhere). `.exists()`
    alone cannot tell a real 183MB Reactome file from its 3-line pointer stub; this can, cheaply
    (read only the first line, never the whole file)."""
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            return handle.readline().strip() == _LFS_POINTER_SIGNATURE
    except OSError:
        return False


def _resolve_reactome_dir(reactome_dir=None):
    if reactome_dir is not None:
        return Path(reactome_dir)
    env_value = os.environ.get("REACTOME_RAW_DIR")
    if env_value:
        return Path(env_value)
    return Path(__file__).resolve().parents[2] / "data" / "external" / "reactome" / "raw"


def load_gene_pathway_membership(reactome_dir=None, cache_dir=None, use_cache=True):
    """Returns (gene_pathway_df, pathway_meta_df) -- same two columns/shape
    rag/build_knowledge_graph.py::load_reactome's first two return values have, parsed
    independently from the raw files. Filters to 'Homo sapiens' immediately, since most of each
    raw file's rows are other organisms, and requires `gene_id` to start with 'ENSG' (Reactome's
    own gene column is already Ensembl -- this excludes any other identifier namespace the raw
    file might carry).

    Caches the parsed result as Parquet under backend/research_service/runtime/ (never inside data/external/)
    so a repeated query against the same raw files does not re-parse them from scratch --
    invalidated whenever either raw file's mtime is newer than the cache's own.
    """
    reactome_dir = _resolve_reactome_dir(reactome_dir)
    e2r_path = reactome_dir / "Ensembl2Reactome.txt"
    pathways_path = reactome_dir / "ReactomePathways.txt"
    if not e2r_path.exists() or not pathways_path.exists():
        raise FileNotFoundError(f"Raw Reactome files not found under {reactome_dir}")
    if _is_lfs_pointer(e2r_path) or _is_lfs_pointer(pathways_path):
        raise FileNotFoundError(
            f"Raw Reactome files under {reactome_dir} are unpulled Git LFS pointer stubs, not "
            f"the real data -- run `git lfs pull` for data/external/reactome/raw/ before this "
            f"fallback can be used."
        )

    cache_dir = Path(cache_dir) if cache_dir is not None else _DEFAULT_CACHE_DIR
    gene_pathway_cache = cache_dir / "gene_pathway.parquet"
    pathway_meta_cache = cache_dir / "pathway_meta.parquet"
    source_mtime = max(e2r_path.stat().st_mtime, pathways_path.stat().st_mtime)

    if use_cache and gene_pathway_cache.exists() and pathway_meta_cache.exists():
        if gene_pathway_cache.stat().st_mtime >= source_mtime:
            return pd.read_parquet(gene_pathway_cache), pd.read_parquet(pathway_meta_cache)

    e2r = pd.read_csv(
        e2r_path,
        sep="\t",
        header=None,
        names=["gene_id", "reactome_id", "url", "pathway_name", "evidence_code", "species"],
        usecols=["gene_id", "reactome_id", "species"],
    )
    human_e2r = e2r[(e2r["species"] == "Homo sapiens") & e2r["gene_id"].str.startswith("ENSG")]
    gene_pathway_df = human_e2r[["gene_id", "reactome_id"]].drop_duplicates().reset_index(drop=True)

    pathways = pd.read_csv(
        pathways_path, sep="\t", header=None, names=["reactome_id", "pathway_name", "species"]
    )
    pathway_meta_df = (
        pathways[pathways["species"] == "Homo sapiens"][["reactome_id", "pathway_name"]]
        .reset_index(drop=True)
    )

    if use_cache:
        cache_dir.mkdir(parents=True, exist_ok=True)
        gene_pathway_df.to_parquet(gene_pathway_cache, engine="pyarrow")
        pathway_meta_df.to_parquet(pathway_meta_cache, engine="pyarrow")

    return gene_pathway_df, pathway_meta_df


def gene_pathway_memberships(ensembl_ids, reactome_dir=None, cache_dir=None, use_cache=True):
    """`{ensembl_id: [{"reactome_id", "pathway_name"}, ...]}` for exactly the requested genes,
    in stable `reactome_id` order. A gene simply absent from the membership table gets an empty
    list -- "checked, none found" -- distinct from the raw files being absent entirely, which
    raises FileNotFoundError before this function can return anything misleading."""
    gene_pathway_df, pathway_meta_df = load_gene_pathway_membership(
        reactome_dir, cache_dir, use_cache
    )
    pathway_names = pathway_meta_df.set_index("reactome_id")["pathway_name"].to_dict()

    result = {}
    for ensembl_id in ensembl_ids:
        rows = gene_pathway_df[gene_pathway_df["gene_id"] == ensembl_id].sort_values("reactome_id")
        result[ensembl_id] = [
            {"reactome_id": row.reactome_id, "pathway_name": pathway_names.get(row.reactome_id)}
            for row in rows.itertuples(index=False)
        ]
    return result
