"""Fixed-cohort quantitative view extraction (CONTRACTS.md C5).

This module prepares the explicit, small matrices that S07 will compare.  It does not calculate
an alternatives score, change suitability, read a genome-wide matrix, or import AVI's separate
model catalogue.  Every percentile is calculated from the complete same-lineage processed spine
recorded in the snapshot, before top-k display, eligibility, or candidate filtering is considered.
"""

import json
import math
import sys
from pathlib import Path

import pandas as pd

_RESEARCH_DIR = str(Path(__file__).resolve().parent)
if _RESEARCH_DIR not in sys.path:
    sys.path.insert(0, _RESEARCH_DIR)

import serialization  # noqa: E402 -- extension public serialization boundary

POLICY_PATH = Path(__file__).with_name("similarity_policy.json")
PROFILE_SCHEMA_VERSION = "fixed-cohort-profiles-v1"


def load_similarity_policy(path=None):
    """Load and validate the single C5 policy source.

    The values are deliberately data, rather than constants repeated across S06/S07.  Keeping
    the AVI reference SHA here records inspiration without importing its incompatible loaders.
    """
    with open(Path(path) if path is not None else POLICY_PATH, encoding="utf-8") as handle:
        policy = json.load(handle)

    required = {
        "schema_version", "method_version", "avi_source_sha", "minimum_reference_n",
        "small_cohort_warning_n", "minimum_joint_fraction", "minimum_available_views", "top_n",
        "views",
    }
    missing = sorted(required - set(policy))
    if missing:
        raise ValueError(f"Similarity policy is missing required keys: {missing}")
    if len(policy["views"]) != 4:
        raise ValueError("C5 requires exactly four quantitative similarity views")
    if len({view["view_id"] for view in policy["views"]}) != 4:
        raise ValueError("Similarity policy view_id values must be unique")
    if any(view["weight"] <= 0 for view in policy["views"]):
        raise ValueError("Similarity policy weights must be positive")
    return policy


def _as_metadata_frame(snapshot):
    """Return the pinned model spine from the immutable snapshot, with strict identifier checks."""
    models = snapshot.get("models")
    if isinstance(models, pd.DataFrame):
        frame = models.copy()
    elif isinstance(models, list):
        frame = pd.DataFrame(models)
    else:
        raise ValueError("Snapshot models must be a list of pinned processed-spine metadata rows")

    required = {"ModelID", "lineage"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Snapshot model metadata is missing required columns: {missing}")
    if frame["ModelID"].isna().any() or frame["ModelID"].duplicated().any():
        raise ValueError("Snapshot model metadata must have one non-null row per ModelID")
    return frame


def _known_lineage(value):
    """Unknown lineage cannot be inferred from a model name or a query filter."""
    if value is None or pd.isna(value):
        return None
    value = str(value).strip()
    return value if value and value.lower() not in {"unknown", "nan", "none"} else None


def _human_reference_rows(metadata):
    """The processed spine is human-only; respect an explicit species column if one is supplied."""
    if "species" not in metadata.columns:
        return metadata.copy()
    species = metadata["species"].fillna("").astype(str).str.strip().str.lower()
    return metadata.loc[species.isin({"homo sapiens", "human"})].copy()


def _query_gene_roles(snapshot):
    """Deduplicate targets and exclusions in stable request order for C5's fixed denominator."""
    query = snapshot.get("query") or snapshot.get("native_result", {}).get("query") or {}
    genes = []
    roles = {}
    for role, field in (("inclusion", "inclusion_genes"), ("exclusion", "exclusion_genes")):
        for gene_id in query.get(field, []) or []:
            if gene_id not in roles:
                genes.append(gene_id)
                roles[gene_id] = role
    if not genes:
        raise ValueError("Snapshot query has no inclusion or exclusion genes for profile extraction")
    return genes, roles


def _eligible_native_ids(native_result):
    """Use only the complete native positive-eligible set, never the visible top-k alone."""
    result = []
    for bucket in ("ranked_cell_lines", "ranked_beyond_top_n"):
        for row in native_result.get(bucket, []) or []:
            score = row.get("D")
            tier = row.get("confidence_tier")
            veto = row.get("veto")
            if (
                isinstance(score, (int, float))
                and math.isfinite(score)
                and score > 0
                and tier in {"High", "Moderate"}
                and not veto
                and row.get("model_id") not in result
            ):
                result.append(row["model_id"])
    return result


def _anchor_pan_essential_flags(native_result, anchor_id, genes):
    """Retain native dependency pan-essential flags without using them to remove observations.

    Pan-essentialness is a suitability/explanation distinction from the frozen scorer, not a
    missing-data state.  S06 therefore carries it beside the dependency percentile table so S07
    can disclose it while still comparing the measured dependency values C5 permits.
    """
    for bucket in (
        "ranked_cell_lines", "ranked_beyond_top_n", "low_confidence_lines",
        "insufficient_evidence_lines", "disqualified_lines",
    ):
        for line in native_result.get(bucket, []) or []:
            if line.get("model_id") != anchor_id:
                continue
            flags = {gene_id: None for gene_id in genes}
            for gene in line.get("per_gene", []) or []:
                gene_id = gene.get("ensembl_id")
                if gene_id in flags:
                    flags[gene_id] = gene.get("layers", {}).get("dependency", {}).get("pan_essential")
            return flags
    return {gene_id: None for gene_id in genes}


def _source_state(row, view):
    """Return (state, numeric value) without turning any missingness into zero."""
    value_column = view["value_column"]
    if not row:
        return "not_observed", None
    if value_column not in row:
        return "missing_value_column", None
    if view["view_id"] == "protein":
        detected = row.get("detected")
        if detected is False or detected == 0:
            return "non_detected", None
        if detected is None or pd.isna(detected):
            return "detection_unknown", None
        if detected is not True and detected != 1:
            return "detection_unknown", None
    value = row.get(value_column)
    if value is None or pd.isna(value):
        return "not_observed", None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "non_numeric", None
    if not math.isfinite(value):
        return "non_finite", None
    return "observed", value


def _view_profile(frame, view, cohort_ids, genes, policy):
    """Build a per-view percentile table for only the requested genes and fixed cohort."""
    value_column = view["value_column"]
    required_columns = {"ModelID", "ensembl_id", value_column}
    if view["view_id"] == "protein":
        required_columns.add("detected")
    if not required_columns.issubset(frame.columns):
        missing = sorted(required_columns - set(frame.columns))
        return {
            "view_id": view["view_id"], "source": view["source"], "weight": view["weight"],
            "processed_layer": view["processed_layer"], "value_column": value_column,
            "status": "unavailable", "reason": f"Processed layer lacks required columns: {missing}",
            "gene_metadata": [], "percentile_table": [], "source_state_exclusions": {},
        }

    row_by_key = {
        (row.ModelID, row.ensembl_id): row._asdict()
        for row in frame.itertuples(index=False)
    }
    percentile_rows = []
    metadata = []
    exclusions = {}

    for gene_id in genes:
        observed = []
        raw_rows = []
        for model_id in cohort_ids:
            state, value = _source_state(row_by_key.get((model_id, gene_id), {}), view)
            raw_rows.append({"model_id": model_id, "gene_id": gene_id, "raw_value": value, "state": state})
            if state == "observed":
                observed.append(value)
            else:
                exclusions[state] = exclusions.get(state, 0) + 1

        observed_n = len(observed)
        distinct_n = len(set(observed))
        if observed_n < policy["minimum_reference_n"]:
            qualifying, reason = False, "observed_n_below_minimum"
        elif distinct_n < 2:
            qualifying, reason = False, "fewer_than_two_distinct_values"
        else:
            qualifying, reason = True, None

        ranks = pd.Series(observed, dtype=float).rank(method="average") if observed else pd.Series(dtype=float)
        observed_ranks = iter(ranks.tolist())
        ranked_values = []
        for row in raw_rows:
            if row["state"] == "observed":
                average_rank = float(next(observed_ranks))
                row["average_rank"] = average_rank
                row["percentile"] = average_rank / observed_n if qualifying else None
                ranked_values.append((row["raw_value"], average_rank))
            else:
                row["average_rank"] = None
                row["percentile"] = None
            percentile_rows.append(row)

        tie_counts = {}
        for value, _ in ranked_values:
            tie_counts[value] = tie_counts.get(value, 0) + 1
        metadata.append({
            "gene_id": gene_id,
            "observed_n": observed_n,
            "distinct_observed_n": distinct_n,
            "qualifying": qualifying,
            "non_discriminating_reason": reason,
            "tie_groups": [
                {"value": value, "count": count}
                for value, count in sorted(tie_counts.items()) if count > 1
            ],
            "small_cohort_warning": observed_n < policy["small_cohort_warning_n"],
        })

    return {
        "view_id": view["view_id"], "source": view["source"], "weight": view["weight"],
        "processed_layer": view["processed_layer"], "value_column": value_column,
        "status": "ok", "reason": None, "gene_metadata": metadata,
        "percentile_table": percentile_rows, "source_state_exclusions": exclusions,
    }


def extract_fixed_cohort_profiles(snapshot, anchor_id, adapter, data_dir=None, policy_path=None):
    """Extract C5 percentile views for one anchor without calculating a similarity score.

    The returned object is deliberately self-contained: S07 can calculate pairwise overlap from
    it without opening another data source or changing the immutable snapshot.  If the anchor's
    lineage is unknown, the function returns an explicit unavailable result rather than guessing.
    """
    policy = load_similarity_policy(policy_path)
    metadata = _as_metadata_frame(snapshot)
    anchor_rows = metadata.loc[metadata["ModelID"].eq(anchor_id)]
    if anchor_rows.empty:
        raise ValueError(f"Anchor model_id {anchor_id!r} is absent from snapshot metadata")
    lineage = _known_lineage(anchor_rows.iloc[0]["lineage"])
    genes, roles = _query_gene_roles(snapshot)
    base = {
        "schema_version": PROFILE_SCHEMA_VERSION,
        "query_id": snapshot["query_id"],
        "anchor_id": anchor_id,
        "query_genes": genes,
        "gene_roles": roles,
        "policy": policy,
    }
    if lineage is None:
        return {
            **base,
            "status": "unavailable",
            "reason": "Anchor lineage is unknown; a fixed same-lineage reference cohort cannot be defined.",
            "anchor_lineage": None,
            "reference_cohort": {"model_ids": [], "count": 0, "hash": None},
            "candidate_universe": {"eligible_model_ids": [], "candidate_model_ids": [], "count": 0},
            "views": {},
            "limitations": ["No lineage was inferred from the model name or query filter."],
        }

    human_models = _human_reference_rows(metadata)
    cohort_ids = human_models.loc[human_models["lineage"].eq(lineage), "ModelID"].tolist()
    cohort_ids = list(dict.fromkeys(cohort_ids))
    if anchor_id not in cohort_ids:
        return {
            **base,
            "status": "unavailable",
            "reason": "Anchor is not a human member of its declared lineage reference cohort.",
            "anchor_lineage": lineage,
            "reference_cohort": {"model_ids": cohort_ids, "count": len(cohort_ids), "hash": None},
            "candidate_universe": {"eligible_model_ids": [], "candidate_model_ids": [], "count": 0},
            "views": {},
            "limitations": [],
        }

    eligible_ids = _eligible_native_ids(snapshot["native_result"])
    cohort_id_set = set(cohort_ids)
    candidate_ids = [model_id for model_id in eligible_ids if model_id in cohort_id_set and model_id != anchor_id]
    dependency_pan_essential = _anchor_pan_essential_flags(snapshot["native_result"], anchor_id, genes)
    views = {}
    limitations = []
    for view in policy["views"]:
        try:
            frame = adapter.read_measurements(
                cohort_ids, genes, [view["processed_layer"]], data_dir=data_dir
            )[view["processed_layer"]]
            view_output = _view_profile(frame, view, cohort_ids, genes, policy)
        except FileNotFoundError as exc:
            view_output = {
                "view_id": view["view_id"], "source": view["source"], "weight": view["weight"],
                "processed_layer": view["processed_layer"], "value_column": view["value_column"],
                "status": "unavailable", "reason": str(exc), "gene_metadata": [],
                "percentile_table": [], "source_state_exclusions": {},
            }
        except Exception as exc:  # noqa: BLE001 -- report source/read failure without inventing values
            view_output = {
                "view_id": view["view_id"], "source": view["source"], "weight": view["weight"],
                "processed_layer": view["processed_layer"], "value_column": view["value_column"],
                "status": "error", "reason": str(exc), "gene_metadata": [],
                "percentile_table": [], "source_state_exclusions": {},
            }
        if view_output["status"] != "ok":
            limitations.append(f"{view['view_id']} view unavailable: {view_output['reason']}")
        if view["view_id"] == "dependency":
            view_output["anchor_pan_essential_flags"] = dependency_pan_essential
        views[view["view_id"]] = view_output

    if len(genes) == 1:
        limitations.append("One-gene query: any later query-profile similarity is intentionally narrow in scope.")
    return {
        **base,
        "status": "ok",
        "reason": None,
        "anchor_lineage": lineage,
        "reference_cohort": {
            "model_ids": cohort_ids,
            "count": len(cohort_ids),
            "hash": serialization.content_hash({"lineage": lineage, "model_ids": cohort_ids}),
        },
        "candidate_universe": {
            "eligible_model_ids": eligible_ids,
            "candidate_model_ids": candidate_ids,
            "count": len(candidate_ids),
        },
        "views": views,
        "limitations": limitations,
    }
