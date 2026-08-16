"""Transparent C5 query-profile alternatives, calculated only from S06 profile extraction.

The scalar produced here is a narrow, lineage-relative query-profile similarity heuristic.  It
does not rescore a model, modify native desirability or confidence, call native RNA Spearman, or
claim that two cell lines are experimentally interchangeable.
"""

import copy
import math

ALTERNATIVES_SCHEMA_VERSION = "query-profile-alternatives-v1"


def _native_eligible_lines(native_result):
    """Return the complete positive native candidate set, including the untruncated tail.

    S06 already makes this partition, but this second defensive check prevents a hand-edited
    profile object from reintroducing a Low, vetoed, zero-score, or insufficient line in S07.
    """
    lines = {}
    for bucket in ("ranked_cell_lines", "ranked_beyond_top_n"):
        for line in native_result.get(bucket, []) or []:
            score = line.get("D")
            if (
                isinstance(score, (int, float))
                and math.isfinite(score)
                and score > 0
                and line.get("confidence_tier") in {"High", "Moderate"}
                and not line.get("veto")
            ):
                lines.setdefault(line.get("model_id"), line)
    return lines


def _native_anchor_line(native_result, anchor_id):
    """Find the anchor in any native partition so legacy RNA-only similarities stay visible."""
    for bucket in (
        "ranked_cell_lines", "ranked_beyond_top_n", "low_confidence_lines",
        "insufficient_evidence_lines", "disqualified_lines",
    ):
        for line in native_result.get(bucket, []) or []:
            if line.get("model_id") == anchor_id:
                return line
    return None


def _empty_result(snapshot, anchor_id, profiles, policy, status, reason, considered=0):
    """Keep every non-success path explicit and schema-compatible for a research UI/export."""
    return {
        "schema_version": ALTERNATIVES_SCHEMA_VERSION,
        "query_id": snapshot.get("query_id"),
        "anchor_id": anchor_id,
        "method": {
            "version": policy.get("method_version"),
            "label": "query-profile similarity",
            "is_probability": False,
            "boundary": "This heuristic does not change native suitability or prove experimental interchangeability.",
        },
        "reference_cohort": copy.deepcopy(profiles.get("reference_cohort", {})),
        "candidate_counts": {"considered": considered, "qualifying": 0, "returned": 0},
        "status": status,
        "reason": reason,
        "top_alternatives": [],
        "legacy_native_similar_lines": [],
        "limitations": copy.deepcopy(profiles.get("limitations", [])),
    }


def _view_rows_by_model(view):
    """Index an S06 long percentile table as model -> gene -> row for one small view."""
    indexed = {}
    for row in view.get("percentile_table", []):
        indexed.setdefault(row["model_id"], {})[row["gene_id"]] = row
    return indexed


def _small_cohort_caveats(profiles):
    """Carry S06's n<30 warnings forward as caveats, never as a significance calculation."""
    caveats = []
    for view_id, view in profiles.get("views", {}).items():
        for gene in view.get("gene_metadata", []):
            if gene.get("small_cohort_warning"):
                caveats.append({
                    "view_id": view_id,
                    "gene_id": gene.get("gene_id"),
                    "observed_n": gene.get("observed_n"),
                    "message": "Reference cohort has fewer than 30 observed values for this gene/view.",
                })
    return caveats


def _score_view(view, anchor_id, candidate_id, query_genes, minimum_joint_fraction):
    """Apply C5's one-view percentile-distance formula without imputing absent values."""
    base = {
        "weight": view.get("weight"),
        "score": None,
        "joint_gene_count": 0,
        "joint_gene_fraction": 0.0,
        "qualifying_gene_count": 0,
        "missing_gene_ids": list(query_genes),
        "status": "unavailable",
        "reason": None,
    }
    if view.get("status") != "ok":
        base["reason"] = view.get("reason") or "source_unavailable"
        return base

    qualifying_genes = {
        item["gene_id"] for item in view.get("gene_metadata", []) if item.get("qualifying")
    }
    indexed = _view_rows_by_model(view)
    anchor_rows = indexed.get(anchor_id, {})
    candidate_rows = indexed.get(candidate_id, {})
    joint = []
    missing = []
    for gene_id in query_genes:
        if gene_id not in qualifying_genes:
            missing.append(gene_id)
            continue
        anchor_value = anchor_rows.get(gene_id, {}).get("percentile")
        candidate_value = candidate_rows.get(gene_id, {}).get("percentile")
        if anchor_value is None or candidate_value is None:
            missing.append(gene_id)
            continue
        joint.append(abs(float(anchor_value) - float(candidate_value)))

    joint_count = len(joint)
    joint_fraction = joint_count / len(query_genes)
    base.update({
        "joint_gene_count": joint_count,
        "joint_gene_fraction": joint_fraction,
        "qualifying_gene_count": len(qualifying_genes),
        "missing_gene_ids": missing,
    })
    if joint_count == 0:
        base["reason"] = "no_joint_qualifying_genes"
        return base
    if joint_fraction < minimum_joint_fraction:
        base["reason"] = "joint_gene_fraction_below_minimum"
        return base

    base["score"] = max(0.0, min(1.0, 1.0 - sum(joint) / joint_count))
    base["status"] = "available"
    return base


def _leave_one_view_out(view_scores):
    """Report the sensitivity range after omitting each available view in turn."""
    available = {view_id: item for view_id, item in view_scores.items() if item["status"] == "available"}
    by_omitted_view = {}
    for omitted_view in available:
        retained = [item for view_id, item in available.items() if view_id != omitted_view]
        total_weight = sum(item["weight"] for item in retained)
        by_omitted_view[omitted_view] = (
            sum(item["weight"] * item["score"] for item in retained) / total_weight
            if total_weight > 0 else None
        )
    values = [value for value in by_omitted_view.values() if value is not None]
    return {
        "minimum": min(values) if values else None,
        "maximum": max(values) if values else None,
        "by_omitted_view": by_omitted_view,
    }


def find_alternatives(snapshot, anchor_id, profiles, policy):
    """Return up to C5's top three transparent query-profile alternatives for an anchor.

    `profiles` must come from S06 for the same immutable snapshot and anchor.  The function
    performs no I/O and leaves both inputs unmodified, so toggling alternatives cannot mutate the
    frozen native result or cause a second score calculation.
    """
    if profiles.get("query_id") != snapshot.get("query_id"):
        raise ValueError("Profiles query_id does not match the supplied snapshot")
    if profiles.get("anchor_id") != anchor_id:
        raise ValueError("Profiles anchor_id does not match the requested anchor")
    if profiles.get("policy") != policy:
        raise ValueError("Profiles were extracted with a different similarity policy")
    if profiles.get("status") != "ok":
        return _empty_result(
            snapshot, anchor_id, profiles, policy, "unavailable",
            profiles.get("reason") or "Fixed-cohort profiles are unavailable.",
        )

    anchor_line = _native_anchor_line(snapshot.get("native_result", {}), anchor_id)
    if anchor_line is None:
        return _empty_result(
            snapshot, anchor_id, profiles, policy, "unavailable",
            "Anchor is absent from every native result partition.",
        )

    native_eligible = _native_eligible_lines(snapshot.get("native_result", {}))
    permitted_by_profiles = set(profiles.get("candidate_universe", {}).get("candidate_model_ids", []))
    cohort_ids = set(profiles.get("reference_cohort", {}).get("model_ids", []))
    candidates = [
        candidate_id for candidate_id in profiles.get("candidate_universe", {}).get("candidate_model_ids", [])
        if candidate_id != anchor_id
        and candidate_id in native_eligible
        and candidate_id in cohort_ids
    ]
    # A duplicate list would not be a second candidate. Stable de-duplication preserves S06's
    # audit order before the mandatory S-then-ModelID final sort.
    candidates = list(dict.fromkeys(candidate_id for candidate_id in candidates if candidate_id in permitted_by_profiles))
    if not candidates:
        result = _empty_result(
            snapshot, anchor_id, profiles, policy, "no_candidates",
            "No same-lineage positive eligible native candidates remain after exclusions.",
        )
        result["legacy_native_similar_lines"] = copy.deepcopy(anchor_line.get("similar_lines", []))
        return result

    query_genes = profiles.get("query_genes", [])
    if not query_genes:
        return _empty_result(snapshot, anchor_id, profiles, policy, "unavailable", "Profiles contain no query genes.")
    total_policy_weight = sum(view["weight"] for view in policy["views"])
    small_cohort = _small_cohort_caveats(profiles)
    qualifying = []
    for candidate_id in candidates:
        view_scores = {
            view["view_id"]: _score_view(
                profiles.get("views", {}).get(view["view_id"], {}), anchor_id, candidate_id,
                query_genes, policy["minimum_joint_fraction"],
            )
            for view in policy["views"]
        }
        available = {view_id: item for view_id, item in view_scores.items() if item["status"] == "available"}
        missing_views = [view_id for view_id in view_scores if view_id not in available]
        available_weight = sum(item["weight"] for item in available.values())
        if len(available) < policy["minimum_available_views"]:
            continue
        similarity = sum(item["weight"] * item["score"] for item in available.values()) / available_weight
        qualifying.append({
            "model_id": candidate_id,
            "S": max(0.0, min(1.0, similarity)),
            "native_D": native_eligible[candidate_id].get("D"),
            "native_confidence_tier": native_eligible[candidate_id].get("confidence_tier"),
            "view_scores": view_scores,
            "missing_views": missing_views,
            "available_weight": available_weight,
            "available_weight_fraction": available_weight / total_policy_weight,
            "small_cohort_caveats": copy.deepcopy(small_cohort),
            "one_gene_caveat": len(query_genes) == 1,
            "leave_one_view_out_S_range": _leave_one_view_out(view_scores),
        })

    if not qualifying:
        result = _empty_result(
            snapshot, anchor_id, profiles, policy, "no_qualifying_candidates",
            "No candidate had the minimum multi-view overlap required by the predeclared policy.",
            considered=len(candidates),
        )
        result["legacy_native_similar_lines"] = copy.deepcopy(anchor_line.get("similar_lines", []))
        return result

    qualifying.sort(key=lambda item: (-item["S"], item["model_id"]))
    returned = qualifying[:policy["top_n"]]
    for rank, item in enumerate(returned, start=1):
        item["rank"] = rank
    return {
        "schema_version": ALTERNATIVES_SCHEMA_VERSION,
        "query_id": snapshot.get("query_id"),
        "anchor_id": anchor_id,
        "method": {
            "version": policy["method_version"],
            "label": "query-profile similarity",
            "is_probability": False,
            "boundary": "This heuristic does not change native suitability or prove experimental interchangeability.",
        },
        "reference_cohort": copy.deepcopy(profiles.get("reference_cohort", {})),
        "candidate_counts": {"considered": len(candidates), "qualifying": len(qualifying), "returned": len(returned)},
        "status": "ok",
        "reason": None,
        "top_alternatives": returned,
        "legacy_native_similar_lines": copy.deepcopy(anchor_line.get("similar_lines", [])),
        "limitations": copy.deepcopy(profiles.get("limitations", [])),
    }
