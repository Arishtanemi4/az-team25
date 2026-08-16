"""Deterministic measured biological context packets (extensions/plans/CONTRACTS.md C3).

`build_measured_context(snapshot, model_id, adapter)` assembles everything already computed by
the frozen scorer for one model, copied verbatim, plus a small amount of additive context read
fresh through the S03 adapter (raw mutation/fusion event rows, supplementary assay availability)
-- never a rescore, never a new interpretive claim. Works entirely offline: no graph, no LLM, and
degrades to "unavailable" rather than crashing when the adapter has no data to read (this frozen
worktree's own condition).
"""

import copy
import sys
from pathlib import Path

_RESEARCH_DIR = str(Path(__file__).resolve().parent)
if _RESEARCH_DIR not in sys.path:
    sys.path.insert(0, _RESEARCH_DIR)

import templates  # noqa: E402 -- sys.path must be set up first

PARTITION_BUCKETS = (
    "ranked_cell_lines",
    "ranked_beyond_top_n",
    "low_confidence_lines",
    "insufficient_evidence_lines",
    "disqualified_lines",
)

MAX_DISPLAYED_EVENTS = 20

# Layer name in the native per_gene["layers"] dict -> adapter.LAYER_FILES key. Only the two
# event-table layers need a fresh read here: every other layer's value/unit/d/state is already
# complete in the snapshot's own native_result, copied verbatim (C3), needing no additive read.
EVENT_LAYER_NAMES = {"mutation": "mutations", "fusion": "fusions"}


def _find_model_result(native_result, model_id):
    """Searches every partition bucket, not just `ranked_cell_lines` -- a model_id may be
    low-confidence, insufficient, or disqualified, and each of those is still a real result
    this project promises never to silently drop (scoring/_.md's Outputs table)."""
    for bucket in PARTITION_BUCKETS:
        for line in native_result.get(bucket, []):
            if line.get("model_id") == model_id:
                return line, bucket
    raise ValueError(
        f"model_id {model_id!r} is not present in any partition bucket of this snapshot "
        f"(checked: {', '.join(PARTITION_BUCKETS)})"
    )


def _warning_message_for_gene(warnings, ensembl_id):
    """Reuses a top-level warning's already-formatted message (e.g. the pan-essential note,
    which needs a fraction_strong value this module has no other access to) rather than
    reformatting the same template a second time from a different data path."""
    for warning in warnings:
        if warning.get("ensembl_id") == ensembl_id:
            return warning.get("message")
    return None


def _missing_evidence_caveats(missing_evidence, ensembl_id):
    """Every already-computed missing_evidence entry for this gene (scoring/explain.py's own
    public missing_evidence_report output, shipped verbatim in native_result) becomes one
    structured caveat -- no severity is recomputed here, only relabelled as {type, message}."""
    caveats = []
    for entry in missing_evidence:
        if entry.get("ensembl_id") != ensembl_id:
            continue
        severity = entry["severity"]
        layer = entry.get("layer")
        if severity == "gene_role_unknown":
            message = templates.GENE_ROLE_UNKNOWN_NOTE.format(layer=layer.replace("_", " ").capitalize())
        elif severity == "cannot_certify_absence":
            message = (
                "This is an exclusion criterion with no evidence in any measured layer -- "
                "absence cannot be certified, which is the dangerous case for an exclusion "
                "gene, not the safe one."
            )
        elif severity == "untested":
            message = "No evidence in any measured layer for this inclusion gene -- untested, not scored."
        elif severity == "measured_unresolvable":
            message = (
                f"{layer.replace('_', ' ').capitalize()} data for this gene was measured but "
                f"left out during identity resolution (its symbol collided with another gene's) "
                f"-- a false 'never looked' would be misleading; this is 'measured, unresolvable'."
            )
        else:
            # An unknown future severity must render explicitly, never silently count as
            # measured or be dropped (C3: "unknown future states must render explicitly").
            message = f"Unrecognised missing-evidence severity {severity!r} for layer {layer!r}."
        caveats.append({"type": severity, "layer": layer, "message": message})
    return caveats


def _layer_caveats(layer_name, layer_info):
    """Per-layer structural caveats derived directly from that layer's own state -- never a
    new claim, only C3's own stated scope limits attached to the layer they describe."""
    caveats = []
    state = layer_info.get("state")

    if state == "not_assayed":
        if layer_name == "fusion":
            caveats.append({"type": "not_assayed", "message": templates.FUSION_STATUS_UNKNOWN_NOTE})
        else:
            caveats.append({
                "type": "not_assayed",
                "message": (
                    f"{layer_name.replace('_', ' ').capitalize()} was not assayed for this gene "
                    f"on this line -- absence of evidence, not evidence of absence."
                ),
            })
    elif state == "non_detected":
        caveats.append({"type": "non_detected", "message": templates.PROTEIN_NON_DETECTED_NOTE})
    elif state == "measured_absent" and layer_name in ("mutation", "fusion"):
        caveats.append({
            "type": "measured_absent_not_confirmed_absence",
            "message": templates.MEASURED_ABSENT_NOT_CONFIRMED_ABSENCE_NOTE,
        })
    elif state == "measured" and layer_name == "copy_number":
        caveats.append({
            "type": "relative_not_absolute",
            "message": templates.COPY_NUMBER_RELATIVE_NOT_ABSOLUTE_NOTE,
        })
    elif state == "measured" and layer_name == "protein":
        caveats.append({"type": "mixed_scale", "message": templates.PROTEIN_MIXED_SCALE_NOTE})

    return caveats


def _read_event_rows(adapter_module, model_id, ensembl_id, layer_name, data_dir):
    """Additive-only: reads the raw event rows for one (model, gene, event-layer) through the
    S03 adapter. Never raises past this function -- a missing table in this frozen worktree, or
    any other adapter failure, becomes an explicit 'unavailable'/'error' status (C2's own
    component-status convention: ok / partial / unavailable / error, with a reason), not a
    crashed context build."""
    adapter_layer = EVENT_LAYER_NAMES[layer_name]
    try:
        result = adapter_module.read_measurements(
            [model_id], [ensembl_id], [adapter_layer], data_dir=data_dir
        )
    except FileNotFoundError as exc:
        return {"status": "unavailable", "reason": str(exc), "rows": [], "total_count": 0}
    except Exception as exc:  # noqa: BLE001 -- any other adapter failure is reported, not raised
        return {"status": "error", "reason": str(exc), "rows": [], "total_count": 0}

    frame = result[adapter_layer]
    total_count = len(frame)
    displayed = frame.head(MAX_DISPLAYED_EVENTS).to_dict(orient="records")
    payload = {"status": "ok", "reason": None, "rows": displayed, "total_count": total_count}
    if total_count > MAX_DISPLAYED_EVENTS:
        payload["truncation_note"] = templates.EVENTS_TRUNCATED_NOTE.format(
            displayed=MAX_DISPLAYED_EVENTS, total=total_count
        )
    return payload


def _read_supplementary_availability(adapter_module, model_id, data_dir):
    """Model-level supplementary assay availability (metabolomics/miRNA/genome_signatures) --
    reported as coverage/context, never as a new score input (C3, root _.md SS4/SS6). Each
    layer's availability is independent: one missing table must not hide another's real
    'available' answer. These three tables are not gene-keyed (no ensembl_id column), so
    adapter.read_model_level_availability is used, not read_measurements."""
    availability = {}
    for layer in ("metabolomics", "mirna", "genome_signatures"):
        try:
            result = adapter_module.read_model_level_availability(
                [model_id], layer, data_dir=data_dir
            )
        except FileNotFoundError as exc:
            availability[layer] = {"status": "unavailable", "reason": str(exc), "available": None}
            continue
        except Exception as exc:  # noqa: BLE001 -- reported, never raised past this point
            availability[layer] = {"status": "error", "reason": str(exc), "available": None}
            continue
        availability[layer] = {"status": "ok", "reason": None, "available": result[model_id]}
    return availability


def build_measured_context(snapshot, model_id, adapter, data_dir=None):
    """Assembles one model's deterministic, measured-only biological context packet from an
    S03 snapshot. Copies native values/narratives/warnings verbatim; reads only extra event
    rows and supplementary-assay availability as additive context (never overwriting anything
    already in the snapshot). Adds no template claim beyond describing scope and disclosing
    caveats already implied by the native evidence -- see backend/research_service/templates.py's
    own module docstring for the exact boundary.
    """
    native_result = snapshot["native_result"]
    line, bucket = _find_model_result(native_result, model_id)
    warnings = line.get("warnings", [])
    missing_evidence = line.get("missing_evidence", [])

    genes = []
    for gene in line.get("per_gene", []):
        ensembl_id = gene["ensembl_id"]
        layers_out = {}
        for layer_name, layer_info in gene["layers"].items():
            layer_copy = copy.deepcopy(layer_info)  # native value/unit/d/state, untouched
            layer_copy["unit_description"] = templates.LAYER_UNIT_DESCRIPTIONS.get(layer_name)
            caveats = _layer_caveats(layer_name, layer_info)
            if layer_name == "dependency" and layer_info.get("pan_essential"):
                pan_essential_message = _warning_message_for_gene(warnings, ensembl_id)
                if pan_essential_message:
                    caveats.append({"type": "pan_essential", "message": pan_essential_message})
            if layer_name in EVENT_LAYER_NAMES and layer_info.get("state") != "not_assayed":
                layer_copy["events"] = _read_event_rows(
                    adapter, model_id, ensembl_id, layer_name, data_dir
                )
                if layer_copy["events"]["total_count"] > MAX_DISPLAYED_EVENTS:
                    caveats.append({
                        "type": "events_truncated",
                        "message": layer_copy["events"]["truncation_note"],
                    })
            layer_copy["caveats"] = caveats
            layers_out[layer_name] = layer_copy

        genes.append({
            "ensembl_id": ensembl_id,
            "symbol": gene["symbol"],
            "role": gene["role"],  # inclusion/exclusion -- the researcher's intent, unchanged
            "gene_class": next(
                (info.get("gene_class") for info in gene["layers"].values() if "gene_class" in info),
                "unknown",
            ),  # "unknown" stays "unknown" here -- never guessed at this or any later stage
            "d_gene": gene["d_gene"],
            "layers": layers_out,
            "narrative": gene.get("narrative"),  # native prose, preserved verbatim
            "evidence_caveats": _missing_evidence_caveats(missing_evidence, ensembl_id),
        })

    return {
        "schema_version": "measured-context-v1",
        "query_id": snapshot["query_id"],
        "model_id": model_id,
        "partition_bucket": bucket,
        "scoring_snapshot": {
            "D": line.get("D"),
            "confidence_tier": line.get("confidence_tier"),
            "veto": line.get("veto"),
            "rho_bar": line.get("rho_bar"),
            "m_eff": line.get("m_eff"),
            "m": line.get("m"),
            "hpa_agreement": line.get("hpa_agreement"),
            "geo_agreement": line.get("geo_agreement"),
            "warnings": copy.deepcopy(warnings),  # verbatim, never re-derived
        },
        "genes": genes,
        "model_context": {
            "is_problematic": line.get("is_problematic"),
            "tumour_representativeness": line.get("tumour_representativeness"),
            "tumour_representativeness_note": line.get("tumour_representativeness_note"),
            "hallmark_tags": line.get("hallmark_tags"),
            "hallmark_tags_note": line.get("hallmark_tags_note"),
            "supplementary_assay_availability": _read_supplementary_availability(
                adapter, model_id, data_dir
            ),
        },
        "boundary_statement": native_result.get("boundary_statement"),
    }
