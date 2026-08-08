import json

RNA_CONSTANTS_PATH = "scoring/resources/desirability_constants.json"


def load_rna_constants(path=RNA_CONSTANTS_PATH):
    with open(path) as f:
        return json.load(f)


EXTENDED_CONSTANTS_PATH = "scoring/resources/desirability_constants_extended.json"


def load_extended_constants(path=EXTENDED_CONSTANTS_PATH):
    with open(path) as f:
        return json.load(f)


ESSENTIALITY_CONSTANTS_PATH = "scoring/resources/common_essential_genes.json"


def load_essentiality_constants(path=ESSENTIALITY_CONSTANTS_PATH):
    with open(path) as f:
        return json.load(f)


def _flip_role(direction):
    return "exclusion" if direction == "inclusion" else "inclusion"


def get_lt(ensembl_id, layer, rna_constants=None, extended_constants=None, lineage=None):
    if layer == "rna":
        if lineage is not None and lineage in rna_constants["per_lineage"]:
            override = rna_constants["per_lineage"][lineage].get(ensembl_id)
            if override is not None:
                return override["L"], override["T"]
        entry = rna_constants["global"].get(ensembl_id)
        return (entry["L"], entry["T"]) if entry else None

    if layer == "protein":

        entry = extended_constants["protein"]["global"]
        return entry["L"], entry["T"]

    if layer in ("dependency", "copy_number"):
        entry = extended_constants[layer]["global"].get(ensembl_id)
        return (entry["L"], entry["T"]) if entry else None

    raise ValueError(f"unknown layer: {layer}")


R_EXPONENT = 1


def desirability_transform(y, L, T, direction, r=R_EXPONENT):
    if T <= L:
        return None  # no discrimination available for this gene; should not occur post-calibration

    fraction = min(max((y - L) / (T - L), 0.0), 1.0) ** r

    if direction == "inclusion":
        if y < L:
            return 0.0
        if y > T:
            return 1.0
        return fraction

    if direction == "exclusion":
        if y < L:
            return 1.0
        if y > T:
            return 0.0
        return 1.0 - fraction

    raise ValueError(f"unknown direction: {direction}")


def score_rna(y, ensembl_id, direction, rna_constants, lineage=None):
    lt = get_lt(ensembl_id, "rna", rna_constants=rna_constants, lineage=lineage)
    if lt is None:
        return None
    return desirability_transform(y, lt[0], lt[1], direction)


def score_protein(zscore, direction, extended_constants, detected=True):
    if not detected:
        return None
    lt = get_lt(None, "protein", extended_constants=extended_constants)
    return desirability_transform(zscore, lt[0], lt[1], direction)


def score_dependency(chronos_score, ensembl_id, direction, extended_constants):
    lt = get_lt(ensembl_id, "dependency", extended_constants=extended_constants)
    if lt is None:
        return None
    dependency_strength = -chronos_score
    return desirability_transform(dependency_strength, lt[0], lt[1], direction)


def is_pan_essential(ensembl_id, essentiality_constants):
    return ensembl_id in essentiality_constants["pan_essential"]


def score_copy_number(copy_number_value, ensembl_id, direction, extended_constants, gene_class="unknown"):
    if gene_class == "unknown":
        return None
    lt = get_lt(ensembl_id, "copy_number", extended_constants=extended_constants)
    if lt is None:
        return None
    effective_direction = direction if gene_class == "oncogene" else _flip_role(direction)
    return desirability_transform(copy_number_value, lt[0], lt[1], effective_direction)


VUS_D = 0.5

MUTATION_D_TABLE = {

    "driver_hotspot": {"inclusion": 1.0, "exclusion": 0.0},

    "high_lof": {"inclusion": 0.0, "exclusion": 1.0},

    "moderate_vus": {"inclusion": VUS_D, "exclusion": VUS_D},

    "low_or_modifier": {"inclusion": 1.0, "exclusion": 1.0},

    "none_sequenced_clean": {"inclusion": 1.0, "exclusion": 1.0},
}

TUMOUR_SUPPRESSOR_MUTATION_D_TABLE = dict(MUTATION_D_TABLE)
TUMOUR_SUPPRESSOR_MUTATION_D_TABLE["high_lof"] = {
    "inclusion": MUTATION_D_TABLE["high_lof"]["exclusion"],
    "exclusion": MUTATION_D_TABLE["high_lof"]["inclusion"],
}


def classify_mutation_rows(rows, has_mutations):
    if not has_mutations:
        return "no_sequencing"
    if rows is None or len(rows) == 0:
        return "none_sequenced_clean"
    if (rows["is_driver"] | rows["is_hotspot"]).any():
        return "driver_hotspot"
    if (rows["vep_impact"] == "HIGH").any():
        return "high_lof"
    if (rows["vep_impact"] == "MODERATE").any():
        return "moderate_vus"
    return "low_or_modifier"


def score_mutation(rows, has_mutations, direction, gene_class="unknown"):
    if gene_class == "unknown":
        return None
    category = classify_mutation_rows(rows, has_mutations)
    if category == "no_sequencing":
        return None
    table = MUTATION_D_TABLE if gene_class == "oncogene" else TUMOUR_SUPPRESSOR_MUTATION_D_TABLE
    return table[category][direction]


FUSION_D_TABLE = {
    "fusion_present": {"inclusion": 1.0, "exclusion": 0.0},
    "assayed_no_fusion": {"inclusion": 0.0, "exclusion": 1.0},
}

FUSION_CONFIDENCE_FACTORS = {"both": 1.0, "one": 0.5, "neither": 0.25}


def classify_fusion_rows(rows, fusions_state):
    if fusions_state == "not_assayed":
        return "no_assay"
    if rows is None or len(rows) == 0:
        return "assayed_no_fusion"
    return "fusion_present"


def fusion_confidence_factor(rows):
    has_high_confidence = bool(rows["confidence_high"].any())
    has_in_frame = bool(rows["in_frame"].any())
    if has_high_confidence and has_in_frame:
        return FUSION_CONFIDENCE_FACTORS["both"]
    if has_high_confidence or has_in_frame:
        return FUSION_CONFIDENCE_FACTORS["one"]
    return FUSION_CONFIDENCE_FACTORS["neither"]
