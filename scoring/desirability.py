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
