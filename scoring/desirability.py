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
