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
