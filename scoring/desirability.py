import json

RNA_CONSTANTS_PATH = "scoring/resources/desirability_constants.json"


def load_rna_constants(path=RNA_CONSTANTS_PATH):
    with open(path) as f:
        return json.load(f)
