import numpy as np


def normalize_hpa(expression_hpa):
    df = expression_hpa.copy()
    df["log2ntpm1"] = np.log2(df["nTPM"] + 1)
    return df
