import numpy as np


def scale_expression_rna(expression_rna):
    return expression_rna


def scale_protein(protein):
    return protein


def scale_metabolomics(metabolomics_presence):
    return metabolomics_presence


def scale_expression_geo(expression_geo):
    df = expression_geo.copy()
    df["log1p_expr"] = np.log1p(df["raw_expr"])
    return df


def scale_mirna(mirna):
    df = mirna.copy()
    df["log1p_value"] = np.log1p(df["raw_value"])
    return df
