def impute_protein_detection_floor(protein, floor_percentile=1.0):
    detected_zscores = protein.loc[protein["detected"], "zscore"]
    floor_value = detected_zscores.quantile(floor_percentile / 100)

    imputed = protein.copy()
    n_imputed = int((~imputed["detected"]).sum())
    imputed.loc[~imputed["detected"], "zscore"] = floor_value

    report = {
        "floor_percentile": floor_percentile,
        "floor_value": round(float(floor_value), 4),
        "n_rows_imputed": n_imputed,
    }
    return imputed, report
