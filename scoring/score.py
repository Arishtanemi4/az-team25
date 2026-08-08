import json

import pandas as pd

import combine
import desirability
import evidence
import explain
import io_utils

DATA_DIR = "data/processed"
GENE_ROLE_REFERENCE_PATH = "preprocessing/resources/gene_role_reference.csv"
UNRESOLVED_SYMBOLS_PATH = "preprocessing/resources/unresolved_gene_symbols.json"

COVERAGE_STATE_COLUMNS = [
    "expression_rna_state", "expression_rna_hpa_state", "expression_rna_geo_state",
    "protein_state", "dependency_state", "copy_number_state",
    "mutations_state", "fusions_state",
]

LAYER_UNITS = {
    "rna": "log2(TPM+1)",
    "protein": "z-score",
    "dependency": "Chronos gene effect (dependency strength = -score)",
    "copy_number": "relative copy number (diploid ~= 1.0)",
    "mutation": "categorical -- VEP impact / driver / hotspot rule (desirability.MUTATION_D_TABLE)",
    "fusion": "categorical -- event presence, binary carve-out (desirability.FUSION_D_TABLE)",
}


def _load_unresolved_symbols(path):
    try:
        with open(path) as f:
            data = json.load(f)
    except FileNotFoundError:
        return {}
    return {layer: set(ids) for layer, ids in data.items() if not layer.startswith("_")}


def resolve_genes(tokens, gene_reference_df):
    resolved, ambiguous, unresolved = {}, [], []
    for token in tokens:
        by_id = gene_reference_df[gene_reference_df["ensembl_id"] == token]
        if len(by_id) == 1:
            resolved[token] = token
            continue

        by_symbol = gene_reference_df[gene_reference_df["symbol"] == token]
        if len(by_symbol) == 1:
            resolved[token] = by_symbol.iloc[0]["ensembl_id"]
        elif len(by_symbol) > 1:
            ambiguous.append(token)
        else:
            unresolved.append(token)
    return resolved, ambiguous, unresolved


def resolve_genes_or_raise(tokens, gene_reference_df):
    resolved, ambiguous, unresolved = resolve_genes(tokens, gene_reference_df)
    if ambiguous:
        raise ValueError(f"Ambiguous gene symbol(s) -- specify by Ensembl ID instead: {ambiguous}")
    if unresolved:
        raise ValueError(f"Gene(s) not found in gene_reference: {unresolved}")
    return resolved

STAGE_LABELS = {
    "resolve_genes": "Resolving requested genes against the reference table",
    "filter_candidates": "Filtering candidate cell lines by lineage, disease, and QC criteria",
    "load_tables": "Loading expression, dependency, mutation, copy-number, and fusion evidence tables",
    "correlation": "Computing inter-gene correlation weights across the reference expression panel",
    "score_candidates": "Scoring each candidate cell line against the requested gene panel",
    "rank": "Ranking candidates and splitting results by confidence tier",
    "similarity": "Identifying similar backup cell lines for top-ranked results",
}


COVERAGE_FILTER_COLUMNS = ["ModelID", "msi_high", "ploidy", "metabolomics_available", "mirna_available"]


def filter_candidates(cell_lines_df, coverage_df=None, lineage=None, primary_disease=None,
                       exclude_problematic=False, msi_high=None, ploidy_min=None, ploidy_max=None,
                       require_metabolomics=False, require_mirna=False):
    needs_coverage = (
        msi_high is not None or ploidy_min is not None or ploidy_max is not None
        or require_metabolomics or require_mirna
    )
    if needs_coverage and coverage_df is None:
        raise ValueError(
            "filter_candidates: msi_high/ploidy_min/ploidy_max/require_metabolomics/"
            "require_mirna need coverage_df (data/processed/coverage.csv) -- none was given."
        )

    subset = cell_lines_df
    if lineage is not None:
        subset = subset[subset["lineage"] == lineage]
    if primary_disease is not None:
        subset = subset[subset["primary_disease"] == primary_disease]
    if exclude_problematic:
        subset = subset[~subset["is_problematic"].astype(bool)]

    if needs_coverage:
        subset = subset.merge(coverage_df[COVERAGE_FILTER_COLUMNS], on="ModelID", how="left")
        if msi_high is not None:
            subset = subset[subset["msi_high"] == msi_high]
        if ploidy_min is not None:
            subset = subset[subset["ploidy"] >= ploidy_min]
        if ploidy_max is not None:
            subset = subset[subset["ploidy"] <= ploidy_max]
        if require_metabolomics:
            subset = subset[subset["metabolomics_available"].astype(bool)]
        if require_mirna:
            subset = subset[subset["mirna_available"].astype(bool)]

    return subset, len(subset)


def _load_filtered(path, genes, model_ids, chunksize=1_000_000):
    model_ids = set(model_ids)
    df = io_utils.read_gene_filtered(path, genes, chunksize=chunksize)
    return df[df["ModelID"].isin(model_ids)].reset_index(drop=True)


def load_tables_for_query(query_genes, model_ids, gene_reference_df, data_dir=DATA_DIR,
                           rna_constants_path=desirability.RNA_CONSTANTS_PATH,
                           extended_constants_path=desirability.EXTENDED_CONSTANTS_PATH,
                           essentiality_constants_path=desirability.ESSENTIALITY_CONSTANTS_PATH,
                           gene_role_path=GENE_ROLE_REFERENCE_PATH,
                           unresolved_symbols_path=UNRESOLVED_SYMBOLS_PATH,
                           coverage_df=None, rna_constants=None, extended_constants=None,
                           essentiality_constants=None, gene_role_reference=None,
                           unresolved_symbols=None, preloaded_layer_frames=None):
    if preloaded_layer_frames is not None:
        def _filter_preloaded(name):
            df = preloaded_layer_frames[name]
            if df.empty:
                return df
            mask = df["ensembl_id"].isin(query_genes) & df["ModelID"].isin(model_ids)
            return df[mask]

        expr = _filter_preloaded("expression_rna")
        hpa = _filter_preloaded("expression_rna_hpa")
        geo = _filter_preloaded("expression_rna_geo")
        protein = _filter_preloaded("protein")
        dependency = _filter_preloaded("dependency")
        copy_number = _filter_preloaded("copy_number")
        mutations = _filter_preloaded("mutations")
        fusions = _filter_preloaded("fusions")
    else:
        expr = _load_filtered(f"{data_dir}/expression_rna.csv", query_genes, model_ids)
        hpa = _load_filtered(f"{data_dir}/expression_rna_hpa.csv", query_genes, model_ids)
        geo = _load_filtered(f"{data_dir}/expression_rna_geo.csv", query_genes, model_ids)
        protein = _load_filtered(f"{data_dir}/protein.csv", query_genes, model_ids)
        dependency = _load_filtered(f"{data_dir}/dependency.csv", query_genes, model_ids)
        copy_number = _load_filtered(f"{data_dir}/copy_number.csv", query_genes, model_ids)
        mutations = _load_filtered(f"{data_dir}/mutations.csv", query_genes, model_ids)
        fusions = _load_filtered(f"{data_dir}/fusions.csv", query_genes, model_ids)

    coverage = coverage_df if coverage_df is not None else pd.read_csv(f"{data_dir}/coverage.csv")
    coverage = coverage[coverage["ModelID"].isin(model_ids)].copy()
    coverage["_layer_count"] = (coverage[COVERAGE_STATE_COLUMNS] != "not_assayed").sum(axis=1)

    protein_lookup = {}
    if not protein.empty:
        for row in protein.itertuples(index=False):
            protein_lookup[(row.ModelID, row.ensembl_id)] = (row.zscore, row.detected)

    return {
        "expression_rna": expr.set_index(["ModelID", "ensembl_id"])["log2tpm1"].to_dict() if not expr.empty else {},
        "expression_rna_hpa": hpa.set_index(["ModelID", "ensembl_id"])["log2ntpm1"].to_dict() if not hpa.empty else {},
        "expression_rna_geo": geo.set_index(["ModelID", "ensembl_id"])["log1p_expr"].to_dict() if not geo.empty else {},
        "protein": protein_lookup,
        "dependency": dependency.set_index(["ModelID", "ensembl_id"])["dependency_score"].to_dict() if not dependency.empty else {},
        "copy_number": copy_number.set_index(["ModelID", "ensembl_id"])["copy_number"].to_dict() if not copy_number.empty else {},
        "mutations": {k: v for k, v in mutations.groupby(["ModelID", "ensembl_id"])},
        "fusions": {k: v for k, v in fusions.groupby(["ModelID", "ensembl_id"])},
        "mutations_state": coverage.set_index("ModelID")["mutations_state"].to_dict(),
        "fusions_state": coverage.set_index("ModelID")["fusions_state"].to_dict(),
        "coverage_layer_count": coverage.set_index("ModelID")["_layer_count"].to_dict(),
        "gene_symbols": gene_reference_df.set_index("ensembl_id")["symbol"].to_dict(),
        "rna_constants": rna_constants if rna_constants is not None
            else desirability.load_rna_constants(rna_constants_path),
        "extended_constants": extended_constants if extended_constants is not None
            else desirability.load_extended_constants(extended_constants_path),
        "essentiality_constants": essentiality_constants if essentiality_constants is not None
            else desirability.load_essentiality_constants(essentiality_constants_path),

        "gene_class": (
            gene_role_reference if gene_role_reference is not None
            else pd.read_csv(gene_role_path).set_index("ensembl_id")["gene_class"].to_dict()
        ),

        "unresolved_symbols": (
            unresolved_symbols if unresolved_symbols is not None
            else _load_unresolved_symbols(unresolved_symbols_path)
        ),
    }


def _score_gene_layers(model_id, ensembl_id, role, tables, lineage, has_mutations, fusions_state_line):
    rna_constants = tables["rna_constants"]
    extended_constants = tables["extended_constants"]
    layers = {}

    y_rna = tables["expression_rna"].get((model_id, ensembl_id))
    if y_rna is not None:
        d = desirability.score_rna(y_rna, ensembl_id, role, rna_constants, lineage=lineage)
        layers["rna"] = {"value": y_rna, "unit": LAYER_UNITS["rna"], "d": d, "state": "measured"}
    else:
        layers["rna"] = {"value": None, "unit": LAYER_UNITS["rna"], "d": None, "state": "not_assayed"}

    protein_row = tables["protein"].get((model_id, ensembl_id))
    if protein_row is not None:
        zscore, detected = protein_row
        d = desirability.score_protein(zscore, role, extended_constants, detected=detected)

        state = "measured" if detected else "non_detected"
        value = zscore if detected else None
        layers["protein"] = {"value": value, "unit": LAYER_UNITS["protein"], "d": d, "state": state}
    else:
        layers["protein"] = {"value": None, "unit": LAYER_UNITS["protein"], "d": None, "state": "not_assayed"}

    y_dep = tables["dependency"].get((model_id, ensembl_id))
    if y_dep is not None:
        essentiality_constants = tables.get("essentiality_constants", {"pan_essential": {}})
        pan_essential = desirability.is_pan_essential(ensembl_id, essentiality_constants)
        if pan_essential:

            d = None
        else:
            d = desirability.score_dependency(y_dep, ensembl_id, role, extended_constants)
        state = "excluded_pan_essential" if pan_essential else "measured"
        layers["dependency"] = {
            "value": y_dep, "unit": LAYER_UNITS["dependency"], "d": d, "state": state,
            "pan_essential": pan_essential,
        }
    else:
        layers["dependency"] = {
            "value": None, "unit": LAYER_UNITS["dependency"], "d": None, "state": "not_assayed",
            "pan_essential": False,
        }

    gene_class = tables.get("gene_class", {}).get(ensembl_id, "unknown")

    y_cn = tables["copy_number"].get((model_id, ensembl_id))
    if y_cn is not None:
        d = desirability.score_copy_number(y_cn, ensembl_id, role, extended_constants, gene_class=gene_class)
        layers["copy_number"] = {
            "value": y_cn, "unit": LAYER_UNITS["copy_number"], "d": d, "state": "measured",
            "gene_class": gene_class,
        }
    else:
        layers["copy_number"] = {
            "value": None, "unit": LAYER_UNITS["copy_number"], "d": None, "state": "not_assayed",
            "gene_class": gene_class,
        }

    mutation_rows = tables["mutations"].get((model_id, ensembl_id))
    mutation_category = desirability.classify_mutation_rows(mutation_rows, has_mutations)
    d_mutation = desirability.score_mutation(mutation_rows, has_mutations, role, gene_class=gene_class)
    mutation_state = {
        "no_sequencing": "not_assayed",
        "none_sequenced_clean": "measured_absent",
    }.get(mutation_category, "measured")
    layers["mutation"] = {
        "value": mutation_category, "unit": LAYER_UNITS["mutation"], "d": d_mutation, "state": mutation_state,
        "gene_class": gene_class,
    }

    fusion_rows = tables["fusions"].get((model_id, ensembl_id))
    d_fusion, fusion_confidence_factor = desirability.score_fusion(fusion_rows, fusions_state_line, role)
    fusion_category = desirability.classify_fusion_rows(fusion_rows, fusions_state_line)
    fusion_state = {
        "no_assay": "not_assayed",
        "assayed_no_fusion": "measured_absent",
        "fusion_present": "measured",
    }[fusion_category]
    layers["fusion"] = {
        "value": 0 if fusion_rows is None else len(fusion_rows), "unit": LAYER_UNITS["fusion"],
        "d": d_fusion, "state": fusion_state, "confidence_factor": fusion_confidence_factor,
    }

    return layers, mutation_rows, fusion_rows


def score_one_line(model_id, inclusion_genes, exclusion_genes, tables, correlation_weights,
                    cell_lines_row, layer_weights=None, tier_params=None, build_narrative=True):
    lineage = cell_lines_row.get("lineage")
    has_mutations = tables["mutations_state"].get(model_id, "not_assayed") != "not_assayed"
    fusions_state_line = tables["fusions_state"].get(model_id, "not_assayed")

    query = [(g, "inclusion") for g in inclusion_genes] + [(g, "exclusion") for g in exclusion_genes]

    per_gene_results = []
    rna_values, hpa_values, geo_values = {}, {}, {}
    fusion_context, mutation_context, pan_essential_context = [], [], []

    for ensembl_id, role in query:
        symbol = tables["gene_symbols"].get(ensembl_id, ensembl_id)
        layers, mutation_rows, fusion_rows = _score_gene_layers(
            model_id, ensembl_id, role, tables, lineage, has_mutations, fusions_state_line
        )

        if layers["rna"]["d"] is not None:
            rna_values[ensembl_id] = layers["rna"]["value"]
        y_hpa = tables["expression_rna_hpa"].get((model_id, ensembl_id))
        if y_hpa is not None:
            hpa_values[ensembl_id] = y_hpa  # corroboration only -- never scored (CONSTRAINTS.md #9)
        y_geo = tables["expression_rna_geo"].get((model_id, ensembl_id))
        if y_geo is not None:
            geo_values[ensembl_id] = y_geo  # corroboration/breadth only (METHOD_DECISION.md SS2)

        if mutation_rows is not None and len(mutation_rows) > 0:
            mutation_context.append({
                "ensembl_id": ensembl_id, "symbol": symbol,
                "category": layers["mutation"]["value"], "rows": mutation_rows,
            })
        if fusion_rows is not None and len(fusion_rows) > 0:
            fusion_context.append({
                "ensembl_id": ensembl_id, "symbol": symbol,
                "events": len(fusion_rows), "rows": fusion_rows,
            })
        if layers["dependency"].get("pan_essential"):
            pan_essential_context.append({
                "ensembl_id": ensembl_id, "symbol": symbol,
                "fraction_strong": tables["essentiality_constants"]["pan_essential"].get(ensembl_id),
            })

        layer_d = {layer: info["d"] for layer, info in layers.items() if info["d"] is not None}
        gene_weights = dict(layer_weights) if layer_weights is not None else dict(evidence.LAYER_WEIGHTS)
        if layers["fusion"]["state"] == "measured":
            gene_weights["fusion"] = evidence.LAYER_WEIGHTS["fusion"] * layers["fusion"]["confidence_factor"]
        fusion_weight_zeroed = False
        if role == "exclusion" and layers["fusion"]["state"] == "measured_absent":
            gene_weights["fusion"] = 0.0
            fusion_weight_zeroed = True

        per_gene_results.append({
            "ensembl_id": ensembl_id,
            "symbol": symbol,
            "role": role,
            "d_gene": evidence.combine_gene_evidence(layer_d, weights=gene_weights),
            "layers": layers,
            "missing_layers": [
                l for l, info in layers.items()
                if info["state"] in ("not_assayed", "non_detected")
            ],
            "fusion_weight_zeroed": fusion_weight_zeroed,
        })

    active = [g for g in per_gene_results if g["d_gene"] is not None]
    if not active:
        D, veto_info = None, None
    else:
        d_gene = {g["ensembl_id"]: g["d_gene"] for g in active}
        weights = {eid: correlation_weights["weights"].get(eid, 1.0) for eid in d_gene}
        D, veto_info = combine.combine_across_genes(d_gene, weights)

    hpa_agreement = explain.check_hpa_corroboration(rna_values, hpa_values, tables["rna_constants"])
    geo_agreement = explain.check_geo_corroboration(rna_values, geo_values, tables["rna_constants"])

    return explain.build_result(
        model_id, D, veto_info, per_gene_results,
        correlation_weights["rho_bar"], correlation_weights["m_eff"], len(active),
        cell_lines_row, hpa_agreement=hpa_agreement, geo_agreement=geo_agreement,
        fusion_context=fusion_context, mutation_context=mutation_context,
        pan_essential_context=pan_essential_context,
        tier_params=tier_params, build_narrative=build_narrative,
        unresolved_symbols=tables.get("unresolved_symbols"),
    )
