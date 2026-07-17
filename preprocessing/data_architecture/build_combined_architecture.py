#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import fmean, median


WORKSPACE = Path(__file__).resolve().parents[2]
BASE = WORKSPACE / "data" / "depmap_24Q4"
RAW = BASE / "raw"
DERIVED = BASE / "derived"
SUPPLEMENTARY = BASE / "supplementary"
ARCHITECTURE = BASE / "architecture"

SPINE = DERIVED / "identity_spine.csv"
CROSSWALK = DERIVED / "profile_crosswalk.csv"
MODEL_MASTER = ARCHITECTURE / "model_master.csv"
GENE_MASTER = ARCHITECTURE / "gene_master.csv"
GENE_ALIASES = ARCHITECTURE / "gene_alias_crosswalk.csv"
CRISPR_CONTEXT = ARCHITECTURE / "model_crispr_context.csv"
FUSION_PAIRS = ARCHITECTURE / "fusion_pair_evidence.csv"
TABLE_REGISTRY = ARCHITECTURE / "table_registry.csv"
JOIN_CONTRACTS = ARCHITECTURE / "join_contracts.csv"
EXTENSION_CONTRACTS = ARCHITECTURE / "extension_table_contracts.csv"
MANIFEST = ARCHITECTURE / "build_manifest.json"

FEATURE_LABEL = re.compile(r"^(.*?)\s+\(([^()]+)\)$")
ENSEMBL_IN_LABEL = re.compile(r"(ENSG\d+)(?:\.\d+)?")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as handle:
        return next(csv.reader(handle))


def first_column_ids(path: Path) -> set[str]:
    values: set[str] = set()
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        for row in reader:
            if row and row[0].strip():
                values.add(row[0].strip())
    return values


def true_value(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes"}


def float_or_none(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_feature_columns(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    by_symbol: dict[str, str] = {}
    by_identifier: dict[str, str] = {}
    for label in header(path)[1:]:
        match = FEATURE_LABEL.match(label)
        if match:
            symbol, identifier = match.groups()
            by_symbol.setdefault(symbol.upper(), label)
            by_identifier.setdefault(identifier.split(".")[0], label)
        elif label.startswith("ENSG"):
            by_identifier.setdefault(label.split(".")[0], label)
    return by_symbol, by_identifier


def fusion_gene(value: str) -> tuple[str, str]:
    match = FEATURE_LABEL.match(value.strip())
    symbol = match.group(1) if match else value.strip().split()[0]
    ensembl_match = ENSEMBL_IN_LABEL.search(value)
    ensembl = ensembl_match.group(1) if ensembl_match else ""
    return symbol, ensembl


def assay_model_sets(crosswalk: list[dict[str, str]]) -> dict[str, set[str]]:
    profile_to_model = {row["profile_id"]: row["model_id"] for row in crosswalk}
    unfiltered_profiles = first_column_ids(RAW / "OmicsFusionUnfilteredProfile.csv")
    return {
        "fusion_assay": {profile_to_model[p] for p in unfiltered_profiles if p in profile_to_model},
        "expression": first_column_ids(DERIVED / "expression_default_rna_model.csv"),
        "mutation_event": first_column_ids(DERIVED / "mutations_default_dna_gene_summary.csv"),
        "fusion_call": first_column_ids(DERIVED / "fusions_filtered_model.csv"),
        "signature": first_column_ids(DERIVED / "signatures_model.csv"),
        "crispr_effect": first_column_ids(RAW / "CRISPRGeneEffect.csv"),
        "crispr_dependency": first_column_ids(RAW / "CRISPRGeneDependency.csv"),
        "copy_number": first_column_ids(RAW / "OmicsAbsoluteCNGene.csv"),
        "condition": first_column_ids(DERIVED / "model_condition_context.csv"),
        "crispr_qc": first_column_ids(DERIVED / "crispr_screen_qc_model.csv"),
        "crispr_efficacy": first_column_ids(RAW / "CRISPRInferredModelEfficacy.csv"),
        "crispr_growth": first_column_ids(RAW / "CRISPRInferredModelGrowthRate.csv"),
        "metabolomics": first_column_ids(SUPPLEMENTARY / "metabolomics_model.csv"),
        "proteomics": first_column_ids(SUPPLEMENTARY / "proteomics_model.csv"),
        "mirna": first_column_ids(SUPPLEMENTARY / "mirna_model.csv"),
    }


def count_by_model(path: Path, model_column: str = "model_id") -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as handle:
        for row in csv.DictReader(handle):
            counts[row[model_column]] += 1
    return counts


def build_model_master(spine: list[dict[str, str]], sets: dict[str, set[str]]) -> list[dict[str, object]]:
    condition_counts = count_by_model(DERIVED / "model_condition_context.csv")
    qc_counts = count_by_model(DERIVED / "crispr_screen_qc_model.csv")
    rows: list[dict[str, object]] = []
    for source in spine:
        model_id = source["model_id"]
        row: dict[str, object] = dict(source)
        row.update(
            {
                "has_expression": int(model_id in sets["expression"]),
                "has_mutation_assay": int(source["has_default_dna"] == "1"),
                "has_any_mutation_event": int(model_id in sets["mutation_event"]),
                "has_fusion_assay": int(model_id in sets["fusion_assay"]),
                "has_any_filtered_fusion_call": int(model_id in sets["fusion_call"]),
                "has_genomic_signatures": int(model_id in sets["signature"]),
                "has_crispr_gene_effect": int(model_id in sets["crispr_effect"]),
                "has_crispr_dependency": int(model_id in sets["crispr_dependency"]),
                "has_absolute_copy_number": int(model_id in sets["copy_number"]),
                "has_model_condition": int(model_id in sets["condition"]),
                "n_model_conditions": condition_counts.get(model_id, 0),
                "has_crispr_screen_qc": int(model_id in sets["crispr_qc"]),
                "n_crispr_screens": qc_counts.get(model_id, 0),
                "has_crispr_efficacy": int(model_id in sets["crispr_efficacy"]),
                "has_crispr_growth_rate": int(model_id in sets["crispr_growth"]),
                "has_supplementary_metabolomics": int(model_id in sets["metabolomics"]),
                "has_supplementary_proteomics": int(model_id in sets["proteomics"]),
                "has_supplementary_mirna": int(model_id in sets["mirna"]),
                "eligible_direct_fusion_route": int(model_id in sets["fusion_assay"]),
                "eligible_codependency_route": int(
                    model_id in sets["crispr_effect"] and model_id in sets["crispr_dependency"]
                ),
                "eligible_conditional_dependency_route": int(
                    model_id in sets["crispr_effect"]
                    and (
                        model_id in sets["expression"]
                        or model_id in sets["copy_number"]
                        or source["has_default_dna"] == "1"
                    )
                ),
                "release": "DepMap 24Q4 Public",
            }
        )
        rows.append(row)
    if len({row["model_id"] for row in rows}) != len(rows):
        raise ValueError("model_master would not be one row per model_id")
    return rows


def build_crispr_context() -> list[dict[str, object]]:
    qc_by_model: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_rows(DERIVED / "crispr_screen_qc_model.csv"):
        qc_by_model[row["model_id"]].append(row)

    library_columns = ["Achilles-Avana-2D", "Achilles-Humagne-CD-2D", "Achilles-KY-2D"]

    def inferred(path: Path) -> dict[str, dict[str, str]]:
        values: dict[str, dict[str, str]] = {}
        for row in read_rows(path):
            values[row["ModelID"]] = {library: row.get(library, "") for library in library_columns}
        return values

    efficacy = inferred(RAW / "CRISPRInferredModelEfficacy.csv")
    growth = inferred(RAW / "CRISPRInferredModelGrowthRate.csv")
    model_ids = sorted(first_column_ids(RAW / "CRISPRGeneEffect.csv"))
    rows: list[dict[str, object]] = []
    for model_id in model_ids:
        screens = qc_by_model.get(model_id, [])

        def values(field: str) -> list[float]:
            return [value for row in screens if (value := float_or_none(row.get(field, ""))) is not None]

        passes = [true_value(row.get("PassesQC", "")) for row in screens]
        eff_by_library = efficacy.get(model_id, {})
        growth_by_library = growth.get(model_id, {})
        eff_observed = [value for value in (float_or_none(eff_by_library.get(lib, "")) for lib in library_columns) if value is not None]
        growth_observed = [value for value in (float_or_none(growth_by_library.get(lib, "")) for lib in library_columns) if value is not None]
        observed_libraries = [
            library for library in library_columns
            if eff_by_library.get(library, "") or growth_by_library.get(library, "")
        ]
        rows.append(
            {
                "model_id": model_id,
                "n_screens": len(screens),
                "n_screens_passing_qc": sum(passes),
                "any_screen_passes_qc": int(any(passes)),
                "all_screens_pass_qc": int(bool(passes) and all(passes)),
                "mean_screen_roc_auc": fmean(values("ScreenROCAUC")) if values("ScreenROCAUC") else "",
                "mean_screen_fpr": fmean(values("ScreenFPR")) if values("ScreenFPR") else "",
                "mean_screen_nnmd": fmean(values("ScreenNNMD")) if values("ScreenNNMD") else "",
                "mean_cas_activity": fmean(values("CasActivity")) if values("CasActivity") else "",
                "median_screen_doubling_time": median(values("ScreenDoublingTime")) if values("ScreenDoublingTime") else "",
                "mean_inferred_model_efficacy": fmean(eff_observed) if eff_observed else "",
                "mean_inferred_model_growth_rate": fmean(growth_observed) if growth_observed else "",
                "n_inference_libraries": len(observed_libraries),
                "inference_libraries": "|".join(observed_libraries),
                "efficacy_achilles_avana_2d": eff_by_library.get("Achilles-Avana-2D", ""),
                "efficacy_achilles_humagne_cd_2d": eff_by_library.get("Achilles-Humagne-CD-2D", ""),
                "efficacy_achilles_ky_2d": eff_by_library.get("Achilles-KY-2D", ""),
                "growth_rate_achilles_avana_2d": growth_by_library.get("Achilles-Avana-2D", ""),
                "growth_rate_achilles_humagne_cd_2d": growth_by_library.get("Achilles-Humagne-CD-2D", ""),
                "growth_rate_achilles_ky_2d": growth_by_library.get("Achilles-KY-2D", ""),
            }
        )
    return rows


def build_fusion_pairs() -> tuple[list[dict[str, object]], set[str], set[str]]:
    summaries: dict[tuple[str, str, str], dict[str, object]] = {}
    fusion_symbols: set[str] = set()
    fusion_ensembl: set[str] = set()
    for row in read_rows(DERIVED / "fusions_filtered_model.csv"):
        left_symbol, left_ensembl = fusion_gene(row["LeftGene"])
        right_symbol, right_ensembl = fusion_gene(row["RightGene"])
        fusion_symbols.update({left_symbol.upper(), right_symbol.upper()})
        fusion_ensembl.update(value for value in (left_ensembl, right_ensembl) if value)
        gene_a, gene_b = sorted((left_symbol, right_symbol), key=str.upper)
        key = (row["model_id"], gene_a.upper(), gene_b.upper())
        summary = summaries.setdefault(
            key,
            {
                "model_id": row["model_id"],
                "gene_a_symbol": gene_a,
                "gene_b_symbol": gene_b,
                "fusion_pair_id": f"{gene_a}--{gene_b}",
                "n_filtered_calls": 0,
                "total_junction_reads": 0,
                "total_spanning_fragments": 0,
                "max_ffpm": 0.0,
                "has_large_anchor_support": 0,
                "orientations": set(),
                "annotations": set(),
            },
        )
        summary["n_filtered_calls"] = int(summary["n_filtered_calls"]) + 1
        summary["total_junction_reads"] = int(summary["total_junction_reads"]) + int(row["JunctionReadCount"] or 0)
        summary["total_spanning_fragments"] = int(summary["total_spanning_fragments"]) + int(row["SpanningFragCount"] or 0)
        summary["max_ffpm"] = max(float(summary["max_ffpm"]), float(row["FFPM"] or 0))
        summary["has_large_anchor_support"] = max(
            int(summary["has_large_anchor_support"]), int(row["LargeAnchorSupport"].upper() == "YES_LDAS")
        )
        summary["orientations"].add(f"{left_symbol}->{right_symbol}")
        if row["annots"].strip():
            summary["annotations"].add(row["annots"].strip())

    output: list[dict[str, object]] = []
    for key in sorted(summaries):
        row = summaries[key]
        row["orientations"] = "|".join(sorted(row["orientations"]))
        row["annotations"] = "|".join(sorted(row["annotations"]))
        output.append(row)
    return output, fusion_symbols, fusion_ensembl


def build_gene_tables(fusion_symbols: set[str], fusion_ensembl: set[str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    expression_symbol, expression_id = parse_feature_columns(DERIVED / "expression_default_rna_model.csv")
    crispr_symbol, crispr_id = parse_feature_columns(RAW / "CRISPRGeneEffect.csv")
    dependency_symbol, dependency_id = parse_feature_columns(RAW / "CRISPRGeneDependency.csv")
    copy_symbol, copy_id = parse_feature_columns(RAW / "OmicsAbsoluteCNGene.csv")
    mutation_symbols: set[str] = set()
    mutation_ensembl: set[str] = set()
    with (DERIVED / "mutations_default_dna_gene_summary.csv").open(
        "r", encoding="utf-8-sig", newline="", errors="replace"
    ) as handle:
        for row in csv.DictReader(handle):
            mutation_symbols.add(row["hugo_symbol"].upper())
            mutation_ensembl.add(row["gene_id"].split(".")[0])
    common_essential = {
        FEATURE_LABEL.match(value).group(1).upper() if FEATURE_LABEL.match(value) else value.upper()
        for value in first_column_ids(RAW / "CRISPRInferredCommonEssentials.csv")
    }

    genes: list[dict[str, object]] = []
    aliases: list[dict[str, object]] = []
    alias_seen: set[tuple[str, str, str]] = set()
    for source in read_rows(RAW / "Gene.csv"):
        symbol = source["symbol"]
        symbol_key = symbol.upper()
        ensembl = source["ensembl_gene_id"].split(".")[0]
        entrez = source["entrez_id"]

        def feature(by_symbol: dict[str, str], by_id: dict[str, str], identifier: str) -> str:
            return by_id.get(identifier, "") or by_symbol.get(symbol_key, "")

        row = {
            "gene_id": source["hgnc_id"] or ensembl or entrez or symbol,
            "hgnc_id": source["hgnc_id"],
            "symbol": symbol,
            "gene_name": source["name"],
            "status": source["status"],
            "locus_group": source["locus_group"],
            "locus_type": source["locus_type"],
            "chromosomal_location": source["location"],
            "entrez_id": entrez,
            "ensembl_gene_id": ensembl,
            "uniprot_ids": source["uniprot_ids"],
            "alias_symbols": source["alias_symbol"],
            "previous_symbols": source["prev_symbol"],
            "is_protein_coding": int(source["locus_group"] == "protein-coding gene"),
            "expression_column": feature(expression_symbol, expression_id, ensembl),
            "crispr_effect_column": feature(crispr_symbol, crispr_id, entrez),
            "crispr_dependency_column": feature(dependency_symbol, dependency_id, entrez),
            "copy_number_column": feature(copy_symbol, copy_id, entrez),
            "has_observed_mutation_event": int(symbol_key in mutation_symbols or ensembl in mutation_ensembl),
            "has_observed_filtered_fusion_call": int(symbol_key in fusion_symbols or ensembl in fusion_ensembl),
            "is_crispr_common_essential": int(symbol_key in common_essential),
        }
        row.update(
            {
                "has_expression_feature": int(bool(row["expression_column"])),
                "has_crispr_effect_feature": int(bool(row["crispr_effect_column"])),
                "has_crispr_dependency_feature": int(bool(row["crispr_dependency_column"])),
                "has_copy_number_feature": int(bool(row["copy_number_column"])),
            }
        )
        row["eligible_core_gene"] = int(
            bool(row["crispr_effect_column"])
            and (bool(row["expression_column"]) or bool(row["copy_number_column"]))
        )
        genes.append(row)

        identifiers = [
            (symbol, "canonical_symbol"),
            (source["hgnc_id"], "hgnc_id"),
            (ensembl, "ensembl_gene_id"),
            (entrez, "entrez_id"),
        ]
        identifiers.extend((value, "alias_symbol") for value in source["alias_symbol"].split("|") if value)
        identifiers.extend((value, "previous_symbol") for value in source["prev_symbol"].split("|") if value)
        for identifier, identifier_type in identifiers:
            if not identifier:
                continue
            key = (identifier.upper(), symbol_key, identifier_type)
            if key in alias_seen:
                continue
            alias_seen.add(key)
            aliases.append(
                {
                    "identifier": identifier,
                    "identifier_normalized": identifier.upper(),
                    "identifier_type": identifier_type,
                    "canonical_symbol": symbol,
                    "gene_id": row["gene_id"],
                }
            )
    return genes, aliases


def registry_rows() -> list[dict[str, object]]:
    return [
        {"table": "model_master", "status": "materialized", "grain": "one row per model_id", "primary_key": "model_id", "path": str(MODEL_MASTER), "role": "Identity, cancer context, coverage, and route eligibility."},
        {"table": "gene_master", "status": "materialized", "grain": "one row per HGNC gene record", "primary_key": "gene_id", "path": str(GENE_MASTER), "role": "Canonical gene identifiers and exact matrix-column routing."},
        {"table": "gene_alias_crosswalk", "status": "materialized", "grain": "one identifier-to-canonical-gene assertion", "primary_key": "identifier_normalized + canonical_symbol + identifier_type", "path": str(GENE_ALIASES), "role": "Resolves symbols, aliases, HGNC, Ensembl, and Entrez identifiers."},
        {"table": "profile_crosswalk", "status": "existing", "grain": "one row per profile_id", "primary_key": "profile_id", "path": str(CROSSWALK), "role": "Audited ProfileID to ModelID bridge and default-profile record."},
        {"table": "expression", "status": "existing", "grain": "model_id by gene matrix", "primary_key": "model_id", "path": str(DERIVED / "expression_default_rna_model.csv"), "role": "Expression and conditional-state evidence."},
        {"table": "mutations", "status": "existing", "grain": "one row per model_id and gene", "primary_key": "model_id + gene_id", "path": str(DERIVED / "mutations_default_dna_gene_summary.csv"), "role": "Gene alteration and loss-of-function evidence."},
        {"table": "fusion_calls", "status": "existing", "grain": "one filtered fusion call per model", "primary_key": "event row", "path": str(DERIVED / "fusions_filtered_model.csv"), "role": "Original direct-fusion evidence with breakpoints and read support."},
        {"table": "fusion_pair_evidence", "status": "materialized", "grain": "one row per model_id and unordered gene pair", "primary_key": "model_id + fusion_pair_id", "path": str(FUSION_PAIRS), "role": "Query-ready aggregation of direct fusion calls."},
        {"table": "crispr_gene_effect", "status": "existing", "grain": "model_id by gene matrix", "primary_key": "model_id", "path": str(RAW / "CRISPRGeneEffect.csv"), "role": "Quantitative knockout effect for co-dependency and conditional dependency."},
        {"table": "crispr_dependency", "status": "existing", "grain": "model_id by gene matrix", "primary_key": "model_id", "path": str(RAW / "CRISPRGeneDependency.csv"), "role": "Dependency probability as corroborating knockout evidence."},
        {"table": "copy_number", "status": "existing", "grain": "model_id by gene matrix", "primary_key": "model_id", "path": str(RAW / "OmicsAbsoluteCNGene.csv"), "role": "Deletion and amplification context."},
        {"table": "genomic_signatures", "status": "existing", "grain": "one row per model_id", "primary_key": "model_id", "path": str(DERIVED / "signatures_model.csv"), "role": "MSI, ploidy, LOH, WGD, CIN, and aneuploidy context."},
        {"table": "crispr_model_context", "status": "materialized", "grain": "one row per CRISPR-covered model_id", "primary_key": "model_id", "path": str(CRISPR_CONTEXT), "role": "Aggregated CRISPR QC, efficacy, growth rate, and library context."},
        {"table": "model_condition", "status": "existing", "grain": "one row per model condition", "primary_key": "ModelConditionID", "path": str(DERIVED / "model_condition_context.csv"), "role": "Media, culture format, and treatment context."},
        {"table": "supplementary_omics", "status": "support_only", "grain": "separate model-level matrices", "primary_key": "model_id", "path": str(SUPPLEMENTARY), "role": "Older miRNA, metabolomics, and proteomics; excluded from the 24Q4 core score."},
    ]


def join_rows() -> list[dict[str, object]]:
    return [
        {"from_table": "profile_crosswalk", "to_table": "model_master", "join_key": "model_id", "cardinality": "many-to-one", "join_type": "validated left", "used_by": "profile attribution", "missing_rule": "Unknown profiles fail the upstream identity audit."},
        {"from_table": "expression", "to_table": "model_master", "join_key": "model_id", "cardinality": "one-to-zero-or-one after gene projection", "join_type": "left", "used_by": "conditional dependency", "missing_rule": "NA means expression was not observed; never replace with zero."},
        {"from_table": "mutations", "to_table": "model_master", "join_key": "model_id", "cardinality": "many-to-one before gene filter", "join_type": "left after filtering Gene A/B", "used_by": "conditional dependency", "missing_rule": "No event becomes zero only when has_mutation_assay=1; otherwise NA."},
        {"from_table": "fusion_pair_evidence", "to_table": "model_master", "join_key": "model_id", "cardinality": "zero-or-one after pair filter", "join_type": "left", "used_by": "direct fusion", "missing_rule": "No pair call becomes zero only when has_fusion_assay=1; otherwise NA."},
        {"from_table": "crispr_gene_effect", "to_table": "model_master", "join_key": "model_id", "cardinality": "one-to-zero-or-one after gene projection", "join_type": "left", "used_by": "co-dependency and conditional dependency", "missing_rule": "NA means no usable knockout measurement."},
        {"from_table": "crispr_dependency", "to_table": "model_master", "join_key": "model_id", "cardinality": "one-to-zero-or-one after gene projection", "join_type": "left", "used_by": "dependency corroboration", "missing_rule": "NA means no usable dependency probability."},
        {"from_table": "copy_number", "to_table": "model_master", "join_key": "model_id", "cardinality": "one-to-zero-or-one after gene projection", "join_type": "left", "used_by": "conditional dependency", "missing_rule": "NA means copy number was not observed; it is not diploid zero."},
        {"from_table": "genomic_signatures", "to_table": "model_master", "join_key": "model_id", "cardinality": "one-to-zero-or-one", "join_type": "left", "used_by": "biological context and confounder audit", "missing_rule": "Retain NA."},
        {"from_table": "crispr_model_context", "to_table": "model_master", "join_key": "model_id", "cardinality": "one-to-zero-or-one", "join_type": "left", "used_by": "confidence and QC", "missing_rule": "Missing QC lowers evidence confidence; it is not a failed screen."},
        {"from_table": "gene_alias_crosswalk", "to_table": "gene_master", "join_key": "gene_id", "cardinality": "many-to-one", "join_type": "validated lookup", "used_by": "input resolution", "missing_rule": "Unknown or ambiguous identifiers stop the query."},
    ]


def extension_rows() -> list[dict[str, object]]:
    return [
        {"module": "population_genomic_validation", "required_table": "patient_master", "grain": "one row per patient/sample", "join_keys": "patient_id", "status": "not_materialized", "purpose": "Cancer type, clinical outcome, and covariates for an ISLE-like validation layer."},
        {"module": "population_genomic_validation", "required_table": "patient_gene_state", "grain": "patient_id by gene", "join_keys": "patient_id + gene_id", "status": "not_materialized", "purpose": "TCGA mutation, copy-number, and expression inactivation states."},
        {"module": "paralog_buffering", "required_table": "paralog_pair_reference", "grain": "one row per paralog pair", "join_keys": "gene_a_id + gene_b_id", "status": "not_materialized", "purpose": "Sequence similarity, conservation, and shared interaction partners."},
        {"module": "knowledge_graph", "required_table": "knowledge_graph_nodes", "grain": "one row per biological entity", "join_keys": "node_id", "status": "not_materialized", "purpose": "Genes, pathways, diseases, drugs, and complexes."},
        {"module": "knowledge_graph", "required_table": "knowledge_graph_edges", "grain": "one typed source-target relation", "join_keys": "source_node_id + target_node_id + relation_type", "status": "not_materialized", "purpose": "Explainable paths and shared biological mechanisms."},
        {"module": "metabolic_simulation", "required_table": "metabolic_reaction_master", "grain": "one row per reaction", "join_keys": "reaction_id", "status": "not_materialized", "purpose": "Reaction bounds, metabolites, compartments, and objective participation."},
        {"module": "metabolic_simulation", "required_table": "gene_reaction_mapping", "grain": "one gene-protein-reaction rule", "join_keys": "gene_id + reaction_id", "status": "not_materialized", "purpose": "Maps a queried gene knockout to metabolic reactions."},
        {"module": "metabolic_simulation", "required_table": "pair_simulation_result", "grain": "model_id by gene pair", "join_keys": "model_id + gene_a_id + gene_b_id", "status": "not_materialized", "purpose": "Single- and double-knockout growth predictions."},
    ]


def main() -> None:
    required = [SPINE, CROSSWALK, RAW / "Gene.csv", RAW / "CRISPRGeneEffect.csv"]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required architecture inputs: {', '.join(missing)}")
    ARCHITECTURE.mkdir(parents=True, exist_ok=True)

    spine = read_rows(SPINE)
    crosswalk = read_rows(CROSSWALK)
    sets = assay_model_sets(crosswalk)
    model_master = build_model_master(spine, sets)
    crispr_context = build_crispr_context()
    fusion_pairs, fusion_symbols, fusion_ensembl = build_fusion_pairs()
    gene_master, aliases = build_gene_tables(fusion_symbols, fusion_ensembl)

    write_rows(MODEL_MASTER, model_master)
    write_rows(GENE_MASTER, gene_master)
    write_rows(GENE_ALIASES, aliases)
    write_rows(CRISPR_CONTEXT, crispr_context)
    write_rows(FUSION_PAIRS, fusion_pairs)
    write_rows(TABLE_REGISTRY, registry_rows())
    write_rows(JOIN_CONTRACTS, join_rows())
    write_rows(EXTENSION_CONTRACTS, extension_rows())

    model_ids = {row["model_id"] for row in model_master}
    unknown_sets = {name: sorted(values - model_ids) for name, values in sets.items() if values - model_ids}
    if unknown_sets:
        raise ValueError(f"Architecture layers contain unknown model IDs: {unknown_sets}")

    counts = {
        "model_master": len(model_master),
        "gene_master": len(gene_master),
        "gene_alias_crosswalk": len(aliases),
        "model_crispr_context": len(crispr_context),
        "fusion_pair_evidence": len(fusion_pairs),
        "table_registry": len(registry_rows()),
        "join_contracts": len(join_rows()),
        "extension_table_contracts": len(extension_rows()),
    }
    coverage = {
        field: sum(int(row[field]) for row in model_master)
        for field in model_master[0]
        if field.startswith("has_") or field.startswith("eligible_")
    }
    manifest = {
        "release": "DepMap 24Q4 Public",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "master_grain": "one row per model_id",
        "counts": counts,
        "model_coverage": coverage,
        "validation": {
            "unique_model_ids": len(model_ids) == len(model_master),
            "unknown_model_ids_across_layers": 0,
            "supplementary_omics_in_core_score": False,
            "missing_measurements_encoded_as_zero": False,
        },
        "outputs": {row["table"]: row["path"] for row in registry_rows() if row["status"] == "materialized"},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    for name, count in counts.items():
        print(f"{name}: {count:,} rows")
    print(f"manifest: {MANIFEST}")


if __name__ == "__main__":
    main()
