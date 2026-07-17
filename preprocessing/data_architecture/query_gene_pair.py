#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


WORKSPACE = Path(__file__).resolve().parents[2]
BASE = WORKSPACE / "data" / "depmap_24Q4"
RAW = BASE / "raw"
DERIVED = BASE / "derived"
ARCHITECTURE = BASE / "architecture"


def resolve_gene(query: str, genes: pd.DataFrame, aliases: pd.DataFrame) -> pd.Series:
    normalized = query.strip().upper()
    canonical = genes[genes["symbol"].str.upper() == normalized]
    if len(canonical) == 1:
        return canonical.iloc[0]
    matches = aliases[aliases["identifier_normalized"] == normalized]
    symbols = sorted(matches["canonical_symbol"].drop_duplicates())
    if not symbols:
        raise ValueError(f"Unknown gene identifier: {query}")
    if len(symbols) > 1:
        raise ValueError(f"Ambiguous gene identifier {query}: {', '.join(symbols)}")
    return genes[genes["symbol"] == symbols[0]].iloc[0]


def project_matrix(path: Path, feature_rows: list[pd.Series], feature_field: str, value_name: str) -> pd.DataFrame:
    first_column = pd.read_csv(path, nrows=0).columns[0]
    selected: list[str] = []
    rename = {first_column: "model_id"}
    for label, gene in zip(("gene_a", "gene_b"), feature_rows):
        feature = str(gene[feature_field]) if pd.notna(gene[feature_field]) else ""
        if feature:
            selected.append(feature)
            rename[feature] = f"{label}_{value_name}"
    if not selected:
        return pd.DataFrame(columns=["model_id"])
    frame = pd.read_csv(path, usecols=[first_column, *dict.fromkeys(selected)])
    frame = frame.rename(columns=rename)
    if feature_rows[0][feature_field] == feature_rows[1][feature_field] and selected:
        frame[f"gene_b_{value_name}"] = frame[f"gene_a_{value_name}"]
    return frame


def mutation_projection(symbol: str, prefix: str, mutation_assayed_models: set[str]) -> pd.DataFrame:
    columns = [
        "model_id", "hugo_symbol", "n_variants", "max_allele_fraction", "has_likely_lof",
        "has_hotspot", "has_oncogene_high_impact", "has_tumor_suppressor_high_impact", "highest_vep_impact",
    ]
    parts = []
    for chunk in pd.read_csv(DERIVED / "mutations_default_dna_gene_summary.csv", usecols=columns, chunksize=100_000):
        parts.append(chunk[chunk["hugo_symbol"].str.upper() == symbol.upper()])
    observed = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=columns)
    observed = observed.drop(columns="hugo_symbol").rename(
        columns={column: f"{prefix}_mutation_{column}" for column in observed.columns if column != "model_id"}
    )
    measured = pd.DataFrame({"model_id": sorted(mutation_assayed_models)})
    result = measured.merge(observed, on="model_id", how="left", validate="one_to_one")
    zero_fields = [
        f"{prefix}_mutation_n_variants", f"{prefix}_mutation_max_allele_fraction",
        f"{prefix}_mutation_has_likely_lof", f"{prefix}_mutation_has_hotspot",
        f"{prefix}_mutation_has_oncogene_high_impact", f"{prefix}_mutation_has_tumor_suppressor_high_impact",
    ]
    result[zero_fields] = result[zero_fields].fillna(0)
    return result


def left_join(base: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    if len(evidence.columns) == 1:
        return base
    before = len(base)
    result = base.merge(evidence, on="model_id", how="left", validate="one_to_one")
    if len(result) != before:
        raise ValueError("A side-table join changed the model-master row count")
    return result


def safe_corr(frame: pd.DataFrame, a: str, b: str, method: str) -> float | None:
    complete = frame[[a, b]].dropna()
    if len(complete) < 3:
        return None
    value = complete[a].corr(complete[b], method=method)
    return None if pd.isna(value) else float(value)


def conditional_summary(frame: pd.DataFrame, loss_column: str, dependency_column: str) -> dict[str, object]:
    eligible = frame[[loss_column, dependency_column]].dropna()
    lost = eligible[eligible[loss_column] == 1][dependency_column]
    intact = eligible[eligible[loss_column] == 0][dependency_column]
    return {
        "eligible_models": len(eligible),
        "loss_models": len(lost),
        "mean_effect_when_partner_lost": None if lost.empty else float(lost.mean()),
        "mean_effect_when_partner_not_lost": None if intact.empty else float(intact.mean()),
        "effect_difference": None if lost.empty or intact.empty else float(lost.mean() - intact.mean()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-a", required=True)
    parser.add_argument("--gene-b", required=True)
    parser.add_argument("--output-dir", type=Path, default=ARCHITECTURE / "pair_queries")
    args = parser.parse_args()

    genes = pd.read_csv(ARCHITECTURE / "gene_master.csv", dtype=str, keep_default_na=False)
    aliases = pd.read_csv(ARCHITECTURE / "gene_alias_crosswalk.csv", dtype=str, keep_default_na=False)
    gene_a = resolve_gene(args.gene_a, genes, aliases)
    gene_b = resolve_gene(args.gene_b, genes, aliases)
    selected = [gene_a, gene_b]

    result = pd.read_csv(ARCHITECTURE / "model_master.csv")
    if result["model_id"].duplicated().any():
        raise ValueError("model_master is not unique on model_id")

    projections = [
        project_matrix(DERIVED / "expression_default_rna_model.csv", selected, "expression_column", "expression_log2_tpm_plus_1"),
        project_matrix(RAW / "CRISPRGeneEffect.csv", selected, "crispr_effect_column", "crispr_gene_effect"),
        project_matrix(RAW / "CRISPRGeneDependency.csv", selected, "crispr_dependency_column", "crispr_dependency_probability"),
        project_matrix(RAW / "OmicsAbsoluteCNGene.csv", selected, "copy_number_column", "absolute_copy_number"),
    ]
    for projection in projections:
        result = left_join(result, projection)

    mutation_assayed = set(result.loc[result["has_mutation_assay"] == 1, "model_id"])
    result = left_join(result, mutation_projection(gene_a["symbol"], "gene_a", mutation_assayed))
    result = left_join(result, mutation_projection(gene_b["symbol"], "gene_b", mutation_assayed))
    result = left_join(result, pd.read_csv(DERIVED / "signatures_model.csv"))
    result = left_join(result, pd.read_csv(ARCHITECTURE / "model_crispr_context.csv"))

    pair_symbols = sorted((gene_a["symbol"], gene_b["symbol"]), key=str.upper)
    fusions = pd.read_csv(ARCHITECTURE / "fusion_pair_evidence.csv")
    pair_fusions = fusions[
        (fusions["gene_a_symbol"].str.upper() == pair_symbols[0].upper())
        & (fusions["gene_b_symbol"].str.upper() == pair_symbols[1].upper())
    ].drop(columns=["gene_a_symbol", "gene_b_symbol", "fusion_pair_id"])
    result = left_join(result, pair_fusions)
    result["direct_fusion_call"] = pd.NA
    assayed = result["has_fusion_assay"] == 1
    result.loc[assayed, "direct_fusion_call"] = result.loc[assayed, "n_filtered_calls"].notna().astype(int)

    for prefix in ("gene_a", "gene_b"):
        expression = f"{prefix}_expression_log2_tpm_plus_1"
        copy_number = f"{prefix}_absolute_copy_number"
        lof = f"{prefix}_mutation_has_likely_lof"
        if expression in result:
            mean = result[expression].mean()
            std = result[expression].std()
            result[f"{prefix}_expression_z"] = (result[expression] - mean) / std if std else pd.NA
        else:
            result[f"{prefix}_expression_z"] = pd.NA
        loss_parts = []
        observed_parts = []
        if lof in result:
            loss_parts.append(result[lof].eq(1))
            observed_parts.append(result[lof].notna())
        if copy_number in result:
            loss_parts.append(result[copy_number].lt(1))
            observed_parts.append(result[copy_number].notna())
        loss_parts.append(result[f"{prefix}_expression_z"].lt(-1.5))
        observed_parts.append(result[f"{prefix}_expression_z"].notna())
        observed = pd.concat(observed_parts, axis=1).any(axis=1)
        result[f"{prefix}_loss_proxy"] = pd.NA
        result.loc[observed, f"{prefix}_loss_proxy"] = pd.concat(loss_parts, axis=1).any(axis=1).astype(int)

    effect_a = "gene_a_crispr_gene_effect"
    effect_b = "gene_b_crispr_gene_effect"
    result["codependency_route_observed"] = (
        result.get(effect_a, pd.Series(index=result.index, dtype=float)).notna()
        & result.get(effect_b, pd.Series(index=result.index, dtype=float)).notna()
    ).astype(int)
    result["conditional_a_loss_to_b_dependency_observed"] = (
        result["gene_a_loss_proxy"].notna()
        & result.get(effect_b, pd.Series(index=result.index, dtype=float)).notna()
    ).astype(int)
    result["conditional_b_loss_to_a_dependency_observed"] = (
        result["gene_b_loss_proxy"].notna()
        & result.get(effect_a, pd.Series(index=result.index, dtype=float)).notna()
    ).astype(int)

    pair_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{gene_a['symbol']}__{gene_b['symbol']}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table_path = args.output_dir / f"{pair_name}_pair_model_evidence.csv"
    summary_path = args.output_dir / f"{pair_name}_summary.json"
    result.to_csv(table_path, index=False)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "release": "DepMap 24Q4 Public",
        "gene_a": {"input": args.gene_a, "symbol": gene_a["symbol"], "gene_id": gene_a["gene_id"]},
        "gene_b": {"input": args.gene_b, "symbol": gene_b["symbol"], "gene_id": gene_b["gene_id"]},
        "rows": len(result),
        "direct_fusion": {
            "assayed_models": int(result["has_fusion_assay"].sum()),
            "models_with_pair_call": int(result["direct_fusion_call"].eq(1).sum()),
        },
        "codependency": {
            "overlap_models": int(result["codependency_route_observed"].sum()),
            "pearson_gene_effect": safe_corr(result, effect_a, effect_b, "pearson") if effect_a in result and effect_b in result else None,
            "spearman_gene_effect": safe_corr(result, effect_a, effect_b, "spearman") if effect_a in result and effect_b in result else None,
        },
        "conditional_dependency": {
            "a_loss_to_b_dependency": conditional_summary(result, "gene_a_loss_proxy", effect_b) if effect_b in result else None,
            "b_loss_to_a_dependency": conditional_summary(result, "gene_b_loss_proxy", effect_a) if effect_a in result else None,
            "loss_proxy_rule": "likely LoF mutation OR absolute copy number < 1 OR expression z-score < -1.5",
        },
        "interpretation": "These are auditable evidence summaries, not a calibrated probability or validated final score.",
        "output_table": str(table_path),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"resolved: {args.gene_a} -> {gene_a['symbol']}; {args.gene_b} -> {gene_b['symbol']}")
    print(f"pair evidence: {len(result):,} rows -> {table_path}")
    print(f"summary: {summary_path}")


if __name__ == "__main__":
    main()
