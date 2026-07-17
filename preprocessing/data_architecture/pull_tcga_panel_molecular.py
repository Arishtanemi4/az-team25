from __future__ import annotations

import argparse
import gzip
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


API = "https://www.cbioportal.org/api"
STUDY_SUFFIX = "_tcga_pan_can_atlas_2018"
PANEL_GENES = [
    "ABL1", "ALK", "ARID1A", "BCR", "BRAF", "BRCA1", "BRCA2", "CDK4",
    "EML4", "ERG", "EWSR1", "EZH2", "FLI1", "KEAP1", "KRAS", "MAPK1",
    "MAT2A", "MDM2", "MTAP", "PARP1", "PRMT5", "RAF1", "SMARCA2",
    "SMARCA4", "STK11", "TMPRSS2",
]


def find_project_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        if (candidate / "data" / "depmap_24Q4" / "architecture" / "gene_master.csv").exists():
            return candidate
    raise FileNotFoundError("Run this script from inside the az-team25 repository.")


def build_session() -> requests.Session:
    retry = Retry(
        total=6,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(("GET", "POST")),
    )
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "az-team25-tcga-panel-audit/1.0"})
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def get_json(session: requests.Session, path: str, **params: Any) -> list[dict[str, Any]] | dict[str, Any]:
    response = session.get(f"{API}{path}", params=params, timeout=180)
    response.raise_for_status()
    return response.json()


def post_json(session: requests.Session, path: str, body: dict[str, Any], **params: Any) -> list[dict[str, Any]]:
    response = session.post(f"{API}{path}", params=params, json=body, timeout=300)
    response.raise_for_status()
    if not response.content:
        return []
    return response.json()


def write_json_gz(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=True, separators=(",", ":"))


def read_json_gz(path: Path) -> Any:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def cached_get(session: requests.Session, path: str, cache_path: Path, refresh: bool = False, **params: Any) -> Any:
    if cache_path.exists() and not refresh:
        return read_json_gz(cache_path)
    value = get_json(session, path, **params)
    write_json_gz(cache_path, value)
    return value


def cached_post(
    session: requests.Session,
    path: str,
    body: dict[str, Any],
    cache_path: Path,
    refresh: bool = False,
    **params: Any,
) -> Any:
    if cache_path.exists() and not refresh:
        return read_json_gz(cache_path)
    value = post_json(session, path, body, **params)
    write_json_gz(cache_path, value)
    return value


def choose_profile(profiles: list[dict[str, Any]], layer: str) -> str | None:
    if layer == "mutation":
        candidates = [p for p in profiles if p.get("molecularAlterationType") == "MUTATION_EXTENDED"]
    elif layer == "cna":
        candidates = [
            p for p in profiles
            if p.get("molecularAlterationType") == "COPY_NUMBER_ALTERATION"
            and p.get("datatype") == "DISCRETE"
        ]
    elif layer == "rna":
        candidates = [
            p for p in profiles
            if p.get("molecularAlterationType") == "MRNA_EXPRESSION"
            and p.get("datatype") == "CONTINUOUS"
            and "RSEM" in p.get("name", "").upper()
        ]
    elif layer == "structural_variant":
        candidates = [p for p in profiles if p.get("molecularAlterationType") == "STRUCTURAL_VARIANT"]
    else:
        raise ValueError(layer)
    if not candidates:
        return None
    candidates.sort(key=lambda p: ("ZSCORE" in p.get("molecularProfileId", "").upper(), p["molecularProfileId"]))
    return candidates[0]["molecularProfileId"]


def choose_sample_list(sample_lists: list[dict[str, Any]], category: str, fallback: str) -> str | None:
    matches = [x for x in sample_lists if x.get("category") == category]
    if matches:
        return sorted(x["sampleListId"] for x in matches)[0]
    fallback_matches = [x for x in sample_lists if x.get("sampleListId") == fallback]
    return fallback_matches[0]["sampleListId"] if fallback_matches else None


def normalize_samples(study_id: str, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "study_id": study_id,
            "sample_id": row.get("sampleId"),
            "patient_id": row.get("patientId"),
            "sample_type": row.get("sampleType"),
            "sample_type_id": row.get("sampleTypeId"),
        }
        for row in rows
    ]


def normalize_numeric(layer: str, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "study_id": row.get("studyId"),
            "sample_id": row.get("sampleId"),
            "patient_id": row.get("patientId"),
            "gene_symbol": (row.get("gene") or {}).get("hugoGeneSymbol"),
            "entrez_gene_id": row.get("entrezGeneId"),
            "value": row.get("value"),
            "molecular_profile_id": row.get("molecularProfileId"),
            "layer": layer,
        }
        for row in rows
    ]


def normalize_mutations(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fields = [
        "studyId", "sampleId", "patientId", "entrezGeneId", "mutationType",
        "proteinChange", "aminoAcidChange", "chromosome", "startPosition",
        "endPosition", "referenceAllele", "variantAllele", "variantType",
        "ncbiBuild", "center", "mutationStatus", "validationStatus",
        "tumorAltCount", "tumorRefCount", "normalAltCount", "normalRefCount",
    ]
    result = []
    for row in rows:
        out = {field: row.get(field) for field in fields}
        out["geneSymbol"] = (row.get("gene") or {}).get("hugoGeneSymbol")
        result.append(out)
    return result


def normalize_structural_variants(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    fields = [
        "studyId", "sampleId", "patientId", "molecularProfileId",
        "site1EntrezGeneId", "site1HugoSymbol", "site1Chromosome", "site1Position",
        "site2EntrezGeneId", "site2HugoSymbol", "site2Chromosome", "site2Position",
        "site2EffectOnFrame", "dnaSupport", "rnaSupport", "tumorPairedEndReadCount",
        "tumorSplitReadCount", "annotation", "breakpointType", "connectionType",
        "eventInfo", "variantClass", "svStatus", "ncbiBuild",
    ]
    return [{field: row.get(field) for field in fields} for row in rows]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true", help="Ignore compressed API response cache.")
    parser.add_argument("--max-studies", type=int, default=None, help="Optional smoke-test limit.")
    args = parser.parse_args()

    root = find_project_root(Path.cwd().resolve())
    base = root / "data" / "external" / "tcga_cbioportal"
    raw = base / "raw"
    derived = base / "derived"
    derived.mkdir(parents=True, exist_ok=True)

    gene_master = pd.read_csv(root / "data" / "depmap_24Q4" / "architecture" / "gene_master.csv", low_memory=False)
    genes = gene_master.set_index("symbol").reindex(PANEL_GENES).reset_index()
    if genes["entrez_id"].isna().any():
        missing = genes.loc[genes["entrez_id"].isna(), "symbol"].tolist()
        raise ValueError(f"Panel genes without Entrez ID: {missing}")
    genes["entrez_id"] = genes["entrez_id"].astype(int)
    entrez_ids = genes["entrez_id"].tolist()
    genes[["symbol", "gene_id", "hgnc_id", "entrez_id", "ensembl_gene_id"]].to_csv(
        derived / "tcga_panel_gene_manifest.csv", index=False
    )

    session = build_session()
    studies_all = cached_get(
        session, "/studies", raw / "api_metadata" / "studies.json.gz",
        args.refresh, projection="SUMMARY", pageSize=10_000_000,
    )
    studies = sorted(
        [study for study in studies_all if study.get("studyId", "").endswith(STUDY_SUFFIX)],
        key=lambda study: study["studyId"],
    )
    if args.max_studies is not None:
        studies = studies[: args.max_studies]

    all_samples: list[dict[str, Any]] = []
    all_mutations: list[dict[str, Any]] = []
    all_cna: list[dict[str, Any]] = []
    all_rna: list[dict[str, Any]] = []
    all_sv: list[dict[str, Any]] = []
    study_manifest: list[dict[str, Any]] = []

    for index, study in enumerate(studies, start=1):
        study_id = study["studyId"]
        print(f"[{index:02d}/{len(studies):02d}] {study_id}", flush=True)
        study_raw = raw / study_id
        profiles = cached_get(
            session, f"/studies/{study_id}/molecular-profiles",
            study_raw / "molecular_profiles.json.gz", args.refresh, projection="SUMMARY",
        )
        sample_lists = cached_get(
            session, f"/studies/{study_id}/sample-lists",
            study_raw / "sample_lists.json.gz", args.refresh, projection="DETAILED",
        )
        samples = cached_get(
            session, f"/studies/{study_id}/samples",
            study_raw / "samples.json.gz", args.refresh,
            projection="DETAILED", pageSize=10_000_000,
        )
        all_samples.extend(normalize_samples(study_id, samples))

        profile_ids = {layer: choose_profile(profiles, layer) for layer in ("mutation", "cna", "rna", "structural_variant")}
        sample_list_ids = {
            "mutation": choose_sample_list(sample_lists, "all_cases_with_mutation_data", f"{study_id}_sequenced"),
            "cna": choose_sample_list(sample_lists, "all_cases_with_cna_data", f"{study_id}_cna"),
            "rna": choose_sample_list(sample_lists, "all_cases_with_mrna_rnaseq_data", f"{study_id}_rna_seq_v2_mrna"),
            "structural_variant": choose_sample_list(sample_lists, "all_cases_with_sv_data", f"{study_id}_sv"),
        }
        list_sizes = {row["sampleListId"]: len(row.get("sampleIds", [])) for row in sample_lists}

        mutation_rows: list[dict[str, Any]] = []
        if profile_ids["mutation"]:
            mutation_rows = cached_post(
                session, "/mutations/fetch",
                {"molecularProfileIds": [profile_ids["mutation"]], "entrezGeneIds": entrez_ids},
                study_raw / "panel_mutations.json.gz", args.refresh,
                projection="DETAILED", pageSize=10_000_000,
            )
            all_mutations.extend(normalize_mutations(mutation_rows))

        cna_rows: list[dict[str, Any]] = []
        if profile_ids["cna"] and sample_list_ids["cna"]:
            cna_rows = cached_post(
                session, f"/molecular-profiles/{profile_ids['cna']}/molecular-data/fetch",
                {"sampleListId": sample_list_ids["cna"], "entrezGeneIds": entrez_ids},
                study_raw / "panel_cna.json.gz", args.refresh, projection="DETAILED",
            )
            all_cna.extend(normalize_numeric("discrete_cna", cna_rows))

        rna_rows: list[dict[str, Any]] = []
        if profile_ids["rna"] and sample_list_ids["rna"]:
            rna_rows = cached_post(
                session, f"/molecular-profiles/{profile_ids['rna']}/molecular-data/fetch",
                {"sampleListId": sample_list_ids["rna"], "entrezGeneIds": entrez_ids},
                study_raw / "panel_rna.json.gz", args.refresh, projection="DETAILED",
            )
            all_rna.extend(normalize_numeric("rsem_rna", rna_rows))

        sv_rows: list[dict[str, Any]] = []
        if profile_ids["structural_variant"]:
            sv_rows = cached_post(
                session, "/structural-variant/fetch",
                {"molecularProfileIds": [profile_ids["structural_variant"]], "entrezGeneIds": entrez_ids},
                study_raw / "panel_structural_variants.json.gz", args.refresh,
            )
            all_sv.extend(normalize_structural_variants(sv_rows))

        study_manifest.append({
            "study_id": study_id,
            "study_name": study.get("name"),
            "all_sample_count_reported": study.get("allSampleCount"),
            "api_sample_rows": len(samples),
            "mutation_profile_id": profile_ids["mutation"],
            "mutation_sample_list_id": sample_list_ids["mutation"],
            "mutation_measured_samples": list_sizes.get(sample_list_ids["mutation"], 0),
            "panel_mutation_events": len(mutation_rows),
            "cna_profile_id": profile_ids["cna"],
            "cna_sample_list_id": sample_list_ids["cna"],
            "cna_measured_samples": list_sizes.get(sample_list_ids["cna"], 0),
            "panel_cna_rows": len(cna_rows),
            "rna_profile_id": profile_ids["rna"],
            "rna_sample_list_id": sample_list_ids["rna"],
            "rna_measured_samples": list_sizes.get(sample_list_ids["rna"], 0),
            "panel_rna_rows": len(rna_rows),
            "sv_profile_id": profile_ids["structural_variant"],
            "sv_sample_list_id": sample_list_ids["structural_variant"],
            "sv_measured_samples": list_sizes.get(sample_list_ids["structural_variant"], 0),
            "panel_structural_variant_events": len(sv_rows),
        })
        time.sleep(0.1)

    samples_df = pd.DataFrame(all_samples).drop_duplicates()
    mutations_df = pd.DataFrame(all_mutations)
    cna_df = pd.DataFrame(all_cna)
    rna_df = pd.DataFrame(all_rna)
    sv_df = pd.DataFrame(all_sv)
    studies_df = pd.DataFrame(study_manifest)

    samples_df.to_csv(derived / "tcga_pancan_panel_samples.csv", index=False)
    mutations_df.to_csv(derived / "tcga_pancan_panel_mutations.csv", index=False)
    cna_df.to_csv(derived / "tcga_pancan_panel_cna.csv", index=False)
    rna_df.to_csv(derived / "tcga_pancan_panel_rna.csv", index=False)
    sv_df.to_csv(derived / "tcga_pancan_panel_structural_variants.csv", index=False)
    studies_df.to_csv(derived / "tcga_pancan_panel_study_manifest.csv", index=False)

    retrieved = datetime.now(timezone.utc).isoformat()
    manifest = {
        "generated_utc": retrieved,
        "source": "cBioPortal public API",
        "api_base": API,
        "study_family": "TCGA PanCancer Atlas 2018",
        "study_suffix": STUDY_SUFFIX,
        "scope": "Current 26-gene EDA panel only; not genome-wide TCGA data",
        "reference_genome": "Study metadata reports hg19/GRCh37",
        "genes": PANEL_GENES,
        "entrez_gene_ids": entrez_ids,
        "studies": len(studies_df),
        "derived_rows": {
            "samples": len(samples_df),
            "mutation_events": len(mutations_df),
            "cna_sample_gene": len(cna_df),
            "rna_sample_gene": len(rna_df),
            "structural_variant_events": len(sv_df),
        },
        "files": {
            "gene_manifest": "derived/tcga_panel_gene_manifest.csv",
            "study_manifest": "derived/tcga_pancan_panel_study_manifest.csv",
            "samples": "derived/tcga_pancan_panel_samples.csv",
            "mutations": "derived/tcga_pancan_panel_mutations.csv",
            "cna": "derived/tcga_pancan_panel_cna.csv",
            "rna": "derived/tcga_pancan_panel_rna.csv",
            "structural_variants": "derived/tcga_pancan_panel_structural_variants.csv",
        },
        "important_limitations": [
            "cBioPortal PanCancer Atlas processing is not identical to current GDC harmonized processing.",
            "The acquisition is panel-scoped and cannot support discovery outside the listed genes.",
            "Mutation and structural-variant tables contain observed events; measured negatives require the layer sample-list denominator.",
            "Discrete CNA values are GISTIC calls (-2, -1, 0, 1, 2), not DepMap absolute copy number.",
            "RSEM RNA values and DepMap expression values require within-source normalization before comparison.",
        ],
    }
    with open(base / "tcga_panel_acquisition_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)

    external_manifest_path = root / "data" / "external" / "external_layer_manifest.json"
    if external_manifest_path.exists():
        with open(external_manifest_path, encoding="utf-8") as handle:
            external_manifest = json.load(handle)
        external_manifest["generated_utc"] = retrieved
        external_manifest.setdefault("downloads", {})["tcga_cbioportal_panel"] = {
            "url": API,
            "target": str(base),
            "release": "TCGA PanCancer Atlas 2018 studies; cBioPortal public API state at retrieval time",
            "status": "downloaded",
            "scope": "26-gene EDA panel, 32 studies, mutation/CNA/RNA/structural variants",
        }
        external_manifest.setdefault("results", {})["tcga_cbioportal_panel"] = manifest["derived_rows"] | {
            "studies": len(studies_df), "genes": len(PANEL_GENES)
        }
        with open(external_manifest_path, "w", encoding="utf-8") as handle:
            json.dump(external_manifest, handle, indent=2)

    print(json.dumps(manifest["derived_rows"] | {"studies": len(studies_df), "genes": len(PANEL_GENES)}, indent=2))


if __name__ == "__main__":
    main()
