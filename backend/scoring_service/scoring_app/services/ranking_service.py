"""Wraps scoring/score.py::score_panel -- the one entry point a researcher's inclusion/exclusion
gene query actually goes through (docs/reference/ALGORITHM_SPEC.md). This class validates tokens
before SSE starts, threads request filters and top-N into the native scorer, assembles a complete
export, and streams progress checkpoints. It does not re-implement, alter, or duplicate scoring
math.
"""

import copy
import os
import queue
import threading
import time

from scoring_app.lib import scoring_path  # noqa: F401 -- import order matters: puts scoring/ on
                                    # sys.path before the `import score`/`import export` below can succeed.
import export as scoring_export
import io_utils as scoring_io_utils
import score as scoring_score

from scoring_app.config import (
    CELL_LINES_CSV,
    DATA_DIR,
    DEFAULT_TOP_K,
    ESSENTIALITY_CONSTANTS_PATH,
    EXTENDED_CONSTANTS_PATH,
    GENE_REFERENCE_CSV,
    GENE_ROLE_PATH,
    RNA_CONSTANTS_PATH,
    UNRESOLVED_SYMBOLS_PATH,
)

# Sentinel pushed onto rank_stream's queue to signal the background thread has finished (either
# with a result or an exception) -- a plain None can't be used since a stage label is also a
# plain Python object and could theoretically collide.
_STREAM_DONE = object()

# 7 checkpoints inside score_panel plus 1 for this class's own export step.
TOTAL_STAGES = len(scoring_score.STAGE_LABELS) + 1

FINAL_STAGE_LABEL = "Assembling the final evidence report"

# Mirrors scoring/score.py::load_tables_for_query's own file list (also independently redeclared
# in scoring/tests/test_parquet_read_path.py::LAYER_FILES) -- kept as a plain inline list here
# too rather than a new shared constant, since 8 filenames used in two or three places doesn't
# justify another module.
_EVIDENCE_TABLES = [
    "expression_rna.csv", "expression_rna_hpa.csv", "expression_rna_geo.csv", "protein.csv",
    "dependency.csv", "copy_number.csv", "mutations.csv", "fusions.csv",
]


class RankingService:
    """Loads the small reference tables once at startup (score_panel itself streams the large
    per-gene evidence tables fresh on every call, since those are tens of millions of rows and
    must never be cached whole in memory)."""

    def __init__(self, result_callback=None):
        self.gene_reference_df = scoring_io_utils.read_table(str(GENE_REFERENCE_CSV))
        self.cell_lines_df = scoring_io_utils.read_table(str(CELL_LINES_CSV))
        # Display metadata only, joined on the same ModelID spine every scoring layer already
        # joins on -- never fed into scoring/ (V7-1c, CONSTRAINTS.md F1's sibling display rule).
        self.cell_line_names = dict(
            zip(self.cell_lines_df["ModelID"], self.cell_lines_df["cell_line_name"])
        )
        self._check_parquet_mirrors()
        # S09's optional extension hook receives only a fully assembled native export. Keeping
        # it optional preserves the disabled service's response contract and scoring maths.
        self.result_callback = result_callback

    def _check_parquet_mirrors(self) -> None:
        """Warns once at startup if any big evidence table is missing its `.parquet` mirror --
        without this, that table's every /rank call would silently fall back to a slow chunked
        CSV scan instead of the fast path, and nothing else in the request/response cycle would
        show it."""
        for filename in _EVIDENCE_TABLES:
            csv_path = os.path.join(DATA_DIR, filename)
            parquet_path = csv_path.replace(".csv", ".parquet")
            if not os.path.exists(parquet_path):
                print(
                    f"[RankingService] WARNING: no Parquet mirror for {csv_path} -- "
                    "this table will use the slow chunked-CSV fallback on every query. "
                    "Run preprocessing/preprocess.py to generate it."
                )

    def check_genes(self, inclusion_genes: list[str], exclusion_genes: list[str]) -> None:
        """Runs score_panel's own gene-resolution check up front, before an SSE stream begins,
        so an ambiguous/unresolved gene token still surfaces as a normal HTTP 422 -- once a
        streaming response's headers are sent (status 200), the status code can no longer
        change. Raises ValueError, same as score_panel's own first stage would. A canonical gene
        cannot appear in both roles: checking resolved IDs catches a symbol/Ensembl spelling mix
        before the stream sends its HTTP 200 headers."""
        if not inclusion_genes:
            raise ValueError(
                "At least one inclusion gene is required; exclusion genes are optional criteria."
            )
        resolved = scoring_score.resolve_genes_or_raise(
            inclusion_genes + exclusion_genes, self.gene_reference_df
        )
        overlapping_ids = sorted(
            set(resolved[token] for token in inclusion_genes)
            & set(resolved[token] for token in exclusion_genes)
        )
        if overlapping_ids:
            raise ValueError(
                "A gene cannot be both an inclusion and exclusion criterion: "
                f"{overlapping_ids}."
            )

    def rank_stream(
        self,
        inclusion_genes: list[str],
        exclusion_genes: list[str],
        lineage: str | None = None,
        primary_disease: str | None = None,
        exclude_problematic: bool = False,
        msi_high: bool | None = None,
        ploidy_min: float | None = None,
        ploidy_max: float | None = None,
        require_metabolomics: bool = False,
        require_mirna: bool = False,
        top_k: int = DEFAULT_TOP_K,
    ):
        """Generator: yields ("stage", label) after each of score_panel's 7 checkpoints plus one
        more for this method's own export step, then exactly one ("result", payload) with
        the same dict rank() returns. Runs score_panel on a background thread so its on_stage
        callback can push labels onto a queue this generator drains in real time -- a plain
        synchronous callback could only report progress after score_panel finished entirely,
        which would defeat the point of streaming a multi-minute computation."""
        stage_queue: queue.Queue = queue.Queue()
        outcome_holder: dict = {}
        start_time = time.monotonic()

        def logged_on_stage(label):
            # TEMP DIAGNOSTIC (see plan): times each score_panel checkpoint to locate where a
            # /rank request stalls. Remove once the hang is root-caused and fixed.
            print(f"[ranking_service] +{time.monotonic() - start_time:.1f}s stage: {label}")
            stage_queue.put(label)

        def run_score_panel():
            try:
                outcome_holder["result"] = scoring_score.score_panel(
                    inclusion_genes,
                    exclusion_genes,
                    self.gene_reference_df,
                    self.cell_lines_df,
                    lineage=lineage,
                    primary_disease=primary_disease,
                    exclude_problematic=exclude_problematic,
                    msi_high=msi_high,
                    ploidy_min=ploidy_min,
                    ploidy_max=ploidy_max,
                    require_metabolomics=require_metabolomics,
                    require_mirna=require_mirna,
                    top_n=top_k,
                    data_dir=str(DATA_DIR),
                    rna_constants_path=str(RNA_CONSTANTS_PATH),
                    extended_constants_path=str(EXTENDED_CONSTANTS_PATH),
                    essentiality_constants_path=str(ESSENTIALITY_CONSTANTS_PATH),
                    gene_role_path=str(GENE_ROLE_PATH),
                    unresolved_symbols_path=str(UNRESOLVED_SYMBOLS_PATH),
                    on_stage=logged_on_stage,
                )
            except Exception as exc:  # noqa: BLE001 -- re-raised on the generator's own thread below
                outcome_holder["error"] = exc
            finally:
                stage_queue.put(_STREAM_DONE)

        worker = threading.Thread(target=run_score_panel, daemon=True)
        worker.start()
        while True:
            item = stage_queue.get()
            if item is _STREAM_DONE:
                break
            yield "stage", item
        worker.join()

        if "error" in outcome_holder:
            raise outcome_holder["error"]

        outcome = outcome_holder["result"]
        print(f"[ranking_service] +{time.monotonic() - start_time:.1f}s stage: {FINAL_STAGE_LABEL}")
        yield "stage", FINAL_STAGE_LABEL

        export = scoring_export.export_query_result(
            outcome,
            inclusion_genes,
            exclusion_genes,
            filters={
                "lineage": lineage,
                "primary_disease": primary_disease,
                "exclude_problematic": exclude_problematic,
                "msi_high": msi_high,
                "ploidy_min": ploidy_min,
                "ploidy_max": ploidy_max,
                "require_metabolomics": require_metabolomics,
                "require_mirna": require_mirna,
            },
        )
        # Compatibility field for the existing UI: unlike the old top-ten-derived value, this is
        # now the complete partition count assembled by the export contract.
        export["total_candidates_scored"] = export["candidate_counts"]["evaluated"]
        # V7-1c: attach each line's cell-line name here, not in scoring/export.py -- this is
        # display metadata, never a scored input. A ModelID missing from cell_lines.csv (or a
        # NaN name) yields None, which sanitize_for_json below turns into a safe JSON null.
        for partition_key in (
            "ranked_cell_lines",
            "ranked_beyond_top_n",
            "low_confidence_lines",
            "insufficient_evidence_lines",
            "disqualified_lines",
        ):
            for line in export[partition_key]:
                line["cell_line_name"] = self.cell_line_names.get(line["model_id"])
        callback = getattr(self, "result_callback", None)
        if callback is not None:
            print(f"[ranking_service] +{time.monotonic() - start_time:.1f}s result_callback: start")
            try:
                # A callback gets copies: extension persistence cannot mutate the native result
                # that will be returned by this service, even accidentally.
                research_result = callback(copy.deepcopy(export), self.cell_lines_df.copy(deep=True))
                export["research_query_id"] = research_result["research_query_id"]
                export["research_status"] = research_result.get("research_status", "ready")
            except Exception:  # noqa: BLE001 -- enrichment must never invalidate a valid rank result
                export["research_query_id"] = None
                export["research_status"] = "unavailable"
                export["research_reason"] = "Research snapshot creation failed; native ranking remains available."
            print(f"[ranking_service] +{time.monotonic() - start_time:.1f}s result_callback: done")
        print(f"[ranking_service] +{time.monotonic() - start_time:.1f}s stage: result (final SSE event)")
        yield "result", scoring_export.sanitize_for_json(export)

    def rank(
        self,
        inclusion_genes: list[str],
        exclusion_genes: list[str],
        lineage: str | None = None,
        primary_disease: str | None = None,
        exclude_problematic: bool = False,
        msi_high: bool | None = None,
        ploidy_min: float | None = None,
        ploidy_max: float | None = None,
        require_metabolomics: bool = False,
        require_mirna: bool = False,
        top_k: int = DEFAULT_TOP_K,
    ) -> dict:
        """Blocking convenience wrapper over rank_stream() for any non-SSE caller (scripts,
        notebooks) -- drains every stage event and returns only the final result, so there is
        one implementation of the pipeline glue, not two that could drift. Raises ValueError on
        an ambiguous or unresolved gene token, same as before."""
        for kind, payload in self.rank_stream(
            inclusion_genes,
            exclusion_genes,
            lineage=lineage,
            primary_disease=primary_disease,
            exclude_problematic=exclude_problematic,
            msi_high=msi_high,
            ploidy_min=ploidy_min,
            ploidy_max=ploidy_max,
            require_metabolomics=require_metabolomics,
            require_mirna=require_mirna,
            top_k=top_k,
        ):
            if kind == "result":
                return payload
        raise RuntimeError("rank_stream() ended without yielding a result")  # unreachable by contract
