"""Wraps rag/'s five entry points -- narrator.narrate, methodology_agent.answer,
query_expansion_agent.expand, literature_agent.find_evidence, and citations.build_registry
(built once here, not per request, and passed into the other three as registry=). This class
adds only what an API layer needs on top: retrieving the C1 doc context narrate() itself expects
the caller to supply, and the one-time registry cache -- it does not re-implement, alter, or
duplicate any of rag/'s own generation or verification logic.
"""

from pathlib import Path

from rag_app.lib import rag_path  # noqa: F401 -- import order matters: puts rag/ on sys.path
                                    # before the `import citations`/`import narrator` below can succeed.

import citations
import display_sanitiser
import literature_agent
import methodology_agent
import narrator
import query_expansion_agent

_REPO_ROOT = Path(__file__).resolve().parents[4]
C1_INDEX_DIR = _REPO_ROOT / "data" / "processed" / "rag" / "c1_index"

_c1_index_cache = None


def _search_c1(query, top_k):
    """Retrieved C1 (project-doc) chunks for narrator.narrate()'s context_chunks -- narrate()
    does no retrieval of its own, the caller must supply it. Same graceful-degrade convention
    rag/methodology_agent.py::search_docs already uses: an unbuilt index (this machine hasn't run
    rag/build_c1_index.py yet) returns no context rather than raising."""
    global _c1_index_cache
    if not (C1_INDEX_DIR / "chunks.jsonl").exists():
        return []
    if _c1_index_cache is None:
        import retrieval

        _c1_index_cache = retrieval.HybridIndex.load(C1_INDEX_DIR)
    return _c1_index_cache.search(query, top_k=top_k)


# These four partitions are unbounded (they can hold every candidate that didn't make the
# ranked top-N -- up to the full cell-line panel), unlike ranked_cell_lines which is already
# capped by the query's top_n. Dumping them verbatim into the narration prompt is what blows the
# LLM's context window on broad queries. The record's own candidate_counts field already carries
# the authoritative counts for each of these, so nothing the narrative needs is lost by
# summarizing them here instead.
_UNBOUNDED_PARTITIONS = (
    "ranked_beyond_top_n",
    "low_confidence_lines",
    "insufficient_evidence_lines",
    "disqualified_lines",
)


def _bound_evidence_record(evidence_record):
    """Returns a copy of evidence_record safe to hand to narrator.narrate(): the unbounded
    partitions are replaced with a count summary instead of their full line-level contents."""
    bounded = dict(evidence_record)
    for key in _UNBOUNDED_PARTITIONS:
        lines = bounded.get(key)
        if isinstance(lines, list):
            bounded[key] = f"[{len(lines)} lines omitted -- see candidate_counts for the summary]"
    return bounded


def _narration_query(evidence_record):
    """A plain retrieval query built from the query's own gene symbols -- narrate()'s own prompt
    says the retrieved context is 'for methodology explanations and citations only'."""
    genes = evidence_record.get("query", {}).get("inclusion_genes", []) + evidence_record.get(
        "query", {}
    ).get("exclusion_genes", [])
    return " ".join(genes) + " methodology"


class RagService:
    """Builds the citation registry once at process startup -- citations.build_registry() is
    cheap per its own docstring (a few dozen small markdown files), but there is no reason to
    rebuild it on every request when every rag/ call in a single process can share one."""

    def __init__(self):
        self.registry = citations.build_registry()

    def narrate_result(self, evidence_record, top_k=5):
        query = _narration_query(evidence_record)
        context_chunks = _search_c1(query, top_k)
        bounded_record = _bound_evidence_record(evidence_record)
        response = narrator.narrate(bounded_record, context_chunks=context_chunks, registry=self.registry)
        return display_sanitiser.sanitise_narration_response(response)

    def answer_methodology_question(self, question):
        response = methodology_agent.answer(question, registry=self.registry)
        return display_sanitiser.sanitise_methodology_response(response)

    def expand_query(self, inclusion_genes):
        response = query_expansion_agent.expand(inclusion_genes, registry=self.registry)
        return display_sanitiser.sanitise_expansion_response(response)

    def find_literature(self, question):
        response = literature_agent.find_evidence(question)
        return display_sanitiser.sanitise_literature_response(response)
