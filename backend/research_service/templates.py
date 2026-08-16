"""Plain-English building blocks for S04's measured-context packets
(extensions/plans/CONTRACTS.md C3). Disclosure/caveat text only -- describing what a native
field means or where its scope ends, never a new interpretive claim (C3/S04's own prohibition:
no 'pathway activated', 'wild type', 'drug sensitive', 'best tumour model', 'confirmed absence').

Where the exact wording already exists as a native disclosure (scoring/explain.py's own
constants -- already surfaced inline in that gene's `narrative` string), it is re-exported here
verbatim rather than paraphrased, so a structured caveat entry and the prose narrative that
already carries the same information never drift into two different wordings of the same fact.
Only the three caveats with no existing native counterpart are new text, each traceable to a
specific C3 sentence.
"""

import sys
from pathlib import Path

_SCORING_DIR = str(Path(__file__).resolve().parents[2] / "scoring")
if _SCORING_DIR not in sys.path:
    sys.path.insert(0, _SCORING_DIR)

import explain as scoring_explain  # noqa: E402 -- sys.path must be set up first

# Re-exported verbatim -- see scoring/explain.py for the source of truth. Do not re-word these
# here; edit scoring/explain.py and this reference follows automatically.
PAN_ESSENTIAL_NOTE = scoring_explain.PAN_ESSENTIAL_NOTE
GENE_ROLE_UNKNOWN_NOTE = scoring_explain.GENE_ROLE_UNKNOWN_NOTE
PROTEIN_NON_DETECTED_NOTE = scoring_explain.PROTEIN_NON_DETECTED_NOTE
PROTEIN_MIXED_SCALE_NOTE = scoring_explain.PROTEIN_MIXED_SCALE_NOTE
FUSION_STATUS_UNKNOWN_NOTE = scoring_explain.FUSION_STATUS_UNKNOWN
RNA_ONLY_CORRELATION_NOTE = scoring_explain.RNA_ONLY_CORRELATION_NOTE

# New here -- each traces to one C3 sentence with no existing native counterpart.

# C3: "Frozen CN is relative copy number, diploid ~1, not AVI absolute copy number or GISTIC."
COPY_NUMBER_RELATIVE_NOT_ABSOLUTE_NOTE = (
    "Copy number here is relative copy number (diploid is approximately 1), not absolute copy "
    "number and not a GISTIC amplification/deletion call. Any amplification/deletion wording "
    "uses this project's own native thresholds, never a copied absolute-CN or GISTIC threshold."
)

# C3: "Assay-wide measured_absent means no recorded event, not proven gene-level wild type/
# callable absence."
MEASURED_ABSENT_NOT_CONFIRMED_ABSENCE_NOTE = (
    "No event was recorded for this gene in this assay. That reflects the assay's own scope "
    "(what it looked for and where), not a proven gene-level wild-type call or a confirmed "
    "biological absence -- a different, stronger claim this evidence does not make."
)

# C3: "Display first 20 events per gene with counts, export full available rows."
EVENTS_TRUNCATED_NOTE = (
    "Showing the first {displayed} of {total} recorded events for this gene. All {total} are "
    "available in the exported record; none are discarded, only not displayed by default."
)

# Plain-English unit descriptions for the six scored layers (ALGORITHM_SPEC.md's own stated
# units/notation, restated for a non-specialist reader -- not a new definition).
LAYER_UNIT_DESCRIPTIONS = {
    "rna": "RNA expression, log2(TPM+1) -- DepMap RNA-seq, higher means more transcript detected.",
    "protein": (
        "Protein abundance, z-score standardised against the measured panel (panel-relative "
        "by construction -- see the mixed-scale note when this layer scores)."
    ),
    "dependency": (
        "CRISPR knockout gene-effect (Chronos score) -- more negative means the line depends "
        "more strongly on this gene; the sign is flipped before scoring so higher desirability "
        "still means stronger evidence."
    ),
    "copy_number": (
        "Relative copy number, diploid is approximately 1 -- see the relative-not-absolute "
        "note; never absolute copy number or a GISTIC call."
    ),
    "mutation": (
        "Categorical mutation classification (e.g. driver/hotspot, VEP impact tier, or "
        "sequenced-and-clean) -- not a continuous measurement."
    ),
    "fusion": (
        "Count of RNA-seq-detected fusion events involving this gene, with a confidence/"
        "reading-frame modifier on the layer's weight, never on this count."
    ),
}
