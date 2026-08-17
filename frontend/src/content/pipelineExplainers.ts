// Static, per-step formula + biological-rationale text for the (i) info icons in
// CellLineDetail.tsx. This text never changes per query, so it is authored once here rather
// than served by an API. Sourced verbatim in spirit from docs/reference/ALGORITHM_SPEC.md --
// if that spec changes, update this file to match.

export interface PipelineStepExplainer {
  step: number;
  title: string;
  formula: string;
  biology: string;
}

export const PIPELINE_STEPS: PipelineStepExplainer[] = [
  {
    step: 1,
    title: "Per-gene desirability",
    formula:
      "Inclusion gene: d = 0 below floor L, d = ((y-L)/(T-L))^r between L and T, d = 1 above T. " +
      "Exclusion gene: the mirror image (d = 1 below L, d = 0 above T). r = 1 (linear), fixed, " +
      "never user-tunable. Mutation and fusion have no continuous value, so they use a " +
      "deterministic rule table instead of this formula.",
    biology:
      "L is the detection floor -- below it, the gene is not credibly expressed. T is the " +
      "confident high-expression reference -- at or above it, the criterion is fully met and " +
      "more expression adds nothing. L and T come from the 10th/90th percentile of measured " +
      "values, per-lineage when enough lines exist, otherwise a global default (Derringer & " +
      "Suich 1980).",
  },
  {
    step: 2,
    title: "Within-gene layer combination",
    formula:
      "d_gene = sum(weight_l * d_l) / sum(weight_l), over whichever of the six evidence layers " +
      "(RNA 4/15 ≈ 0.267, dependency 4/15 ≈ 0.267, protein 2/15 ≈ 0.133, mutation 2/15 ≈ 0.133, " +
      "copy number 2/15 ≈ 0.133, fusion 1/15 ≈ 0.067) are actually present for this gene. A " +
      "fusion layer's effective weight is further discounted (x1.0 / x0.5 / x0.25) by how " +
      "confident/in-frame its events are.",
    biology:
      "This is a weighted arithmetic mean -- deliberately compensatory. An assay that was " +
      "never run must not veto a gene that another layer strongly supports; the veto belongs " +
      "one step later. Weights are fixed and published, never tuned per query.",
  },
  {
    step: 3,
    title: "Across-gene combination (the veto)",
    formula:
      "D = exp( sum(w_i * ln d_i) / sum(w_i) ) -- a weighted geometric mean, with a hard rule: " +
      "if any gene's d_i = 0 (and its weight > 0), then D = 0 for the whole cell line, no " +
      "matter how well every other gene scored.",
    biology:
      "This non-compensatory veto is the entire reason a geometric mean was chosen over an " +
      "arithmetic one: a confidently expressed exclusion gene, or an inclusion gene nowhere " +
      "near its floor, must disqualify a cell line outright rather than just docking it a few " +
      "points (Derringer & Suich 1980; veto pattern also used in Segall 2012).",
  },
  {
    step: 4,
    title: "Correlation discount (CAMERA)",
    formula:
      "VIF_i = 1 + (m-1) * rho_bar_i, w_i = 1/VIF_i, m_eff = sum(w_i). rho_bar_i is gene i's " +
      "mean Spearman correlation with the other query genes, measured across the reference " +
      "panel and clamped to >= 0. These w_i are the weights layer 3 actually uses.",
    biology:
      "Two co-regulated genes are largely one biological signal counted twice. This discounts " +
      "a gene's influence the more it overlaps with the rest of the query, without ever letting " +
      "it neutralise a veto (w_i is always > 0). m_eff -- the effective number of independent " +
      "gene-evidence pieces -- is reported because CAMERA was validated on 15-580 gene sets, " +
      "and a typical 5-10 gene query is smaller than its demonstrated range (Wu & Smyth 2012).",
  },
];

export const BOUNDARY_STATEMENT_FALLBACK =
  "High confidence means the available evidence is consistent. It does not mean the result " +
  "has been experimentally proven.";
