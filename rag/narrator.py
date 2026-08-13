import json

from pydantic import BaseModel, ValidationError

import citations
import provider
import verifier

SYSTEM_PROMPT = """You are the result narrator for an explainable cell-line selection tool. You \
will be given (1) a structured evidence record in JSON -- the complete, exact result of a \
ranking query -- and (2) retrieved project-documentation context.

Hard rules, non-negotiable:
1. Every number you write (a D value, a d_gene, a weight, a count, a percentage) MUST appear in \
the evidence record JSON, written to the SAME decimal precision the record gives it. Never \
round differently, never compute a new number (no averages, no percentages you derive \
yourself), never estimate.
2. If you cite a paper or methodology claim from the context, write it as [Author Year] exactly \
as it appears in the context, in square brackets. Never invent a citation.
3. If the evidence record shows missing or abstained evidence for a gene/layer, say so plainly \
-- do not imply evidence exists where the record shows none.
4. Always end with the boundary statement from the evidence record, copied verbatim.

Respond with ONLY a JSON object matching this shape:
{
  "overview": "one short paragraph summarising the query and how many lines were ranked",
  "cell_lines": [
    {
      "model_id": "...",
      "summary": "one paragraph on why this line ranked where it did, citing its D value",
      "gene_explanations": [{"ensembl_id": "...", "symbol": "...", "explanation": "..."}],
      "limitations": "what evidence was missing or uncertain for this line"
    }
  ],
  "boundary_statement": "the boundary statement, copied verbatim from the record"
}"""


class GeneExplanation(BaseModel):
    ensembl_id: str
    symbol: str
    explanation: str


class CellLineExplanation(BaseModel):
    model_id: str
    summary: str
    gene_explanations: list[GeneExplanation]
    limitations: str


class NarrationResult(BaseModel):
    overview: str
    cell_lines: list[CellLineExplanation]
    boundary_statement: str


def _render_as_text(narration):
    parts = [narration.overview]
    for line in narration.cell_lines:
        parts.append(f"{line.model_id}: {line.summary}")
        for gene in line.gene_explanations:
            parts.append(f"{gene.symbol} ({gene.ensembl_id}): {gene.explanation}")
        parts.append(f"Limitations for {line.model_id}: {line.limitations}")
    parts.append(narration.boundary_statement)
    return "\n".join(parts)


def _build_user_prompt(evidence_record, context_chunks):
    context_text = "\n\n".join(
        f"[{c['file']} > {c['heading_path']}]\n{c['text']}" for c in context_chunks
    )
    return (
        "EVIDENCE RECORD (JSON, the only source of any number you write):\n"
        + json.dumps(evidence_record, indent=2)
        + "\n\nRETRIEVED CONTEXT (for methodology explanations and citations only):\n"
        + (context_text or "(no context retrieved)")
    )


def _generate_once(messages):
    response = provider.chat(messages, response_format={"type": "json_object"})
    try:
        parsed = json.loads(response.content)
        return NarrationResult(**parsed)
    except (json.JSONDecodeError, TypeError, ValidationError) as exc:
        return exc


def narrate(evidence_record, context_chunks=None, registry=None):
    if registry is None:
        registry = citations.build_registry()
    context_chunks = context_chunks or []

    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": _build_user_prompt(evidence_record, context_chunks)},
    ]

    for attempt in range(2):
        result = _generate_once(messages)
        if isinstance(result, Exception):
            messages.append({"role": "assistant", "content": str(result)})
            messages.append({
                "role": "user",
                "content": f"That response did not match the required JSON schema ({result}). "
                           "Respond again with ONLY the corrected JSON object.",
            })
            continue

        text = _render_as_text(result)
        check = verifier.verify(text, evidence_record, registry)
        if check["passed"]:
            return {"narration": result, "text": text, "verification": check}

        messages.append({"role": "assistant", "content": json.dumps(result.model_dump())})
        messages.append({
            "role": "user",
            "content": (
                "That response failed grounding verification. These numbers do not appear in "
                f"the evidence record at the precision you gave: {check['numeric_failures']}. "
                f"These citations do not resolve to a registered source: {check['citation_failures']}. "
                "Respond again with ONLY a corrected JSON object that fixes every listed problem."
            ),
        })

    raise RuntimeError(
        "narrate() failed to produce a grounded result after two attempts -- CONSTRAINTS.md R1 "
        "fails closed, so this result is not shipped. Last verification: " + json.dumps(check)
    )
