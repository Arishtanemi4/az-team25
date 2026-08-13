import re

# Ordered, most-specific-first: a named document gets its tailored phrase before the generic
# fallbacks run. Each name's repo-relative prefix ("docs/plan/", "docs/reference/", "docs/") is
# optional in the pattern, since this codebase's own comments cite these documents both ways
# (e.g. "CONSTRAINTS.md S1" and "docs/PROJECT_ARCHITECTURE.md SS0"), and generated prose could
# echo either form.
_NAMED_REPLACEMENTS = [
    (re.compile(r"(?:docs/)?PROJECT_ARCHITECTURE\.md"), "the project's data architecture"),
    (re.compile(r"(?:docs/reference/)?ALGORITHM_SPEC\.md"), "the algorithm specification"),
    (re.compile(r"(?:docs/plan/)?PARAMETERS\.md"), "the parameter register"),
    (re.compile(r"(?:docs/plan/)?METHOD_DECISION\.md"), "the method decision record"),
    (
        re.compile(r"(?:docs/plan/)?(?:CONSTRAINTS|STATUS|PLAN|EXECUTE)\.md"),
        "the project's design record",
    ),
    (re.compile(r"(?:docs/plan/)?SENSITIVITY\.md"), "the sensitivity analysis"),
    (re.compile(r"(?:docs/plan/)?ESTIMAND\.md"), "the estimand definition"),
    (re.compile(r"(?:docs/)?AZ_REQUIREMENTS\.md"), "the client requirements"),
]

# The project-guide wrapper convention -- every directory's own "_.md" -- matched as a standalone
# token so it never eats part of a longer word.
_PROJECT_GUIDE_RE = re.compile(r"(?<![\w/])_\.md(?![\w])")

# Anything else under docs/ (a bare path, or a *.md this list doesn't name individually).
_BARE_DOCS_PATH_RE = re.compile(r"\bdocs/[\w./-]*")

# Fallback: any other repo-internal *.md reference not already caught above.
_OTHER_MD_RE = re.compile(r"\b[A-Za-z_][\w-]*\.md\b")


def sanitise_text(text):
    if not text:
        return text
    for pattern, replacement in _NAMED_REPLACEMENTS:
        text = pattern.sub(replacement, text)
    text = _PROJECT_GUIDE_RE.sub("the project guide", text)
    text = _BARE_DOCS_PATH_RE.sub("the project documentation", text)
    text = _OTHER_MD_RE.sub("the project documentation", text)
    return text


def sanitise_narration_response(response):
    narration = response.get("narration")
    if narration is not None:
        narration.overview = sanitise_text(narration.overview)
        for line in narration.cell_lines:
            line.summary = sanitise_text(line.summary)
            line.limitations = sanitise_text(line.limitations)
            for gene in line.gene_explanations:
                gene.explanation = sanitise_text(gene.explanation)
        narration.boundary_statement = sanitise_text(narration.boundary_statement)
    if response.get("text") is not None:
        response["text"] = sanitise_text(response["text"])
    return response


def sanitise_methodology_response(response):
    if response.get("answer") is not None:
        response["answer"] = sanitise_text(response["answer"])
    if response.get("context") is not None:
        response["context"] = sanitise_text(response["context"])
    return response


def sanitise_expansion_response(response):
    for key in ("suggestions", "missing_provenance"):
        for suggestion in response.get(key) or []:
            if suggestion.get("reason") is not None:
                suggestion["reason"] = sanitise_text(suggestion["reason"])
    return response


def sanitise_literature_response(response):
    for finding in response.get("findings") or []:
        if finding.get("relevance") is not None:
            finding["relevance"] = sanitise_text(finding["relevance"])
    return response
