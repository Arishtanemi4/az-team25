import re

import citations

# Citation markers are pulled out FIRST so a bracketed year ("[Soman 2024]") is never also
# treated as a bare number needing numeric-fidelity verification -- a citation year is not
# evidence, it identifies a source.
_CITATION_MARKER_RE = re.compile(r"\[([^\[\]]+)\]")

# Digit runs with optional thousands-commas, a decimal point, and a trailing percent sign --
# deliberately permissive (it will also catch date fragments like the "2026"/"08"/"12" inside a
# plain "2026-08-12" mentioned in prose). That is a disclosed over-approximation, not a bug: a
# few date fragments failing to match the evidence record would show up as false failures
# (fail CLOSED, safe direction), never as a false pass. Bare four-digit numbers in the
# 1900-2100 "looks like a year" range are excluded outright, since no real evidence value in
# this project's scoring path (docs/plan/PARAMETERS.md) ever lands in that range.
#
# The lookaround pair matters: without it, this project's own identifiers -- ModelID
# "ACH-000444", Ensembl "ENSG00000146648" -- would be sliced into spurious numeric tokens
# ("-000444", "146648") that fail against the evidence record on every single result, since
# narrator/agent prose constantly repeats these IDs verbatim. Excluding a match glued to a
# letter/digit/underscore/hyphen on either side keeps those identifiers out entirely, at the
# cost of also skipping a genuine number glued directly to text with no separating space or
# punctuation (e.g. "0.5-fold") -- rare in this project's prose and a missed check, not a
# false pass, so it stays on the safe side of R1's fail-closed rule.
_IDENTIFIER_BOUNDARY = r"[A-Za-z0-9_\-]"
_NUMBER_RE = re.compile(
    rf"(?<!{_IDENTIFIER_BOUNDARY})-?\d[\d,]*(?:\.\d+)?%?(?!{_IDENTIFIER_BOUNDARY})"
)
_YEAR_RE = re.compile(r"^(19|20)\d{2}$")


def extract_citation_keys(text):
    return _CITATION_MARKER_RE.findall(text)


def extract_numbers_raw(text):
    tokens = []
    for match in _NUMBER_RE.finditer(text):
        token = match.group()
        bare = token.rstrip("%").replace(",", "")
        if _YEAR_RE.match(bare):
            continue
        tokens.append(token)
    return tokens


def extract_numeric_tokens(text):
    stripped = _CITATION_MARKER_RE.sub(" ", text)
    tokens = []
    for match in _NUMBER_RE.finditer(stripped):
        token = match.group()
        bare = token.rstrip("%").replace(",", "")
        if _YEAR_RE.match(bare):
            continue
        tokens.append(token)
    return tokens


def _flatten_numeric_leaves(obj, out):
    if isinstance(obj, dict):
        for value in obj.values():
            _flatten_numeric_leaves(value, out)
    elif isinstance(obj, (list, tuple)):
        for value in obj:
            _flatten_numeric_leaves(value, out)
    elif isinstance(obj, bool):
        return
    elif isinstance(obj, (int, float)):
        out.append(float(obj))


def verify_numeric_fidelity(text, evidence_record):
    leaves = []
    _flatten_numeric_leaves(evidence_record, leaves)

    failures = []
    for raw_token in extract_numeric_tokens(text):
        is_percent = raw_token.endswith("%")
        cleaned = raw_token.rstrip("%").replace(",", "")
        try:
            token_value = float(cleaned)
        except ValueError:
            failures.append(raw_token)
            continue
        if is_percent:
            token_value = token_value / 100.0

        decimals = len(cleaned.split(".")[1]) if "." in cleaned else 0
        matched = any(round(leaf, decimals) == round(token_value, decimals) for leaf in leaves)
        if not matched:
            failures.append(raw_token)
    return failures


def verify_citations(text, registry):
    return [key for key in extract_citation_keys(text) if not citations.resolve(key, registry)]


def verify(text, evidence_record, registry=None):
    if registry is None:
        registry = citations.build_registry()
    numeric_failures = verify_numeric_fidelity(text, evidence_record)
    citation_failures = verify_citations(text, registry)
    return {
        "passed": not numeric_failures and not citation_failures,
        "numeric_failures": numeric_failures,
        "citation_failures": citation_failures,
    }
