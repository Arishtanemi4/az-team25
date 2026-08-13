import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# One row per (author, year) fragment inside a PAPERS.md "Paper" cell, e.g.
# "**KG-RAG / SPOKE** (Soman et al. 2024)" -> "soman 2024". Handles both "et al." and a lone
# second author ("X & Y Year") loosely -- it only needs the first surname and the year, since
# that is what a generator would naturally write when citing a source it was handed.
_AUTHOR_YEAR_RE = re.compile(r"\(([A-Z][A-Za-z\-]+)(?:[^()]*?)\s(\d{4})\)")

# A markdown table row: leading "|", then cells separated by "|". Skips the "|---|---|" separator
# row (a cell of only dashes/colons/spaces).
_TABLE_ROW_RE = re.compile(r"^\|(.+)\|\s*$")


def _split_row(line):
    cells = [c.strip() for c in line.strip().strip("|").split("|")]
    return cells


def _is_separator_row(cells):
    return all(re.fullmatch(r":?-+:?", c) for c in cells if c)


def _normalize_key(key):
    return re.sub(r"\s+", " ", key.strip().lower())


def parse_paper_citations(papers_md_path):
    entries = {}
    text = papers_md_path.read_text(encoding="utf-8")
    for line in text.splitlines():
        match = _TABLE_ROW_RE.match(line)
        if not match:
            continue
        cells = _split_row(line)
        if _is_separator_row(cells) or not cells:
            continue
        row_text = " | ".join(cells)
        for author, year in _AUTHOR_YEAR_RE.findall(row_text):
            key = _normalize_key(f"{author} {year}")
            entries[key] = {
                "kind": "paper",
                "key": key,
                "source_file": str(papers_md_path.relative_to(REPO_ROOT)),
                "row_text": row_text,
            }
    return entries


def parse_data_source_citations(data_sources_md_path):
    entries = {}
    text = data_sources_md_path.read_text(encoding="utf-8")
    lines = text.splitlines()
    in_location_table = False
    prev_cells = None
    for line in lines:
        match = _TABLE_ROW_RE.match(line)
        if not match:
            in_location_table = False
            prev_cells = None
            continue
        cells = _split_row(line)
        if _is_separator_row(cells):
            # The row just before a separator is this table's header.
            in_location_table = bool(prev_cells) and prev_cells[0].strip().lower() == "location"
            prev_cells = cells
            continue
        if in_location_table and len(cells) >= 2:
            location = cells[0].strip("`")
            if location:
                key = _normalize_key(location)
                entries[key] = {
                    "kind": "data_source",
                    "key": key,
                    "source_file": str(data_sources_md_path.relative_to(REPO_ROOT)),
                    "row_text": " | ".join(cells),
                }
        prev_cells = cells
    return entries


def build_registry(repo_root=REPO_ROOT):
    registry = {}
    for papers_md in sorted((repo_root / "docs" / "research").glob("*/PAPERS.md")):
        registry.update(parse_paper_citations(papers_md))
    data_sources_md = repo_root / "docs" / "reference" / "DATA_SOURCES.md"
    if data_sources_md.exists():
        registry.update(parse_data_source_citations(data_sources_md))
    return registry


def resolve(citation_text, registry):
    return registry.get(_normalize_key(citation_text))
