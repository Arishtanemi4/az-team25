import ast
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CHUNK_WORDS = 512
DEFAULT_OVERLAP_WORDS = 64

_HEADING_RE = re.compile(r"^(#{1,6})\s+(.*)$")


def _iter_heading_sections(markdown_text):
    stack = []  # list of (level, title)
    current_body_lines = []

    def current_path():
        return tuple(title for _, title in stack)

    for line in markdown_text.splitlines():
        match = _HEADING_RE.match(line)
        if match:
            if current_body_lines:
                yield current_path(), "\n".join(current_body_lines).strip()
                current_body_lines = []
            level = len(match.group(1))
            title = match.group(2).strip()
            while stack and stack[-1][0] >= level:
                stack.pop()
            stack.append((level, title))
        else:
            current_body_lines.append(line)

    if current_body_lines:
        yield current_path(), "\n".join(current_body_lines).strip()


def _split_into_word_chunks(text, chunk_words, overlap_words):
    paragraphs = [p.strip() for p in text.split("\n\n") if p.strip()]
    if not paragraphs:
        return []

    chunks = []
    current_words = []
    for paragraph in paragraphs:
        paragraph_words = paragraph.split()
        if current_words and len(current_words) + len(paragraph_words) > chunk_words:
            chunks.append(" ".join(current_words))
            current_words = current_words[-overlap_words:] if overlap_words else []
        current_words.extend(paragraph_words)
    if current_words:
        chunks.append(" ".join(current_words))
    return chunks


def chunk_markdown(markdown_text, source_file, chunk_words=DEFAULT_CHUNK_WORDS,
                    overlap_words=DEFAULT_OVERLAP_WORDS):
    chunks = []
    for heading_path, body in _iter_heading_sections(markdown_text):
        for piece in _split_into_word_chunks(body, chunk_words, overlap_words):
            chunks.append({
                "file": source_file,
                "heading_path": " > ".join(heading_path) if heading_path else source_file,
                "section": heading_path[-1] if heading_path else source_file,
                "text": piece,
            })
    return chunks


def chunk_markdown_file(path, repo_root=REPO_ROOT, chunk_words=DEFAULT_CHUNK_WORDS,
                         overlap_words=DEFAULT_OVERLAP_WORDS):
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    source_file = str(path.resolve().relative_to(repo_root)).replace("\\", "/")
    return chunk_markdown(text, source_file, chunk_words, overlap_words)


def chunk_module_docstring(py_path, repo_root=REPO_ROOT):
    path = Path(py_path)
    source_file = str(path.resolve().relative_to(repo_root)).replace("\\", "/")
    tree = ast.parse(path.read_text(encoding="utf-8"))
    docstring = ast.get_docstring(tree)
    if not docstring:
        return []
    return [{
        "file": source_file,
        "heading_path": f"{source_file} (module docstring)",
        "section": "module docstring",
        "text": docstring.strip(),
    }]


# The C1 corpus scope, RAG_ARCHITECTURE.md SS2: "docs/**.md ..., the ~63 files under
# docs/research/*/summary/, module docstrings from preprocessing/*.py and scoring/*.py, and
# data/processed/build_manifest.json." docs/archive/ is excluded -- it holds superseded
# documents the project's own docs/_.md says to treat as history, and indexing them would let a
# retrieved chunk cite a decision that no longer holds.
MARKDOWN_GLOBS = ["docs/**/*.md"]
MARKDOWN_EXCLUDE_PREFIX = "docs/archive/"
PYTHON_DOCSTRING_GLOBS = ["preprocessing/*.py", "scoring/*.py", "rag/*.py"]


def build_c1_chunks(repo_root=REPO_ROOT, chunk_words=DEFAULT_CHUNK_WORDS,
                     overlap_words=DEFAULT_OVERLAP_WORDS):
    repo_root = Path(repo_root)
    chunks = []

    md_paths = sorted(
        p for glob in MARKDOWN_GLOBS for p in repo_root.glob(glob)
        if not str(p.relative_to(repo_root)).replace("\\", "/").startswith(MARKDOWN_EXCLUDE_PREFIX)
    )
    for path in md_paths:
        chunks.extend(chunk_markdown_file(path, repo_root, chunk_words, overlap_words))

    py_paths = sorted(p for glob in PYTHON_DOCSTRING_GLOBS for p in repo_root.glob(glob))
    for path in py_paths:
        chunks.extend(chunk_module_docstring(path, repo_root))

    manifest_path = repo_root / "data" / "processed" / "build_manifest.json"
    if manifest_path.exists():
        chunks.append({
            "file": "data/processed/build_manifest.json",
            "heading_path": "data/processed/build_manifest.json",
            "section": "build manifest",
            "text": manifest_path.read_text(encoding="utf-8"),
        })

    return chunks
