from pathlib import Path

import verifier

REPO_ROOT = Path(__file__).resolve().parent.parent
SOURCE_FILES = [
    REPO_ROOT / "docs" / "plan" / "PARAMETERS.md",
    REPO_ROOT / "docs" / "reference" / "ALGORITHM_SPEC.md",
    REPO_ROOT / "docs" / "plan" / "METHOD_DECISION.md",
]


def flatten_numbers_from_files(source_files=SOURCE_FILES):
    values = []
    per_file = {}
    for path in source_files:
        if not path.exists():
            continue
        # PARAMETERS.md writes some negative values with a real Unicode minus sign (U+2212,
        # e.g. row 6's protein floor), not the ASCII hyphen verifier.extract_numeric_tokens's
        # regex matches -- normalize before extraction so the sign survives.
        text = path.read_text(encoding="utf-8").replace("−", "-")
        file_values = []
        for token in verifier.extract_numeric_tokens(text):
            cleaned = token.rstrip("%").replace(",", "")
            try:
                value = float(cleaned)
            except ValueError:
                continue
            if token.endswith("%"):
                value = value / 100.0
            file_values.append(value)
        values.extend(file_values)
        per_file[str(path.relative_to(REPO_ROOT)).replace("\\", "/")] = file_values

    return {"values": values, "per_file": per_file}


# evaluate.py has called this function under both names across edits -- keep both live rather
# than chase a single canonical name back and forth.
flatten_all_parameter_values = flatten_numbers_from_files
