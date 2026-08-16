"""Makes the scoring/ package importable. scoring/score.py and its sibling modules use bare
imports (`import combine`, `import correlation`, ...), not a proper installed package, so they
only work once scoring/ itself is on sys.path -- this mirrors the same sys.path.insert bootstrap
scoring/tests/conftest.py already uses, so scoring_service imports scoring/ exactly the way its
own test suite does, without modifying a single line of scoring/ itself.
"""

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[4]
_SCORING_DIR = str(_REPO_ROOT / "scoring")

if _SCORING_DIR not in sys.path:
    sys.path.insert(0, _SCORING_DIR)
