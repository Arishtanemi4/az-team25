"""Makes the rag/ package importable. rag/narrator.py and its sibling modules use bare
imports (`import citations`, `import provider`, `import verifier`, ...), not a proper installed
package, so they only work once rag/ itself is on sys.path -- this mirrors the same
sys.path.insert bootstrap backend/scoring_service/app/lib/scoring_path.py already uses for
scoring/, so rag_service imports rag/ exactly the way scoring_service imports scoring/.
"""

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[4]
_RAG_DIR = str(_REPO_ROOT / "rag")

if _RAG_DIR not in sys.path:
    sys.path.insert(0, _RAG_DIR)
