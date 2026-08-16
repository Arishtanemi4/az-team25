"""Paths and settings for rag_service. rag/'s own modules (narrator.py, methodology_agent.py,
etc.) already resolve their own data paths (the C1 index, the knowledge graph parquet files,
docs/**.md) as absolutes relative to rag/'s own location -- this module does not re-declare any
of those, only what the API layer itself needs: CORS origins, same convention as
gene_service/scoring_service's own app/config.py.
"""

import os

CORS_ORIGINS = os.environ.get(
    "RAG_SERVICE_CORS_ORIGINS", "http://localhost:5173,http://127.0.0.1:5173"
).split(",")
