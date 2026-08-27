#!/bin/bash
# Runs once, before backend/frontend start (docker-compose.yml's `pipeline` service), turning a
# fresh clone's data/raw/ + data/augmented/ (the only two data/ paths git-lfs tracks) into
# everything scoring_service and rag_service need at request time. Idempotent -- each step is
# skipped if its output already exists, so a second `docker compose up` is fast.
#
# Steps 1-2 are required: scoring_service can't produce correct D-values without them, so a
# failure here fails the whole pipeline (set -e) and backend never starts (depends_on:
# condition: service_completed_successfully in docker-compose.yml).
#
# Steps 3-4 are optional: rag_service already degrades gracefully without their output
# (rag_app/services/rag_service.py::_search_c1 returns no context; graph_service.py reports
# available: false) per this project's own "no blank-canvas-by-crash" rule, so a failure there is
# logged and skipped rather than blocking the deploy.
set -euo pipefail
cd /app

echo "=== [1/4] preprocessing: data/raw + data/augmented -> data/processed ==="
if [ -f data/processed/cell_lines.csv ]; then
  echo "data/processed/cell_lines.csv already exists -- skipping preprocessing/preprocess.py"
else
  python preprocessing/preprocess.py
fi

echo "=== [2/4] scoring: desirability calibration constants ==="
if [ -f scoring/resources/desirability_constants.json ]; then
  echo "scoring/resources/desirability_constants.json already exists -- skipping scoring/build_desirability_constants.py"
else
  python scoring/build_desirability_constants.py
fi

echo "=== [3/4] rag: C1 documentation index (optional) ==="
if [ -f data/processed/rag/c1_index/chunks.jsonl ]; then
  echo "C1 index already exists -- skipping rag/build_c1_index.py"
elif python rag/build_c1_index.py; then
  echo "C1 index built."
else
  echo "WARNING: rag/build_c1_index.py failed -- continuing without it. Methodology Q&A and" \
       "narration will just retrieve no documentation context, not fail."
fi

echo "=== [4/4] rag: knowledge graph (optional) ==="
if [ -f data/augmented/rag/knowledge_graph/derived/knowledge_graph_edges.parquet ]; then
  echo "Knowledge graph already built -- skipping."
elif [ -d data/augmented/rag/reactome/raw ] && [ -d data/augmented/rag/string_v12/raw ] \
     && [ -f data/augmented/rag/biogrid/raw/BIOGRID-ORGANISM-5.0.259.tab3.zip ]; then
  if python rag/build_knowledge_graph.py; then
    echo "Knowledge graph built."
  else
    echo "WARNING: rag/build_knowledge_graph.py failed -- continuing without it. The" \
         "/graph/neighborhood endpoint will report available: false."
  fi
else
  echo "data/augmented/rag/ raw inputs (Reactome/STRING/BioGRID) not present on this machine --" \
       "skipping. The /graph/neighborhood endpoint will report available: false until this is built."
fi

echo "=== Pipeline complete ==="
