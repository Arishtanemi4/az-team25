# Build context is the repo root (see docker-compose.yml) so this can COPY scoring/, rag/, docs/,
# validation/, and preprocessing/resources/, which live outside backend/ -- gene_service,
# scoring_service, rag_service, and data_service now run together as one process
# (backend/main.py) instead of four.
FROM python:3.12-slim

WORKDIR /app

# Without this, Python fully buffers stdout when it isn't attached to a TTY (true of every
# `docker compose`/`docker logs` stream) -- preprocess.py's progress prints (and uvicorn's own
# startup logs) would sit invisible in an internal buffer until the process exits, making a long
# first-run pipeline step look hung rather than working.
ENV PYTHONUNBUFFERED=1

COPY backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY backend ./backend
COPY scoring ./scoring
COPY rag ./rag
COPY docs ./docs
COPY validation ./validation
COPY deploy ./deploy
# The pipeline (see deploy/pipeline_entrypoint.sh, docker-compose.yml's `pipeline` service) runs
# preprocess.py from this same image to turn data/raw/ + data/augmented/ into data/processed/ --
# needs preprocess.py's own sibling modules (bare `import clean` etc., same convention rag/ and
# each backend/*_service use), not just resources/.
COPY preprocessing/preprocess.py preprocessing/clean.py preprocessing/joins.py \
     preprocessing/impute.py preprocessing/normalize.py preprocessing/scale.py \
     ./preprocessing/
COPY preprocessing/resources ./preprocessing/resources

# data/processed (gene/scoring/data reference tables, including the optional eda_aggregates.json
# built by backend/data_service/build_eda_aggregates.py) and data/ (rag's PubMed cache, C1 index,
# knowledge graph) are mounted as volumes at runtime, not baked into the image -- gitignored,
# per-machine generated data (see docker-compose.yml). validation/*.csv is small and tracked, so
# it is copied into the image rather than mounted.
ENV GENE_SERVICE_DATA_DIR=/app/data/processed
ENV SCORING_SERVICE_DATA_DIR=/app/data/processed
ENV SCORING_SERVICE_RNA_CONSTANTS_PATH=/app/scoring/resources/desirability_constants.json
ENV SCORING_SERVICE_EXTENDED_CONSTANTS_PATH=/app/scoring/resources/desirability_constants_extended.json
ENV SCORING_SERVICE_ESSENTIALITY_CONSTANTS_PATH=/app/scoring/resources/common_essential_genes.json
ENV SCORING_SERVICE_GENE_ROLE_PATH=/app/preprocessing/resources/gene_role_reference.csv
ENV SCORING_SERVICE_UNRESOLVED_SYMBOLS_PATH=/app/preprocessing/resources/unresolved_gene_symbols.json
ENV DATA_SERVICE_DATA_DIR=/app/data/processed
ENV DATA_SERVICE_VALIDATION_DIR=/app/validation

EXPOSE 8000
CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0", "--port", "8000"]
