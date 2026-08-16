# Build context is the repo root (see docker-compose.yml) so this can COPY scoring/, rag/, docs/,
# validation/, and preprocessing/resources/, which live outside backend/ -- gene_service,
# scoring_service, rag_service, and data_service now run together as one process
# (backend/main.py) instead of four.
FROM python:3.12-slim

WORKDIR /app

COPY backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY backend ./backend
COPY scoring ./scoring
COPY rag ./rag
COPY docs ./docs
COPY validation ./validation
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
