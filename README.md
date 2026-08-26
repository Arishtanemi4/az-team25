# CellLineSelector

A tool built for AstraZeneca that takes a researcher's genes of interest and hands back a short,
ranked list of cancer cell lines to actually work with — along with the evidence behind every
pick, so nobody has to take the ranking on faith.

Give it a set of genes you want included or excluded, and it walks through expression, protein,
mutation, copy-number, dependency, and fusion evidence across ~2,100 cell lines, scores each
candidate, and explains its reasoning in plain language. An AI assistant on top can narrate the
result set, answer questions about the methodology, and pull supporting literature.

## Quick start

The whole thing — preprocessing, scoring, backend, frontend — comes up with one command:

```bash
cp .env.example .env        # add your NVIDIA_API_KEY (used for the AI assistant features)
docker compose up --build
```

First run takes a while: it's turning the raw data into the tables everything else depends on.
After that, `docker compose up` is fast, since it skips anything already built.

Once it's up:
- Frontend: [http://localhost:5173](http://localhost:5173)
- Backend API: [http://localhost:8000](http://localhost:8000) (`/health`, `/docs` for the OpenAPI schema)

> Preprocessing loads some large files into memory and can get killed by Docker Desktop's default
> memory limit on the very first run. If the `pipeline` service exits with `Killed` in its logs,
> bump Docker Desktop's memory allocation (Settings → Resources) to 12GB+ and try again — every
> run after the first is unaffected.

## Manual setup

Prefer working outside Docker (notebooks, debugging, faster iteration)? Set up a local environment
instead.

**Conda (recommended):**

```bash
conda env create -f environment.yaml
conda activate az
```

**Pip:**

```bash
conda create -n az python=3.14
conda activate az
pip install -r requirements.txt
```

**Jupyter kernel**, if you want to run the notebooks:

```bash
python -m ipykernel install --user --name=az --display-name="AstraZeneca"
```

From there, `python preprocessing/preprocess.py` builds `data/processed/` from `data/raw/` +
`data/augmented/`, and each backend service can be run standalone (see `backend/main.py` for the
combined app, or run any `backend/*_service` on its own for its own test suite).

## Project structure

A two-line tour of each part of the repo:

**`preprocessing/`** — Turns the raw data drops into the clean, joined tables everything else
reads from. Handles identity resolution, missingness, and unit scaling so nothing downstream has
to guess.

**`scoring/`** — The actual ranking algorithm: desirability scoring per gene/layer, correlation
weighting, and the calibration that keeps scores comparable across lineages.

**`rag/`** — The AI layer. Narrates a result set in plain English, answers methodology questions,
searches literature, and answers knowledge-graph queries — every claim grounded in real evidence,
nothing hallucinated.

**`backend/`** — Five FastAPI services (genes, scoring, RAG, data, research) that run together as
one process. Each one also works standalone for its own tests.

**`frontend/`** — The React + TypeScript UI researchers actually use: build a query, see the
ranked results, and open the AI assistant for narration or a literature search.

**`data/`** — Where everything reads from and writes to. Only `raw/` and `augmented/` are checked
in (via Git LFS); `processed/` and everything else is generated locally or by the Docker pipeline.

**`eda/`** — Exploratory notebooks used to understand the data before any of the above got built.
Nothing here is imported by the actual pipeline.

**`validation/`** — Cross-checks the pipeline's output against independent published datasets, so
the ranking isn't just internally consistent but externally sane.

**`evaluation/`** — A critical look back at the finished system: what it gets right, where it's
weaker, and what we'd do differently.

**`deploy/`** — The scripts that make `docker compose up` a one-command setup, running
preprocessing and calibration automatically before anything else starts.

## Data

`data/` is gitignored except `raw/` and `augmented/`, which ship via Git LFS — everything else
(`processed/`, the RAG indexes, calibration constants) is regenerable output, built locally by
`preprocessing/preprocess.py` and friends, or automatically by the Docker pipeline above.
