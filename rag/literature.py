import hashlib
import json
import os
from pathlib import Path

from Bio import Entrez

REPO_ROOT = Path(__file__).resolve().parent.parent
CACHE_DIR = REPO_ROOT / "data" / "processed" / "rag" / "literature_cache"

Entrez.email = os.environ.get("ENTREZ_EMAIL", "az-team25-rag@example.invalid")
_entrez_api_key = os.environ.get("ENTREZ_API_KEY")
if _entrez_api_key:
    Entrez.api_key = _entrez_api_key


def _cache_path(key):
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()
    return CACHE_DIR / f"{digest}.json"


def cache_lookup(key):
    path = _cache_path(key)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return None


def _cache_store(key, value):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    _cache_path(key).write_text(json.dumps(value), encoding="utf-8")


def pubmed_search(terms, retmax=10):
    cache_key = f"search:{terms}:{retmax}"
    cached = cache_lookup(cache_key)
    if cached is not None:
        return cached
    handle = Entrez.esearch(db="pubmed", term=terms, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    pmids = list(record.get("IdList", []))
    _cache_store(cache_key, pmids)
    return pmids


def pmc_fetch(pmid):
    cache_key = f"fetch:{pmid}"
    cached = cache_lookup(cache_key)
    if cached is not None:
        return cached or None
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
    text = handle.read()
    handle.close()
    if not text.strip():
        _cache_store(cache_key, None)
        return None
    result = {"pmid": pmid, "text": text.strip()}
    _cache_store(cache_key, result)
    return result


def rerank(query, passages, top_k=None):
    import retrieval

    return retrieval.rerank(query, passages, model_name=retrieval.C4_RERANKER_MODEL, top_k=top_k)
