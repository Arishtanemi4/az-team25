from pathlib import Path

import numpy as np
from rank_bm25 import BM25Okapi

REPO_ROOT = Path(__file__).resolve().parent.parent

C1_DENSE_MODEL = "BAAI/bge-small-en-v1.5"
C1_RERANKER_MODEL = "BAAI/bge-reranker-v2-m3"
C4_DENSE_QUERY_MODEL = "ncbi/MedCPT-Query-Encoder"
C4_DENSE_ARTICLE_MODEL = "ncbi/MedCPT-Article-Encoder"
C4_RERANKER_MODEL = "ncbi/MedCPT-Cross-Encoder"

_model_cache = {}


def _get_sentence_transformer(model_name):
    if model_name not in _model_cache:
        from sentence_transformers import SentenceTransformer
        _model_cache[model_name] = SentenceTransformer(model_name)
    return _model_cache[model_name]


def _tokenize(text):
    return [t for t in "".join(c if c.isalnum() else " " for c in text.lower()).split() if t]


class HybridIndex:

    def __init__(self, chunks, dense_model_name=C1_DENSE_MODEL):
        self.chunks = chunks
        self.dense_model_name = dense_model_name
        self._bm25 = BM25Okapi([_tokenize(c["text"]) for c in chunks])
        self._dense_vectors = None
        self._faiss_index = None

    def build_dense_index(self):
        import faiss

        encoder = _get_sentence_transformer(self.dense_model_name)
        vectors = encoder.encode(
            [c["text"] for c in self.chunks], normalize_embeddings=True, show_progress_bar=False
        )
        vectors = np.asarray(vectors, dtype="float32")
        index = faiss.IndexFlatIP(vectors.shape[1])
        index.add(vectors)
        self._dense_vectors = vectors
        self._faiss_index = index
        return self

    def _bm25_ranking(self, query, top_k):
        scores = self._bm25.get_scores(_tokenize(query))
        order = np.argsort(scores)[::-1][:top_k]
        return [int(i) for i in order if scores[i] > 0]

    def _dense_ranking(self, query, top_k):
        if self._faiss_index is None:
            return []
        encoder = _get_sentence_transformer(self.dense_model_name)
        query_vector = encoder.encode([query], normalize_embeddings=True)
        query_vector = np.asarray(query_vector, dtype="float32")
        _, indices = self._faiss_index.search(query_vector, min(top_k, len(self.chunks)))
        return [int(i) for i in indices[0] if i >= 0]

    def search(self, query, top_k=10, rrf_k=60):
        bm25_ranking = self._bm25_ranking(query, top_k * 4)
        dense_ranking = self._dense_ranking(query, top_k * 4)

        fused_scores = {}
        for ranking in (bm25_ranking, dense_ranking):
            for rank, chunk_index in enumerate(ranking):
                fused_scores[chunk_index] = fused_scores.get(chunk_index, 0.0) + 1.0 / (rrf_k + rank + 1)

        ranked_indices = sorted(fused_scores, key=fused_scores.get, reverse=True)[:top_k]
        return [self.chunks[i] for i in ranked_indices]

    def save(self, out_dir):
        import faiss
        import json

        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(out_dir / "chunks.jsonl", "w", encoding="utf-8") as f:
            for chunk in self.chunks:
                f.write(json.dumps(chunk) + "\n")
        if self._faiss_index is not None:
            faiss.write_index(self._faiss_index, str(out_dir / "dense.faiss"))
        with open(out_dir / "meta.json", "w", encoding="utf-8") as f:
            json.dump({"dense_model_name": self.dense_model_name, "n_chunks": len(self.chunks)}, f)

    @classmethod
    def load(cls, index_dir):
        import faiss
        import json

        index_dir = Path(index_dir)
        with open(index_dir / "meta.json", encoding="utf-8") as f:
            meta = json.load(f)
        chunks = []
        with open(index_dir / "chunks.jsonl", encoding="utf-8") as f:
            for line in f:
                chunks.append(json.loads(line))
        instance = cls(chunks, dense_model_name=meta["dense_model_name"])
        faiss_path = index_dir / "dense.faiss"
        if faiss_path.exists():
            instance._faiss_index = faiss.read_index(str(faiss_path))
        return instance


def rerank(query, candidates, model_name=C1_RERANKER_MODEL, top_k=None):
    from sentence_transformers import CrossEncoder

    if model_name not in _model_cache:
        _model_cache[model_name] = CrossEncoder(model_name)
    encoder = _model_cache[model_name]

    pairs = [(query, c["text"]) for c in candidates]
    scores = encoder.predict(pairs)
    order = np.argsort(scores)[::-1]
    reranked = [candidates[i] for i in order]
    return reranked[:top_k] if top_k else reranked
