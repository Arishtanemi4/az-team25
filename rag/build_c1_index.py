from pathlib import Path

import chunking
import retrieval

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "data" / "processed" / "rag" / "c1_index"


def main():
    print("Chunking the C1 corpus (docs/**.md, module docstrings, build_manifest.json)...")
    chunks = chunking.build_c1_chunks()
    print(f"  {len(chunks)} chunks from {len({c['file'] for c in chunks})} files")

    print(f"Encoding with {retrieval.C1_DENSE_MODEL} and building the FAISS index...")
    index = retrieval.HybridIndex(chunks, dense_model_name=retrieval.C1_DENSE_MODEL)
    index.build_dense_index()

    print(f"Saving to {OUT_DIR} ...")
    index.save(OUT_DIR)
    print("Done.")


if __name__ == "__main__":
    main()
