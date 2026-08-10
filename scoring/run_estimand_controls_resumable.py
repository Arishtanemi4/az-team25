import json
import os
import pickle
import time

import pandas as pd

import estimand
import sensitivity

DATA_DIR = estimand.DATA_DIR
RESOURCES_DIR = estimand.RESOURCES_DIR
CHECKPOINT_PATH = f"{RESOURCES_DIR}/.estimand_checkpoint.json"  # gitignored, ephemeral
LAYER_FRAMES_CACHE_PATH = f"{RESOURCES_DIR}/.estimand_layer_frames_cache.pkl"  # gitignored, ephemeral
FINAL_OUT_PATH = f"{RESOURCES_DIR}/estimand_results.json"


def load_checkpoint():
    if os.path.exists(CHECKPOINT_PATH):
        with open(CHECKPOINT_PATH) as f:
            return json.load(f)
    return {}


def save_checkpoint(ckpt):
    tmp = CHECKPOINT_PATH + ".tmp"
    with open(tmp, "w") as f:
        json.dump(ckpt, f)
    os.replace(tmp, CHECKPOINT_PATH)


def load_or_build_layer_frames(battery_gene_ids, data_dir):
    if os.path.exists(LAYER_FRAMES_CACHE_PATH):
        with open(LAYER_FRAMES_CACHE_PATH, "rb") as f:
            return pickle.load(f)
    layer_frames = sensitivity.load_layer_frames(battery_gene_ids, data_dir)
    tmp = LAYER_FRAMES_CACHE_PATH + ".tmp"
    with open(tmp, "wb") as f:
        pickle.dump(layer_frames, f)
    os.replace(tmp, LAYER_FRAMES_CACHE_PATH)
    return layer_frames


def main():
    t0 = time.time()
    gene_reference = pd.read_csv(f"{DATA_DIR}/gene_reference.csv")
    cell_lines = pd.read_csv(f"{DATA_DIR}/cell_lines.csv")
    coverage = pd.read_csv(f"{DATA_DIR}/coverage.csv")
    battery = sensitivity.load_battery(f"{RESOURCES_DIR}/sensitivity_battery.json")
    print(f"[{time.time() - t0:.1f}s] reference tables + battery loaded", flush=True)

    ckpt = load_checkpoint()
    if ckpt:
        print(f"[{time.time() - t0:.1f}s] resuming from checkpoint, already have: "
              f"{sorted(ckpt.keys())}", flush=True)

    battery_gene_ids = sensitivity.resolve_battery_genes(battery, gene_reference)
    layer_frames = load_or_build_layer_frames(battery_gene_ids, DATA_DIR)
    print(f"[{time.time() - t0:.1f}s] real-battery Tier-1 layer_frames ready", flush=True)

    if "known_negative_pairs" not in ckpt:
        query = next(q for q in battery if q["name"] == estimand.KNOWN_NEGATIVE_PAIR_QUERY_NAME)
        ckpt["known_negative_pairs"] = estimand.run_known_negative_pair_control(
            query, gene_reference, cell_lines, coverage, layer_frames=layer_frames
        )
        save_checkpoint(ckpt)
        print(f"[{time.time() - t0:.1f}s] known_negative_pairs done", flush=True)

    if "known_negative_pairs_lung" not in ckpt:
        ckpt["known_negative_pairs_lung"] = estimand.run_known_negative_pair_control(
            estimand.KNOWN_NEGATIVE_PAIR_QUERY_LUNG, gene_reference, cell_lines, coverage,
            layer_frames=layer_frames,
        )
        save_checkpoint(ckpt)
        print(f"[{time.time() - t0:.1f}s] known_negative_pairs_lung done", flush=True)

    if "real_baseline" not in ckpt:
        real_baseline = estimand._run_real_battery_baseline(
            battery, gene_reference, cell_lines, coverage, layer_frames=layer_frames
        )

        ckpt["real_baseline"] = [{"confidence_tier": r["confidence_tier"], "D": r["D"]} for r in real_baseline]
        save_checkpoint(ckpt)
        print(f"[{time.time() - t0:.1f}s] real_baseline done ({len(real_baseline)} lines)", flush=True)

    if "random_gene_lists" not in ckpt:
        shadow_battery = estimand.build_shadow_battery(battery, gene_reference, exclude=battery_gene_ids)
        ckpt["random_gene_lists"] = estimand.run_random_control(
            shadow_battery, ckpt["real_baseline"], gene_reference, cell_lines, coverage
        )
        save_checkpoint(ckpt)
        print(f"[{time.time() - t0:.1f}s] random_gene_lists done", flush=True)

    if "shuffled_lineage_null" not in ckpt:
        ckpt["shuffled_lineage_null"] = estimand.run_shuffled_lineage_null(
            battery, gene_reference, cell_lines, coverage, layer_frames=layer_frames
        )
        save_checkpoint(ckpt)
        print(f"[{time.time() - t0:.1f}s] shuffled_lineage_null done", flush=True)

    if "leave_one_lineage_out" not in ckpt:
        holdout_lineages = estimand.select_holdout_lineages(cell_lines)
        ckpt["leave_one_lineage_out"] = estimand.run_loo_control(
            battery, holdout_lineages, gene_reference, cell_lines, coverage, layer_frames=layer_frames
        )
        save_checkpoint(ckpt)
        print(f"[{time.time() - t0:.1f}s] leave_one_lineage_out done", flush=True)

    results = {
        "battery_size": len(battery),
        "known_negative_pairs": ckpt["known_negative_pairs"],
        "known_negative_pairs_lung": ckpt["known_negative_pairs_lung"],
        "random_gene_lists": ckpt["random_gene_lists"],
        "shuffled_lineage_null": ckpt["shuffled_lineage_null"],
        "leave_one_lineage_out": ckpt["leave_one_lineage_out"],
        "any_control_passed": any([
            ckpt["known_negative_pairs"]["pass"], ckpt["known_negative_pairs_lung"]["pass"],
            ckpt["random_gene_lists"]["pass"], ckpt["shuffled_lineage_null"]["pass"],
            ckpt["leave_one_lineage_out"]["pass"],
        ]),
    }
    with open(FINAL_OUT_PATH, "w") as f:
        json.dump(results, f, indent=2)
    os.remove(CHECKPOINT_PATH)
    if os.path.exists(LAYER_FRAMES_CACHE_PATH):
        os.remove(LAYER_FRAMES_CACHE_PATH)
    print(f"[{time.time() - t0:.1f}s] wrote {FINAL_OUT_PATH}", flush=True)


if __name__ == "__main__":
    main()
