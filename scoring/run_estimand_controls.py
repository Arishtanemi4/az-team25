import time

import pandas as pd

import estimand

DATA_DIR = estimand.DATA_DIR


def main():
    t0 = time.time()
    gene_reference = pd.read_csv(f"{DATA_DIR}/gene_reference.csv")
    cell_lines = pd.read_csv(f"{DATA_DIR}/cell_lines.csv")
    coverage = pd.read_csv(f"{DATA_DIR}/coverage.csv")
    print(f"[{time.time() - t0:.1f}s] reference tables loaded", flush=True)

    results = estimand.run_all_controls(gene_reference, cell_lines, coverage)
    print(f"[{time.time() - t0:.1f}s] all four controls complete", flush=True)

    print("\n--- Summary ---")
    print(f"known_negative_pairs: pass={results['known_negative_pairs']['pass']}, "
          f"confident_fraction={results['known_negative_pairs']['confident_fraction']:.4f}")
    print(f"random_gene_lists: pass={results['random_gene_lists']['pass']}, "
          f"collapse_fraction={results['random_gene_lists']['low_insufficient_collapse_fraction']:.4f}")
    print(f"shuffled_lineage_null: pass={results['shuffled_lineage_null']['pass']}, "
          f"n_passed={results['shuffled_lineage_null']['n_passed']}/"
          f"{results['shuffled_lineage_null']['n_informative']}")
    print(f"leave_one_lineage_out: pass={results['leave_one_lineage_out']['pass']}, "
          f"median_membership_change={results['leave_one_lineage_out']['median_membership_change']}")
    print(f"\nany_control_passed: {results['any_control_passed']}")
    print(f"\nWrote scoring/resources/estimand_results.json ({time.time() - t0:.1f}s total)")


if __name__ == "__main__":
    main()
