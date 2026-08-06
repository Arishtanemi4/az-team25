
import time

import build_desirability_constants
import build_essentiality_constants
import build_extended_desirability_constants

BUILDERS = [
    ("RNA desirability constants (L/T)", build_desirability_constants),
    ("Extended desirability constants (protein/dependency/copy-number L/T)", build_extended_desirability_constants),
    ("Essentiality constants (pan-essential gene flags)", build_essentiality_constants),
]


def main() -> None:
    overall_t0 = time.time()
    for label, module in BUILDERS:
        print(f"\n=== {label} ===")
        t0 = time.time()
        module.main()
        print(f"({time.time() - t0:.1f}s)")
    print(f"\nAll calibration constants built in {time.time() - overall_t0:.1f}s")


if __name__ == "__main__":
    main()
