#!/usr/bin/env python3

import argparse
import sys
from collections import Counter
from _sample_name import canonical_sample_name


def main():
    parser = argparse.ArgumentParser(
        description="Validate that order files contain exactly the expected genomes")
    parser.add_argument("--selection", required=True,
                        help="Selection file (ground truth genome list)")
    parser.add_argument("--out", required=True,
                        help="Output validation stamp file")
    parser.add_argument("order_files", nargs="+",
                        help="Order files to validate")
    args = parser.parse_args()

    with open(args.selection) as f:
        expected = sorted(canonical_sample_name(line) for line in f if line.strip())

    ok = True
    for fn in args.order_files:
        with open(fn) as f:
            names = [canonical_sample_name(line) for line in f if line.strip()]

        # Check for duplicates within the file
        counts = Counter(names)
        dupes = {k: v for k, v in counts.items() if v > 1}
        if dupes:
            print(f"ERROR: {fn} has duplicate genomes: {dupes}", file=sys.stderr)
            ok = False

        # Check set equality with expected
        actual_sorted = sorted(set(names))
        if actual_sorted != expected:
            missing = set(expected) - set(names)
            extra = set(names) - set(expected)
            if missing:
                print(f"ERROR: {fn} missing genomes: {missing}", file=sys.stderr)
            if extra:
                print(f"ERROR: {fn} extra genomes: {extra}", file=sys.stderr)
            ok = False

    if not ok:
        sys.exit(1)

    with open(args.out, "w") as f:
        f.write("OK\n")


if __name__ == "__main__":
    main()
