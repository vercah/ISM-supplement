#!/usr/bin/env python3
import sys
import lzma
import argparse

from _sample_name import canonical_sample_name


def open_text(path):
    return lzma.open(path, "rt") if path.endswith(".xz") else open(path, "rt")


def main():
    parser = argparse.ArgumentParser(
        description="Convert pairwise distances to PHYLIP distance matrix.")
    parser.add_argument("dists", help="Pairwise distance file, plain or .xz")
    parser.add_argument("selection",
                        help="Selection file determining genome order")
    parser.add_argument("-o", "--output", default="-")
    args = parser.parse_args()

    names = []
    with open(args.selection) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                names.append(canonical_sample_name(line))

    name_to_idx = {name: i for i, name in enumerate(names)}
    n = len(names)
    matrix = [[0.0] * n for _ in range(n)]

    with open_text(args.dists) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue

            a, b, d = parts
            a = canonical_sample_name(a)
            b = canonical_sample_name(b)

            if a not in name_to_idx or b not in name_to_idx:
                continue

            i = name_to_idx[a]
            j = name_to_idx[b]
            matrix[i][j] = float(d)
            matrix[j][i] = float(d)

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    out.write(f"    {n}\n")
    for i, name in enumerate(names):
        row = "  ".join(f"{matrix[i][j]:.10f}" for j in range(n))
        out.write(f"{name}  {row}\n")

    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
