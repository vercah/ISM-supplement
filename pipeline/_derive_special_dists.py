#!/usr/bin/env python3

import argparse
import math
import sys

from xopen import xopen
from _sample_name import canonical_sample_name


VARIANTS = ("logdiff", "avglog")


def load_counts(fn):
    counts = {}
    with open(fn) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            name, value = line.split("\t")
            counts[canonical_sample_name(name)] = int(value)
    return counts


def derive(variant, H, N_a, N_b):
    s = N_a + N_b
    denom = s - H
    if denom < 0:
        raise ValueError(
            f"N_a + N_b - H = {denom} is negative (N_a={N_a}, N_b={N_b}, H={H}); "
            "this violates H = |A △ B| ≤ N_a + N_b — likely a unit mismatch "
            "(kmer distances must pair with kmer counts, unitig with unitig)."
        )
    # denom == 0 means the two genomes share zero features (empty intersection),
    # for which the log-Jaccard formula is -∞. Clamp to 1 so these pairs get
    # the maximum finite distance the formula can express, mirroring how Mash
    # caps non-overlapping sketches at distance 1.0.
    denom = max(denom, 1)
    if variant == "logdiff":
        return math.log(s) - math.log(denom)
    if variant == "avglog":
        return -math.log(denom)
    raise ValueError(f"unknown variant {variant!r}")


def main():
    parser = argparse.ArgumentParser(
        description="Derive log-Jaccard-style distances from hamming + per-genome counts.")
    parser.add_argument("ham_dists", help="pairwise hamming distance file (may be .xz)")
    parser.add_argument("counts", help="per-genome counts TSV (name<TAB>N)")
    parser.add_argument("variant", choices=VARIANTS)
    parser.add_argument("out_path", help="output pairwise TSV (plain, compress externally)")
    args = parser.parse_args()

    counts = load_counts(args.counts)

    with xopen(args.ham_dists) as fin, open(args.out_path, "w") as fout:
        for line in fin:
            a, b, h = line.strip().split()
            a = canonical_sample_name(a)
            b = canonical_sample_name(b)
            H = int(h)
            N_a = counts[a]
            N_b = counts[b]
            d = derive(args.variant, H, N_a, N_b)
            fout.write(f"{a}\t{b}\t{d:.10f}\n")


if __name__ == "__main__":
    main()
