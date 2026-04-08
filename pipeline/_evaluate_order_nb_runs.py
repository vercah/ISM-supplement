#!/usr/bin/env python3

import argparse
import multiprocessing
import sys
import os

import numpy as np

from xopen import xopen
import progressbar

from joblib import Parallel, delayed
from _sample_name import canonical_sample_name
# BitVector slower than Numpy
#from BitVector import BitVector

n_cores = multiprocessing.cpu_count()


def load_distances(fn):
    """Load pairwise distances into a dict, storing only one orientation (a<b)."""
    d = {}
    with xopen(fn) as f:
        for line in f:
            a, b, dist = line.strip().split()
            a = canonical_sample_name(a)
            b = canonical_sample_name(b)
            dist = int(dist)
            if a == b:
                continue
            # enforce ordering
            key = (a, b) if a < b else (b, a)
            d[key] = dist
    return d


def get_distance(d, a, b):
    """Lookup distance between a and b, order-independent."""
    if a == b:
        return 0
    key = (a, b) if a < b else (b, a)
    return d[key]


def load_order(fn):
    """Load genome order from file, stripping extensions if present."""
    order = []
    with xopen(fn) as f:
        for line in f:
            name = canonical_sample_name(line)
            order.append(name)
    return order


def compute_changes(d, order, num_rows):
    n_changes_path = 0
    for i in range(len(order) - 1):
        a, b = order[i], order[i + 1]
        n_changes_path += get_distance(d, a, b)

    n_runs = n_changes_path + num_rows
    n_bit_changes = n_changes_path
    n_cylinder_runs = n_changes_path + get_distance(d, order[0],
                                                    order[len(order) - 1])
    return n_bit_changes, n_runs, n_cylinder_runs


def main():
    parser = argparse.ArgumentParser(
        description="Compute runs from distance matrix and genome order")
    parser.add_argument("dist_fn", help="distance file (pairwise distances)")
    parser.add_argument("order_fn", help="genome order file")
    parser.add_argument("--num-rows",
                        type=int,
                        required=True,
                        help="Number of rows in the original matrix")

    parser.add_argument("-o", "--out", help="output file (default: stdout)")
    args = parser.parse_args()

    d = load_distances(args.dist_fn)
    order = load_order(args.order_fn)
    n_changes, n_runs, n_cylinder_runs = compute_changes(
        d, order, args.num_rows)

    result = f"n_bit_changes\t{n_changes}\nn_runs\t{n_runs}\nn_cylinder_runs\t{n_cylinder_runs}\n"

    if args.out:
        with open(args.out, "w") as fo:
            fo.write(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
