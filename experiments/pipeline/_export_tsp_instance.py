#!/usr/bin/env python3
#

import argparse
import collections
import os
import re
import sys

from pathlib import Path
from xopen import xopen
from pprint import pprint

import numpy as np
import pandas as pd


def write_tsp_instance(arr, keys, fo, without_return=True):
    if without_return:
        # add a dummy city in zero distance from everything
        keys = ["_DUMMY_CITY_SEPARATOR_"] + keys
        arr = np.insert(arr, 0, 0, axis=0)
        arr = np.insert(arr, 0, 0, axis=1)

    dim = arr.shape[0]

    header = f"""NAME : {dim} points
TYPE : TSP
DIMENSION : {dim}
COMMENT : keys:{",".join(keys)}
EDGE_WEIGHT_TYPE : EXPLICIT
EDGE_WEIGHT_FORMAT : FULL_MATRIX
NODE_COORD_TYPE : NO_COORDS
DISPLAY_DATA_TYPE : NO_DISPLAY
EDGE_WEIGHT_SECTION
"""

    footer = """EOF
"""

    df = pd.DataFrame(arr)
    matrix = df.to_csv(sep=" ", index=False, header=False)

    fo.write(header)
    fo.write(matrix)
    fo.write(footer)


def create_matrix_from_file(fn, sel, worst):
    d = {}
    keys = set()
    with open(sel, "r") as f:
        selected_names = set(
            os.path.splitext(line.strip().split("/")[-1])[0] for line in f
        )
    # print(selected_names)
    with xopen(fn) as fo:
        for x in fo:
            x = x.strip()
            a, b, dist = x.split("\t")
            a = os.path.splitext(a)[0]
            b = os.path.splitext(b)[0]
            #           print("just read", a, b)
            if a not in selected_names or b not in selected_names:
                continue
            keys.update([a, b])
            dist = int(dist)
            d[(a, b)] = int(dist)
            d[(b, a)] = int(dist)
    keys = sorted(list(keys))
    n_real = len(keys)

    extra = 35 if n_real < 35 else 0
    keys += [keys[0]] * extra
    dim = n_real + extra

    if worst:
        orig_max = max(d.values())
        for k in d:
            d[k] = orig_max - d[k]

    max_dist = max(d.values())
    scale_factor = 1
    threshold = 10000

    if max_dist > threshold:  # concorde threshold
        print("scaling applied due to large distances")
        scale_factor = max_dist / threshold  # integer factor
        for k in d:
            d[k] = max(1, int(round(d[k] / scale_factor)))

    arr = np.zeros((dim, dim), dtype=np.int32)
    for i, x in enumerate(keys):
        for j, y in enumerate(keys):
            if x == y:
                arr[i, j] = 0
            else:
                arr[i, j] = d[(x, y)]

    while len(keys) < arr.shape[0]:
        keys.append(f"_DUMMY_FILLER_{len(keys)+1}_")

    return keys, arr

def process(fn, out_optimal, out_worst, sel):
    keys, arr = create_matrix_from_file(fn, sel, False)
    with open(f"{out_optimal}", "w") as fo:
        write_tsp_instance(arr, keys, fo)

    keys, worst_arr = create_matrix_from_file(fn, sel, True)
    with open(f"{out_worst}", "w") as fo:
        write_tsp_instance(worst_arr, keys, fo)

def main():
    parser = argparse.ArgumentParser(description="generate optimal and worst case TSP instances")
    parser.add_argument("fn", metavar="input.txt", help="distance file")
    parser.add_argument("sel", help="selection file")
    parser.add_argument(
        "-o",
        "--out-optimal",
        required=True,
        help="filename.tsp for optimal",
    )
    parser.add_argument(
        "-w",
        "--out-worst",
        required=True,
        help="filename.tsp for worst",
    )

    args = parser.parse_args()
    process(args.fn, args.out_optimal, args.out_worst, args.sel)


if __name__ == "__main__":
    main()
