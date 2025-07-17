#! /mnt/lustre/helios-home/hendrver/miniconda3/envs/stromecky/bin/python
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


def create_matrix_from_files(fns, sel):
    d = {}
    keys = set()
    with open(sel, "r") as f:
        selected_names = set(os.path.splitext(line.strip().split("/")[-1])[0] for line in f)
    #print(selected_names)
    for fn in fns:
        with xopen(fn) as fo:
            for x in fo:
                x = x.strip()
                a, b, dist = x.split("\t")
     #           print("just read", a, b)
                if (a not in selected_names or b not in selected_names):
      #              print("one of them not in selected")
                    continue
                keys.update([a, b])
                dist = int(dist)
                # double-check there're no conflicting values
                try:
                    assert d[a, b] == dist
                except KeyError:
                    pass
                try:
                    assert d[b, a] == dist
                except KeyError:
                    pass
                d[(a, b)] = int(dist)
                d[(b, a)] = int(dist)
    keys = sorted(list(keys))
    n_real = len(keys)

    if n_real < 35:
        extra = 35
    else:
        extra = 0
    keys += [keys[0]] * extra

    dim = n_real+extra
    arr = np.zeros((dim, dim), dtype=np.int32)
    for i, x in enumerate(keys):
        for j, y in enumerate(keys):
            arr[i, j] = d[(x, y)]
    #print(arr)
    while len(keys) < arr.shape[0]:
        keys.append(f"_DUMMY_FILLER_{len(keys)+1}_")
    return keys, arr

    #print(df.to_csv(sep="\t"))


def process(fns, fo, sel):
    keys, arr = create_matrix_from_files(fns, sel)
    #print(arr)
    #df = pd.DataFrame(arr)
    #print(df.to_csv(sep=" ", index=False, header=False))
    arr = arr.max() - arr #reverting the distances
    write_tsp_instance(arr, keys, fo)


def main():

    parser = argparse.ArgumentParser(description="")

    parser.add_argument('-o',
                        metavar='output.txt',
                        dest='fo',
                        type=argparse.FileType('w'),
                        default='-',
                        help='output file [-]')

    parser.add_argument(
        'fns',
        metavar='input.txt.xz',
        nargs="+",
        help='',
    )

    parser.add_argument(
        'sel',
    )

    args = parser.parse_args()
    process(args.fns, args.fo, args.sel)

if __name__ == "__main__":
    main()
