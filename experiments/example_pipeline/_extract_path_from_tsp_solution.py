#! /mnt/lustre/helios-home/hendrver/miniconda3/envs/stromecky/bin/python

import argparse
import collections
import os
import re
import sys

from pathlib import Path
from xopen import xopen
from pprint import pprint


def process(inst_fn, sol_fn):
    with open(inst_fn) as fo:
        for x in fo:
            if x.find("keys") != -1:
                keys = x.strip().partition("keys:")[2].split(",")
                #print(keys)
    order = []
    with open(sol_fn) as fo:
        for i, x in enumerate(fo):
            if i == 0:
                dim = int(x)
                continue
            x = x.strip()
            order.extend(map(int, x.split()))
    #print(order)
    assert dim == len(order), (dim, len(order))
    assert dim == len(keys), (dim, len(keys))
    #if len(keys) < dim:
    #    padding = [f"_DUMMY_FILLER_{i}_" for i in range(dim - len(keys))]
    #    keys.extend(padding)
    #assert len(keys) == dim == len(order), (
    #        f"dimension/order/keys mismatch: dim={dim},"
    #        f"order={len(order)}, keys={len(keys)}")
    sorted_keys = [keys[i] for i in order]
    assert sorted_keys[0] == "_DUMMY_CITY_SEPARATOR_", sorted_keys
    
    real_tour = [
        name for name in sorted_keys if not name.startswith("_DUMMY_")
    ]
    
    print(*real_tour, sep="\n")


def main():

    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        'inst_fn',
        metavar='instance.tsp',
        help='',
    )

    parser.add_argument(
        'sol_fn',
        metavar='instance.sol',
        help='',
    )

    args = parser.parse_args()

    process(args.inst_fn, args.sol_fn)


if __name__ == "__main__":
    main()
