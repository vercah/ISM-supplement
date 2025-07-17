#! /mnt/lustre/helios-home/hendrver/miniconda3/envs/stromecky/bin/python

import argparse
import collections
import itertools
import multiprocessing
import sys
import os

import numpy as np

from xopen import xopen
import progressbar

from joblib import Parallel, delayed
# BitVector slower than Numpy
#from BitVector import BitVector

n_cores = multiprocessing.cpu_count()

def _convert_line(x):
    name, binstring = x.strip().split("\t")
    binstring2 = (np.frombuffer(binstring.encode('utf-8'), dtype='int8') - 48).astype("bool")
    return name, binstring2


def load_binstrings(fn):
    tr = []
    with xopen(fn, "rt") as fo:
        tr = Parallel(n_jobs=n_cores)(delayed(_convert_line)(x)
                                      for x in progressbar.progressbar(fo))

    d = {}
    for n, b in tr:
        d[n] = b
    return d


def get_number_of_runs(d, keys):
    # assert len(d) == len(keys), (len(d), len(keys))
    existing_keys = [k for k in keys if k in d]
    missing = [k for k in keys if k not in d]
    if missing:
        print("Warning: these keys are missing and will be skipped:", ", ".join(missing))
    if not existing_keys:
        raise ValueError("No valid keys found in the matrix")
    
    n_runs = 0
    n_changes = 0
    for i, x in enumerate(keys):
        if i == 0:
            n_runs = len(d[x])
        else:
            diff = np.sum(d[x] ^ d[prev_x])
            n_runs += diff
            n_changes += diff
        prev_x = x
    result = f"n_bit_changes\t{n_changes}\nn_runs\t{n_runs}\n"
    return result

def process(order_fn, bs_fn, outfiles):
    #print("Loading binstrings and keys", file=sys.stderr)
    if isinstance(order_fn, str):
        order_fn = [order_fn]
    d = load_binstrings(bs_fn)
    for i, order_file in enumerate(order_fn):
        with open(order_file) as fo:
            keys = [x.strip() for x in fo]
        result = get_number_of_runs(d, keys)
        if outfiles and i < len(outfiles):
            out_path = outfiles[i]
            out_dir = os.path.dirname(out_path)
            if out_dir and not os.path.exists(out_dir):
                os.makedirs(out_dir)
            with open(outfiles[i], 'w') as out:
                out.write(result)
        else:
            print("Results done for", order_file)
            #print(result)

def main():

    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        'bs_fn',
        metavar='bitstrings.txt[.xz]',
        help='',
    )

    parser.add_argument(
        'order_fn',
        nargs='+',
        metavar='order.txt',
        help='',
    )
    
    parser.add_argument(
        '--outfiles',
        nargs='*',
        metavar='output.txt',
        help='Optional output files',
    )

    args = parser.parse_args()

    process(args.order_fn, args.bs_fn, args.outfiles)


if __name__ == "__main__":
    main()
