#!/mnt/lustre/helios-home/hendrver/miniconda3/envs/stromecky/bin/python

import numpy as np
import sys
import itertools
import argparse
import progressbar
from xopen import xopen
from joblib import Parallel, delayed

n_cores = 4

def _convert_line(i, x):
    name, binstring = x.strip().split("\t")
    binstring2 = (np.frombuffer(binstring.encode('utf-8'), dtype='int8') - 48).astype("bool")
    return i, name, binstring2


def load_binstrings(fn, n, l, r):
    tr = []
    total_lines = sum(1 for _ in xopen(fn))
    bar = progressbar.ProgressBar(max_value=total_lines)
    with xopen(fn) as fo:
        tr = Parallel(n_jobs=n_cores)(
            delayed(_convert_line)(i, x)                                    for i, x in enumerate(bar(fo))                                  if (i % n == l or i % n == r))
 #       new_idx = 0
 #           if line.strip().split("\t")[0] in selected_names:
 #               if new_idx % n == l or new_idx % n == r:
 #                  tr.append(_convert_line(new_idx, line))
 #                new_idx += 1
    tr_l = filter(lambda x: x[0] % n == l, tr)
    tr_r = filter(lambda x: x[0] % n == r, tr)
    return list(tr_l), list(tr_r)

def _compare_pair(x, y):
    _, x_n, x_bin = x
    _, y_n, y_bin = y
    dist = np.sum(x_bin ^ y_bin)
    return x_n, y_n, dist

def calculate_all_hamdist(l_bs, r_bs):
    res = Parallel(n_jobs=n_cores)(
        delayed(_compare_pair)(x, y)
        for (x, y) in progressbar.progressbar(itertools.product(l_bs, r_bs)))
    for x in res:
        print(*x, sep="\t")

def process(fn, n, l, r):
    assert n > l >= 0, (n, l)
    assert n > r >= 0, (n, r)
    print("Loading binstrings", file=sys.stderr)
    l_bs, r_bs = load_binstrings(fn, n, l, r)
    print("Calculating all distances", file=sys.stderr)
    calculate_all_hamdist(l_bs, r_bs)

def main():

    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        '-n',
        metavar='int',
        dest='n',
        type=int,
        default='',
        required=True,
        help='number of groups',
    )

    parser.add_argument(
        '-i',
        metavar='int',
        dest='l',
        type=int,
        default='',
        required=True,
        help='left group {0, ..., #number-1}',
    )

    parser.add_argument(
        '-j',
        metavar='int',
        dest='r',
        type=int,
        default='',
        required=True,
        help='right group {0, ..., #number-1}',
    )

    parser.add_argument(
        'fn',
        metavar='binstrings.txt[.xz]',
        help='filename',
    )

    #parser.add_argument(
    #    'sel_fn',
    #    help='selection'
    #)
    
    args = parser.parse_args()
    process(args.fn, args.n, args.l, args.r)

if __name__ == "__main__":
    main()

