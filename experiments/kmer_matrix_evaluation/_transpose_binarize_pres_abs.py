#!/usr/bin/env python3
import numpy as np
import sys

if len(sys.argv) < 2:
    sys.exit("Usage: python script.py input.tsv")

filename = sys.argv[1]
with open(filename, 'r') as f:
    header = f.readline().strip().split('\t')
    print("opened file", file=sys.stderr)
    ncols = len(header)
    print("loaded ncols", file=sys.stderr)
    data = np.loadtxt(filename,
                      delimiter='\t',
                      skiprows=1,
                      usecols=range(1, ncols),
                      dtype=int)
    print("loaded data", file=sys.stderr)
    data[data != 0] = 1
    data_T = data.T
    print("binarized", file=sys.stderr)
    for i, row in enumerate(data_T):
        name = header[i+1].strip()
        print(name, file=sys.stderr)
        if name.endswith('.fa'):
            name = name[:-3]
        binary_str = ''.join(str(int(x)) for x in row)
        print(name + "\t" + binary_str)
