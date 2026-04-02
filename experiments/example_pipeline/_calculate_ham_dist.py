#!/mnt/lustre/helios-home/hendrver/miniconda3/envs/stromecky/bin/python
from tqdm import tqdm
import numpy as np
import sys
import argparse
from xopen import xopen
from collections import defaultdict
from joblib import Parallel, delayed
import multiprocessing

n_cores = multiprocessing.cpu_count() 

def process_chunk(lines, num_genomes):
    local_dists = np.zeros((num_genomes, num_genomes), dtype=np.int64)
    for line in lines:
        bits = np.frombuffer(line.strip().encode('ascii'), dtype=np.uint8) - 48
        xormatr = bits[:, None] ^ bits[None, :]
        local_dists += np.triu(xormatr, k=1).astype(np.int64)
    return local_dists

def chunked_iterable(iterable, size):
    chunk = []
    for line in iterable:
        chunk.append(line)
        if len(chunk) >= size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def compute_streaming_distances_parallel(fn, n_cores, chunk_size=1000):
    with xopen(fn) as fo:
        header = fo.readline().strip()
        names = header.split("\t")
        num_genomes = len(names)

        merged = np.zeros((num_genomes, num_genomes), dtype=np.int64)
        #process chunks in parallel and merge them as they come back

        #chunk_iter = chunked_iterable(fo, chunk_size)
        #for local in Parallel(n_jobs=n_cores, return_as="generator_unordered")(
        #        delayed(process_chunk)(chunk, num_genomes)
        #        for chunk in tqdm(chunk_iter, desc="Processing chunks")):
        #    merged += local

        #for chunk in tqdm(chunked_iterable(fo, chunk_size), desc="Processing chunks"):
        #    locals_list = Parallel(n_jobs=n_cores)(
        #        delayed(process_chunk)([line], num_genomes) for line in chunk
        #    )
        #    # merge immediately and discard locals
        #    for local in locals_list:
        #        merged += local

        results = Parallel(n_jobs=n_cores, return_as="generator")(
            delayed(process_chunk)(chunk, num_genomes)
            for chunk in chunked_iterable(fo, chunk_size)
        )

        for local in tqdm(results, desc="Processing chunks"):
            merged += local
            del local

    tqdm.write("Matrices merged, started writing", file=sys.stderr)
    for i in range(num_genomes):
        for j in range(i + 1, num_genomes):
            print(f"{names[i]}\t{names[j]}\t{merged[i, j]}")

def main():
    parser = argparse.ArgumentParser(description="Compute all pairwise Hamming distances between columns")
    parser.add_argument("fn", help="matrix file (rows=k-mers, cols=genomes)")
    args = parser.parse_args()

#    mat, names = load_matrix(args.fn)

    #compute_all_pairs(mat, names)
    compute_streaming_distances_parallel(args.fn, n_cores, 10000)
    

if __name__ == "__main__":
    main()

