#!/usr/bin/env python3

import os, sys, argparse
import numpy as np
from joblib import Parallel, delayed
import time, datetime

from _sample_name import canonical_sample_name


def load_names(filenames_file):
    # returns 0-based list: names[0] is color_id=0
    names = []
    with open(filenames_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                names.append(canonical_sample_name(line))
    return names


def load_color_sets(color_sets_file):
    color_sets = {}
    with open(color_sets_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(maxsplit=2)
            # first token looks like: color_set_id=0
            color_set_id = int(parts[0].split("=")[1])

            if len(parts) > 2:
                # parse the rest of the line (color IDs) directly into a NumPy array
                colors = np.fromstring(parts[2], dtype=int, sep=" ")
            else:
                colors = []
            color_sets[color_set_id] = colors
    return color_sets


def iter_unitigs_info(unitigs_file, k):
    with open(unitigs_file) as f:
        csid = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # > unitig_id=... color_set_id=...
                parts = line[1:].split()
                csid = int(parts[1].split("=", 1)[1])
            else:
                n = max(0, len(line) - k + 1)
                if n > 0:
                    yield csid, n


def process_batch(batch, num_colors):
    #pid = os.getpid()
    #start_time = datetime.datetime.now().strftime("%H:%M:%S")
    #print(f"{start_time} | Process {pid} started", file=sys.stderr, flush=True)

    #time.sleep(5)
    local_unitig = np.zeros((num_colors, num_colors), dtype=np.int64)
    local_kmer = np.zeros((num_colors, num_colors), dtype=np.int64)

    for ones, n in batch:
        row = np.zeros(num_colors, dtype=np.uint8)
        row[ones] = 1
        zeros = np.flatnonzero(row == 0)

        local_unitig[np.ix_(ones, zeros)] += 1
        local_unitig[np.ix_(zeros, ones)] += 1
        local_kmer[np.ix_(ones, zeros)] += n
        local_kmer[np.ix_(zeros, ones)] += n

    #end_time = datetime.datetime.now().strftime("%H:%M:%S")
    #print(f"{end_time} | Process {pid} finished", file=sys.stderr, flush=True)
    return local_unitig, local_kmer


def count_kmers(unitigs_file, k, metadata_fn):
    total_kmers = 0
    for csid, n in iter_unitigs_info(unitigs_file, k):
        total_kmers += n

    with open(metadata_fn, "r", newline="") as f:
        exists = any(line.startswith("num_kmers=") for line in f)
    if not exists:
        with open(metadata_fn, "a", newline="") as f:
            f.write(f"num_kmers={total_kmers}\n")
    return total_kmers


def subbatch_iterable(batch):
    subbatch_size = 10000
    for i in range(0, len(batch), subbatch_size):
        yield batch[i:i + subbatch_size]


def preprocess_color_sets(batch, color_sets):
    member_batch = []
    for csid, n in batch:
        ones = list(set(color_sets.get(csid, [])))
        member_batch.append((ones, n))
    return member_batch


def unique_raw_batches(color_sets, batch_size):
    # unique_raw counts each dumped color set exactly once, unlike the
    # unitig- and k-mer-weighted matrices above.
    batch = []
    for colors in color_sets.values():
        batch.append((colors, 1))
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def build_and_compute(color_sets, unitigs_file, k, filenames_file, out_path,
                      dataset, total_kmers, cores, batch_size):
    names = load_names(filenames_file)
    num_colors = len(names)
    kmer_dists = np.zeros((num_colors, num_colors), dtype=np.int64)
    unitig_dists = np.zeros((num_colors, num_colors), dtype=np.int64)
    unique_raw_dists = np.zeros((num_colors, num_colors), dtype=np.int64)

    batch = []
    #    pbar = tqdm(desc="Processing kmers", unit="rows", total=total_kmers)
    processed_n = 0
    for csid, n in iter_unitigs_info(unitigs_file, k):
        processed_n += n
        batch.append((csid, n))
        if len(batch) >= batch_size:
            #print(f"Started processing color sets of a batch", file=sys.stderr, flush=True)
            member_batch = preprocess_color_sets(batch, color_sets)
            #print(f"Starting parallel", file=sys.stderr, flush=True)
            for u_matrix, k_matrix in Parallel(
                    n_jobs=cores, backend="threading", return_as="generator")(
                        delayed(process_batch)(subbatch, num_colors)
                        for subbatch in subbatch_iterable(member_batch)):
                unitig_dists += u_matrix
                kmer_dists += k_matrix

            #local_matrices = [process_batch(member_batch, num_colors)]
            #for local_unitig, local_kmer in local_matrices:
            #    unitig_dists += local_unitig
            #    kmer_dists   += local_kmer
            #print(f"Merging local matrices", file=sys.stderr, flush=True)
            #matrices_u = [u for u, _ in local_matrices]
            #if matrices_u:  # only sum if non-empty
            #    unitig_dists += np.sum(np.stack(matrices_u), axis=0)
            #matrices_k = [k for _, k in local_matrices]
            #if matrices_k:
            #    kmer_dists += np.sum(np.stack(matrices_k), axis=0)
            member_batch = []
            batch = []
            #pbar.update(processed_n)
            processed_n = 0
            #print(f"Everything is reset, going for a new batch", file=sys.stderr, flush=True)

    if batch:
        member_batch = preprocess_color_sets(batch, color_sets)
        #local_matrices = [process_batch(member_batch, num_colors)]
        #print(f"Processed last color sets of a batch", file=sys.stderr, flush=True)
        for u_matrix, k_matrix in Parallel(
                n_jobs=cores, backend="threading", return_as="generator")(
                    delayed(process_batch)(subbatch, num_colors)
                    for subbatch in subbatch_iterable(member_batch)):
            unitig_dists += u_matrix
            kmer_dists += k_matrix
        #local_matrices = Parallel(n_jobs=cores, backend="threading")(
        #    delayed(process_batch)(subbatch, num_colors)
        #        for subbatch in subbatch_iterable(member_batch))
        #matrices_u = [u for u, _ in local_matrices]
        #if matrices_u:  # only sum if non-empty
        #    unitig_dists += np.sum(np.stack(matrices_u), axis=0)
        #matrices_k = [k for _, k in local_matrices]
        #if matrices_k:
        #    kmer_dists += np.sum(np.stack(matrices_k), axis=0)
        batch = []

    for member_batch in unique_raw_batches(color_sets, batch_size):
        for u_matrix, _ in Parallel(
                n_jobs=cores, backend="threading", return_as="generator")(
                    delayed(process_batch)(subbatch, num_colors)
                    for subbatch in subbatch_iterable(member_batch)):
            unique_raw_dists += u_matrix

    #print(f"Matrices merged, started writing", file=sys.stderr, flush=True)

    with open(f"{out_path}/{dataset}_k{k}_unitig.dists.txt", "w") as f:
        for i in range(num_colors):
            for j in range(i + 1, num_colors):
                f.write(f"{names[i]}\t{names[j]}\t{unitig_dists[i, j]}\n")

    with open(f"{out_path}/{dataset}_k{k}_kmer.dists.txt", "w") as f:
        for i in range(num_colors):
            for j in range(i + 1, num_colors):
                f.write(f"{names[i]}\t{names[j]}\t{kmer_dists[i, j]}\n")

    with open(f"{out_path}/{dataset}_k{k}_uniqraw.dists.txt", "w") as f:
        for i in range(num_colors):
            for j in range(i + 1, num_colors):
                f.write(f"{names[i]}\t{names[j]}\t{unique_raw_dists[i, j]}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dumpdir",
                        help="directory with fulgor dumps (e.g. 05_dumps)")
    parser.add_argument("dataset", help="dataset name (e.g. ngono)")
    parser.add_argument("k", type=int, help="k-mer size")
    parser.add_argument("out_prefix", help="output matrix path")
    parser.add_argument("n_cores",
                        type=int,
                        help="number of cores to parallelize")

    args = parser.parse_args()

    k = args.k
    dataset = args.dataset
    dumpdir = args.dumpdir
    out_prefix = args.out_prefix
    cores = args.n_cores

    color_sets_file = f"{dumpdir}/{dataset}_k{k}.color_sets.txt"
    metadata_file = f"{dumpdir}/{dataset}_k{k}.metadata.txt"
    unitigs_file = f"{dumpdir}/{dataset}_k{k}.unitigs.fa"
    filenames_file = f"01_datasets/{dataset}.txt"

    kmer_count = count_kmers(unitigs_file, k, metadata_file)
    color_sets = load_color_sets(color_sets_file)
    build_and_compute(color_sets, unitigs_file, k, filenames_file, out_prefix,
                      dataset, kmer_count, cores, 200000)
