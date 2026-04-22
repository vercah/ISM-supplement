#!/usr/bin/env python3

# ChatGPT rewrite of the original distance-building script.
# This is the optimized rewrite (v3).

import argparse
import os
from pathlib import Path

import numpy as np

from _sample_name import canonical_sample_name


# Target default: 10 GB RAM, 10 available threads, ~10k samples.
#
# The old implementation created dense N x N matrices inside each worker.
# For 10k samples, one int64 matrix is ~800 MB, so 10 workers are impossible.
#
# This implementation computes distances blockwise:
#   block memory ~= block_size x num_samples x 8 bytes
#
# For 10k samples and block_size=256:
#   one dense block ~= 256 x 10,000 x 8 ~= 20 MB
#
# Larger blocks are faster but use more RAM.
DEFAULT_BLOCK_SIZE = 256


class FulgorDump:
    """Paths, metadata, and parsers for one Fulgor dump."""

    def __init__(self, dumpdir, dataset: str, k: int, outdir):
        self.dumpdir = Path(dumpdir)
        self.dataset = dataset
        self.k = k
        self.outdir = Path(outdir)

        self.color_sets_file = self.dumpdir / f"{dataset}_k{k}.color_sets.txt"
        self.metadata_file = self.dumpdir / f"{dataset}_k{k}.metadata.txt"
        self.unitigs_file = self.dumpdir / f"{dataset}_k{k}.unitigs.fa"
        self.filenames_file = Path("01_datasets") / f"{dataset}.txt"

        self.validate_files()
        self.outdir.mkdir(parents=True, exist_ok=True)

        self.names = self._load_names()
        self.num_colors = len(self.names)
        self.color_sets = self._load_color_sets()

    def validate_files(self) -> None:
        for path in [
            self.color_sets_file,
            self.metadata_file,
            self.unitigs_file,
            self.filenames_file,
        ]:
            if not path.exists():
                raise FileNotFoundError(path)

    def output_file(self, kind: str) -> Path:
        return self.outdir / f"{self.dataset}_k{self.k}_{kind}.dists.txt"

    def _load_names(self) -> list[str]:
        names = []

        with self.filenames_file.open() as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    names.append(canonical_sample_name(line))

        return names

    def _load_color_sets(self) -> dict[int, np.ndarray]:
        """Load color_set_id -> sample/color ids.

        Color sets are preloaded because they are much smaller than unitigs.
        Duplicate colors are removed here; duplicates would have no effect on
        the old row[colors] = 1 logic either.
        """
        color_sets = {}

        with self.color_sets_file.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split(maxsplit=2)
                color_set_id = int(parts[0].split("=")[1])

                if len(parts) > 2:
                    colors = np.fromstring(parts[2], dtype=np.int32, sep=" ")
                    colors = np.unique(colors)
                else:
                    colors = np.array([], dtype=np.int32)

                color_sets[color_set_id] = colors

        return color_sets

    def iter_unitigs(self):
        """Stream unitigs; do not load the unitig file into memory."""
        with self.unitigs_file.open() as f:
            color_set_id = None

            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith(">"):
                    parts = line[1:].split()
                    color_set_id = int(parts[1].split("=", 1)[1])
                    continue

                n_kmers = max(0, len(line) - self.k + 1)
                if n_kmers > 0:
                    yield color_set_id, n_kmers

    def append_num_kmers_if_missing(self, total_kmers: int) -> None:
        with self.metadata_file.open("r", newline="") as f:
            already_present = any(line.startswith("num_kmers=") for line in f)

        if not already_present:
            with self.metadata_file.open("a", newline="") as f:
                f.write(f"num_kmers={total_kmers}\n")


class ColorSetTable:
    """Sparse color-set incidence table plus per-color-set weights."""

    def __init__(self, dump: FulgorDump):
        self.dump = dump
        self.color_sets = dump.color_sets
        self.num_color_sets = len(self.color_sets)
        self.num_samples = dump.num_colors

        # Map arbitrary color_set_id values to dense sparse-matrix row ids.
        self.color_set_ids = list(self.color_sets.keys())
        self.row_of_color_set_id = {
            color_set_id: row
            for row, color_set_id in enumerate(self.color_set_ids)
        }

        self.uniqrow_weights = np.ones(self.num_color_sets, dtype=np.int64)
        self.unitig_weights = np.zeros(self.num_color_sets, dtype=np.int64)
        self.kmer_weights = np.zeros(self.num_color_sets, dtype=np.int64)

        self.total_kmers = 0

        self._stream_unitigs_and_accumulate_weights()
        self.dump.append_num_kmers_if_missing(self.total_kmers)

        self.incidence = self._build_sparse_incidence_matrix()

    def _stream_unitigs_and_accumulate_weights(self) -> None:
        """Stream unitigs once and aggregate weights by color_set_id.

        This preserves the original memory behavior for large unitig files:
        unitigs are never stored.
        """
        for color_set_id, n_kmers in self.dump.iter_unitigs():
            self.total_kmers += n_kmers

            row = self.row_of_color_set_id.get(color_set_id)
            if row is None:
                # Matches old behavior: missing color set means empty membership,
                # which contributes nothing to pairwise distances.
                continue

            self.unitig_weights[row] += 1
            self.kmer_weights[row] += n_kmers

    def _build_sparse_incidence_matrix(self):
        """Build X where X[color_set, sample] = 1.

        Uses scipy.sparse to avoid dense color_set x sample storage.
        """
        from scipy.sparse import csr_matrix

        rows = []
        cols = []

        for row, color_set_id in enumerate(self.color_set_ids):
            colors = self.color_sets[color_set_id]
            if len(colors) == 0:
                continue

            rows.extend([row] * len(colors))
            cols.extend(colors.tolist())

        data = np.ones(len(rows), dtype=np.int8)

        return csr_matrix(
            (data, (rows, cols)),
            shape=(self.num_color_sets, self.num_samples),
            dtype=np.int8,
        )


class BlockwiseDistanceWriter:
    """Compute weighted binary Hamming distances and write old text format."""

    def __init__(
        self,
        dump: FulgorDump,
        table: ColorSetTable,
        block_size: int = DEFAULT_BLOCK_SIZE,
    ):
        self.dump = dump
        self.table = table
        self.block_size = block_size

    def write_all(self) -> None:
        self.write_distances("unitig", self.table.unitig_weights)
        self.write_distances("kmer", self.table.kmer_weights)
        self.write_distances("uniqrow", self.table.uniqrow_weights)

    def write_distances(self, kind: str, weights: np.ndarray) -> None:
        """Write exactly the previous format:

        sample_i<TAB>sample_j<TAB>distance
        with i increasing and j > i.
        """
        X = self.table.incidence

        # weighted_X[r, i] = weights[r] * X[r, i]
        weighted_X = X.multiply(weights[:, None]).tocsr()

        # marginal[i] = sum_r weights[r] * X[r, i]
        marginal = np.asarray(weighted_X.sum(axis=0)).ravel().astype(np.int64)

        out_path = self.dump.output_file(kind)

        with out_path.open("w") as out:
            for i0 in range(0, self.table.num_samples, self.block_size):
                i1 = min(i0 + self.block_size, self.table.num_samples)

                # cooc[a, b] = sum_r weights[r] * X[r, a] * X[r, b]
                # Shape: (i1 - i0) x num_samples.
                cooc = X[:, i0:i1].T @ weighted_X
                cooc = cooc.toarray().astype(np.int64, copy=False)

                # XOR identity:
                # dist(a, b) = marginal[a] + marginal[b] - 2 * cooc[a, b]
                dist_block = (
                    marginal[i0:i1, None]
                    + marginal[None, :]
                    - 2 * cooc
                )

                self._write_upper_triangle_block(out, dist_block, i0, i1)

    def _write_upper_triangle_block(self, out, dist_block: np.ndarray, i0: int, i1: int):
        names = self.dump.names
        n = self.table.num_samples

        for local_i, i in enumerate(range(i0, i1)):
            for j in range(i + 1, n):
                out.write(f"{names[i]}\t{names[j]}\t{dist_block[local_i, j]}\n")


def configure_threads(n_threads: int) -> None:
    """Best-effort cap for numerical libraries.

    This implementation does not create Python worker-local dense matrices.
    Sparse operations may still use BLAS/OpenMP threads depending on scipy.
    """
    os.environ.setdefault("OMP_NUM_THREADS", str(n_threads))
    os.environ.setdefault("OPENBLAS_NUM_THREADS", str(n_threads))
    os.environ.setdefault("MKL_NUM_THREADS", str(n_threads))
    os.environ.setdefault("NUMEXPR_NUM_THREADS", str(n_threads))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dumpdir", help="directory with fulgor dumps, e.g. 05_dumps")
    parser.add_argument("dataset", help="dataset name, e.g. ngono")
    parser.add_argument("k", type=int, help="k-mer size")
    parser.add_argument("out_prefix", help="output matrix directory")
    parser.add_argument("n_cores", type=int, help="thread cap for numerical libraries")

    parser.add_argument(
        "--block-size",
        type=int,
        default=DEFAULT_BLOCK_SIZE,
        help=(
            "number of sample rows computed at once; default 256 is conservative "
            "for ~10 GB RAM and ~10k samples"
        ),
    )

    args = parser.parse_args()

    configure_threads(args.n_cores)

    dump = FulgorDump(
        dumpdir=args.dumpdir,
        dataset=args.dataset,
        k=args.k,
        outdir=args.out_prefix,
    )

    table = ColorSetTable(dump)

    writer = BlockwiseDistanceWriter(
        dump=dump,
        table=table,
        block_size=args.block_size,
    )

    writer.write_all()


if __name__ == "__main__":
    main()
