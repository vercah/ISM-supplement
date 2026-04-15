#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
from joblib import Parallel, delayed

from _sample_name import canonical_sample_name


# Number of unitig/color-set rows accumulated before dispatching work.
# Larger values reduce scheduling overhead but do not fix dense-matrix memory use.
# For ~1k samples, 200k is reasonable. For many thousands, lower this.
DEFAULT_BATCH_SIZE = 200_000

# Number of rows per parallel joblib task.
# Should be large enough to amortize overhead, but small enough to give several
# tasks per worker. With 10 cores, 10k gives 20 tasks per 200k batch.
DEFAULT_SUBBATCH_SIZE = 10_000


class FulgorDump:
    """File-layout and parsing object for one dumped dataset/k-mer size."""

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
        """Fail early if an expected input file is missing."""
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
        """Load sample names in color-id order."""
        names = []

        with self.filenames_file.open() as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    names.append(canonical_sample_name(line))

        return names

    def _load_color_sets(self) -> dict[int, np.ndarray]:
        """Load color_set_id -> unique sample/color ids."""
        color_sets = {}

        with self.color_sets_file.open() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split(maxsplit=2)
                color_set_id = int(parts[0].split("=")[1])

                if len(parts) > 2:
                    colors = np.fromstring(parts[2], dtype=int, sep=" ")
                else:
                    colors = np.array([], dtype=int)

                # Deduplicate once here, not repeatedly during batch processing.
                color_sets[color_set_id] = np.unique(colors)

        return color_sets

    def iter_unitigs(self):
        """Yield (color_set_id, number_of_kmers_in_unitig)."""
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

    def count_kmers(self) -> int:
        """Count k-mers implied by unitigs and append num_kmers metadata if absent."""
        total_kmers = sum(n_kmers for _, n_kmers in self.iter_unitigs())

        with self.metadata_file.open("r", newline="") as f:
            already_present = any(line.startswith("num_kmers=") for line in f)

        if not already_present:
            with self.metadata_file.open("a", newline="") as f:
                f.write(f"num_kmers={total_kmers}\n")

        return total_kmers


class DistanceBuilder:
    """Build unitig-, kmer-, and unique-row-weighted distance matrices."""

    def __init__(
        self,
        dump: FulgorDump,
        cores: int,
        batch_size: int = DEFAULT_BATCH_SIZE,
        subbatch_size: int = DEFAULT_SUBBATCH_SIZE,
    ):
        self.dump = dump
        self.cores = cores
        self.batch_size = batch_size
        self.subbatch_size = subbatch_size

        # Dense matrices dominate memory:
        # memory per matrix = num_colors^2 * 8 bytes.
        # For 1102 samples: ~9.7 MB per matrix.
        # For 10k samples: ~800 MB per matrix, so this design becomes heavy.
        n = dump.num_colors
        self.unitig_dists = np.zeros((n, n), dtype=np.int64)
        self.kmer_dists = np.zeros((n, n), dtype=np.int64)
        self.uniqrow_dists = np.zeros((n, n), dtype=np.int64)

    def _subbatches(self, batch):
        """Split one large batch into joblib tasks."""
        for i in range(0, len(batch), self.subbatch_size):
            yield batch[i:i + self.subbatch_size]

    @staticmethod
    def _process_batch(batch, num_colors):
        """Compute local unitig/kmer contributions for one subbatch."""
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

        return local_unitig, local_kmer

    def _preprocess_unitig_batch(self, batch):
        """Replace color_set_id by the corresponding sample/color membership array."""
        return [
            (self.dump.color_sets.get(color_set_id, np.array([], dtype=int)), n_kmers)
            for color_set_id, n_kmers in batch
        ]

    def _parallel_batches(self, member_batch):
        """Run subbatches in parallel and stream local matrices back."""
        return Parallel(
            n_jobs=self.cores,
            backend="threading",
            return_as="generator",
        )(
            delayed(self._process_batch)(subbatch, self.dump.num_colors)
            for subbatch in self._subbatches(member_batch)
        )

    def _accumulate_weighted_batch(self, batch):
        """Accumulate unitig- and kmer-weighted distances for one batch."""
        member_batch = self._preprocess_unitig_batch(batch)

        for unitig_matrix, kmer_matrix in self._parallel_batches(member_batch):
            self.unitig_dists += unitig_matrix
            self.kmer_dists += kmer_matrix

    def compute_weighted_matrices(self):
        """Compute unitig and kmer distance matrices from unitig records."""
        batch = []

        for color_set_id, n_kmers in self.dump.iter_unitigs():
            batch.append((color_set_id, n_kmers))

            if len(batch) >= self.batch_size:
                self._accumulate_weighted_batch(batch)
                batch = []

        if batch:
            self._accumulate_weighted_batch(batch)

    def _accumulate_uniqrow_batch(self, batch):
        """Accumulate distances where each unique dumped color set counts once."""
        for uniqrow_matrix, _ in self._parallel_batches(batch):
            self.uniqrow_dists += uniqrow_matrix

    def compute_uniqrow_matrix(self):
        """Compute unique-row distances from unique dumped color sets."""
        batch = []

        for colors in self.dump.color_sets.values():
            batch.append((colors, 1))

            if len(batch) >= self.batch_size:
                self._accumulate_uniqrow_batch(batch)
                batch = []

        if batch:
            self._accumulate_uniqrow_batch(batch)

    def write_matrix(self, kind: str, matrix: np.ndarray):
        """Write upper-triangular pairwise distances."""
        path = self.dump.output_file(kind)

        with path.open("w") as f:
            for i in range(self.dump.num_colors):
                for j in range(i + 1, self.dump.num_colors):
                    f.write(
                        f"{self.dump.names[i]}\t"
                        f"{self.dump.names[j]}\t"
                        f"{matrix[i, j]}\n"
                    )

    def write_all(self):
        self.write_matrix("unitig", self.unitig_dists)
        self.write_matrix("kmer", self.kmer_dists)
        self.write_matrix("uniqrow", self.uniqrow_dists)

    def run(self):
        self.dump.count_kmers()
        self.compute_weighted_matrices()
        self.compute_uniqrow_matrix()
        self.write_all()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dumpdir", help="directory with fulgor dumps, e.g. 05_dumps")
    parser.add_argument("dataset", help="dataset name, e.g. ngono")
    parser.add_argument("k", type=int, help="k-mer size")
    parser.add_argument("out_prefix", help="output matrix directory")
    parser.add_argument("n_cores", type=int, help="number of cores to parallelize")

    # Optional tuning knobs. Defaults are reasonable for ~1k samples.
    # For 10k samples, dense matrices dominate memory; reducing cores matters
    # more than reducing batch size.
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--subbatch-size", type=int, default=DEFAULT_SUBBATCH_SIZE)

    args = parser.parse_args()

    dump = FulgorDump(
        dumpdir=args.dumpdir,
        dataset=args.dataset,
        k=args.k,
        outdir=args.out_prefix,
    )

    builder = DistanceBuilder(
        dump=dump,
        cores=args.n_cores,
        batch_size=args.batch_size,
        subbatch_size=args.subbatch_size,
    )

    builder.run()


if __name__ == "__main__":
    main()
