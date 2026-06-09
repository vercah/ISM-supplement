#!/usr/bin/env python3
import sys
import lzma
import argparse
from collections import OrderedDict
from _sample_name import canonical_sample_name


def main():
  parser = argparse.ArgumentParser(
      description="Convert pairwise distances (.txt.xz) to PHYLIP distance matrix"
  )
  parser.add_argument("dists", help="Pairwise distance file (.txt.xz or .txt)")
  parser.add_argument("selection", help="Selection file with genome names (determines order)")
  parser.add_argument("-o", "--output", default="-", help="Output file (default: stdout)")
  args = parser.parse_args()

  # Load genome names from selection file (preserves order). Use canonical
  # form so they match what `_build_dists` writes into the dists file.
  names = []
  with open(args.selection) as f:
      for line in f:
          line = line.strip()
          if line and not line.startswith("#"):
              names.append(canonical_sample_name(line))

  name_to_idx = {n: i for i, n in enumerate(names)}
  n = len(names)

  # Initialize symmetric distance matrix with zeros on diagonal
  matrix = [[0.0] * n for _ in range(n)]

  # Read pairwise distances
  opener = lzma.open if args.dists.endswith(".xz") else open
  with opener(args.dists, "rt") as f:
      for line in f:
          parts = line.strip().split("\t")
          if len(parts) != 3:
              continue
          a, b, d = parts[0], parts[1], float(parts[2])
          if a in name_to_idx and b in name_to_idx:
              i, j = name_to_idx[a], name_to_idx[b]
              matrix[i][j] = d
              matrix[j][i] = d

  # Strict PHYLIP truncates labels to 10 chars, but quicktree v2.5 reads
  # names whitespace-delimited up to 100 chars (see distancemat.c), so we
  # can emit full canonical names with two spaces before the distances.
  out = sys.stdout if args.output == "-" else open(args.output, "w")
  out.write(f"    {n}\n")
  for i, name in enumerate(names):
      row = "  ".join(f"{matrix[i][j]:.10f}" for j in range(n))
      out.write(f"{name}  {row}\n")

  if out is not sys.stdout:
      out.close()


if __name__ == "__main__":
    main()
