#! /usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

find "$(realpath "$SCRIPT_DIR/../data/minigono")" -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' -o -name '*.fa.gz' -o -name '*.fasta.gz' -o -name '*.fna.gz' \) \
  | sort \
  > "$(realpath "$SCRIPT_DIR/01_datasets")/minigono.txt"
