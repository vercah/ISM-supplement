#!/usr/bin/env bash
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <dataset_name>"
    echo "  dataset_name: name of a folder in experiments/data/"
    echo "  Example: $0 minigono"
    exit 1
fi

DATASET="$1"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATASETS_DIR="$(cd "$SCRIPT_DIR/../../data" && pwd)"

if [ ! -d "$DATASETS_DIR/$DATASET" ]; then
    echo "Error: dataset directory '$DATASETS_DIR/$DATASET' does not exist."
    echo "Available datasets:"
    ls "$DATASETS_DIR"
    exit 1
fi

find "$DATASETS_DIR/$DATASET" -type f -name '*.fa' -o -name '*.fasta' -o -name '*.fna' | sort > "$SCRIPT_DIR/$DATASET.txt"
echo "Created $SCRIPT_DIR/$DATASET.txt with $(wc -l < "$SCRIPT_DIR/$DATASET.txt") entries."
