#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATASETS_DIR="$(cd "$SCRIPT_DIR/../../data" && pwd)"

usage() {
    echo "Usage: $0 <dataset_name | directory | fasta_paths...>"
    echo "  dataset_name: name of a folder in data/"
    echo "  directory: path containing FASTA files"
    echo "  fasta_paths: one or more .fa/.fasta/.fna files, optionally .gz"
    echo "  Examples:"
    echo "    $0 minigono"
    echo "    $0 ~/tmp/neisseria_gonorrhoeae__01"
    echo "    $0 ~/tmp/neisseria_gonorrhoeae__01/SAM*.fa"
    exit 1
}

is_fasta_path() {
    case "$1" in
        *.fa|*.fasta|*.fna|*.fa.gz|*.fasta.gz|*.fna.gz) return 0 ;;
        *) return 1 ;;
    esac
}

collect_from_dir() {
    find "$1" -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' -o -name '*.fa.gz' -o -name '*.fasta.gz' -o -name '*.fna.gz' \) | sort
}

derive_dataset_name() {
    local first_path="$1"
    local common_dir

    if [ -d "$first_path" ]; then
        common_dir="$(cd "$first_path" && pwd)"
    else
        common_dir="$(cd "$(dirname "$first_path")" && pwd)"
    fi

    shift || true
    for path in "$@"; do
        local dir
        if [ -d "$path" ]; then
            dir="$(cd "$path" && pwd)"
        else
            dir="$(cd "$(dirname "$path")" && pwd)"
        fi

        while [ "$dir" != "$common_dir" ] && [ "$dir" != "/" ]; do
            common_dir="$(dirname "$common_dir")"
        done
    done

    basename "$common_dir"
}

[ $# -ge 1 ] || usage

declare -a INPUT_FILES=()

if [ $# -eq 1 ] && [ -d "$DATASETS_DIR/$1" ]; then
    DATASET="$1"
    while IFS= read -r file; do
        INPUT_FILES+=("$file")
    done < <(collect_from_dir "$DATASETS_DIR/$DATASET")
else
    DATASET="$(derive_dataset_name "$@")"
    for path in "$@"; do
        if [ -d "$path" ]; then
            while IFS= read -r file; do
                INPUT_FILES+=("$file")
            done < <(collect_from_dir "$path")
        elif [ -f "$path" ]; then
            if ! is_fasta_path "$path"; then
                echo "Error: unsupported file type '$path'." >&2
                exit 1
            fi
            INPUT_FILES+=("$(cd "$(dirname "$path")" && pwd)/$(basename "$path")")
        else
            echo "Error: '$path' does not exist." >&2
            exit 1
        fi
    done
fi

if [ ${#INPUT_FILES[@]} -eq 0 ]; then
    echo "Error: no FASTA files found." >&2
    exit 1
fi

printf '%s\n' "${INPUT_FILES[@]}" | sort -u > "$SCRIPT_DIR/$DATASET.txt"
echo "Created $SCRIPT_DIR/$DATASET.txt with $(wc -l < "$SCRIPT_DIR/$DATASET.txt") entries."
