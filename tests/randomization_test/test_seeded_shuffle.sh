#!/usr/bin/env bash

set -euo pipefail

test_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "$test_dir/../.." && pwd)"
shuffle_script="$repo_root/pipeline/_shuffle_lines_with_seed.py"
input_file="$test_dir/input.txt"
tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

same_a="$tmp_dir/same_a.txt"
same_b="$tmp_dir/same_b.txt"
different_seed="$tmp_dir/different_seed.txt"

"$shuffle_script" \
    --seed 123 \
    "$input_file" \
    "$same_a"
"$shuffle_script" \
    --seed 123 \
    "$input_file" \
    "$same_b"
"$shuffle_script" \
    --seed 124 \
    "$input_file" \
    "$different_seed"

cmp -s "$same_a" "$same_b"

if cmp -s "$same_a" "$different_seed"; then
    echo "Different seeds produced the same shuffled output." >&2
    exit 1
fi

if ! diff -u <(LC_ALL=C sort "$input_file") <(LC_ALL=C sort "$same_a"); then
    echo "Shuffle output is not a permutation of the input." >&2
    exit 1
fi

echo "Seeded shuffle helper passed reproducibility checks."
