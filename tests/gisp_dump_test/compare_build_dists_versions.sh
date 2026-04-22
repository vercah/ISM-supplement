#!/usr/bin/env bash

set -euo pipefail

test_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "$test_dir/../.." && pwd)"

if [ -d "$repo_root/pipeline/05_dumps" ]; then
    work_dir="$repo_root/pipeline"
else
    work_dir="$repo_root"
fi

required_files=(
    "$work_dir/01_datasets/gisp.txt"
    "$work_dir/05_dumps/gisp_k31.color_sets.txt"
    "$work_dir/05_dumps/gisp_k31.metadata.txt"
    "$work_dir/05_dumps/gisp_k31.unitigs.fa"
)

for path in "${required_files[@]}"; do
    if [ ! -f "$path" ]; then
        echo "Missing required input: $path" >&2
        exit 1
    fi
done

tmp_parent="$(mktemp -d)"
trap 'rm -rf "$tmp_parent"' EXIT

mkdir -p "$tmp_parent/v1" "$tmp_parent/v2" "$tmp_parent/v3"

(
    cd "$work_dir"
    python3 ./_build_dists_v1.py 05_dumps gisp 31 "$tmp_parent/v1" 1
    python3 ./_build_dists_v2.py 05_dumps gisp 31 "$tmp_parent/v2" 1
    python3 ./_build_dists_v3.py 05_dumps gisp 31 "$tmp_parent/v3" 1
)

compare_pair() {
    local left_dir="$1"
    local right_dir="$2"
    local kind="$3"
    local left_file="$left_dir/gisp_k31_${kind}.dists.txt"
    local right_file="$right_dir/gisp_k31_${kind}.dists.txt"

    if cmp -s "$left_file" "$right_file"; then
        return 0
    fi

    echo "Mismatch for ${kind}: ${left_dir##*/} != ${right_dir##*/}" >&2
    diff -u "$left_file" "$right_file" >&2 || true
    exit 1
}

for kind in unitig kmer uniqrow; do
    compare_pair "$tmp_parent/v1" "$tmp_parent/v2" "$kind"
    compare_pair "$tmp_parent/v1" "$tmp_parent/v3" "$kind"
done

echo "gisp dist-builder outputs match across v1, v2, and v3."
