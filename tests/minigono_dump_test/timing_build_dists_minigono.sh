#!/usr/bin/env bash

set -euo pipefail

test_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "$test_dir/../.." && pwd)"
pipeline_dir="$repo_root/pipeline"

required_files=(
    "$pipeline_dir/05_dumps/minigono_k31.color_sets.txt"
    "$pipeline_dir/05_dumps/minigono_k31.metadata.txt"
    "$pipeline_dir/05_dumps/minigono_k31.unitigs.fa"
)

if [ ! -f "$pipeline_dir/01_datasets/minigono.txt" ]; then
    (
        cd "$pipeline_dir"
        ./prepare_minigono.sh
    )
fi

for path in "${required_files[@]}"; do
    if [ ! -f "$path" ]; then
        echo "Missing required input: $path" >&2
        echo "Prepare the minigono dump first, e.g.:" >&2
        echo "  cd $pipeline_dir" >&2
        echo "  snakemake -p -c1 05_dumps/minigono_k31.unitigs.fa" >&2
        exit 1
    fi
done

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

(
    cd "$pipeline_dir"
    "$repo_root/bin/galitime" \
        -l stderr \
        -n build_dists_minigono_k31 \
        "./_build_dists \
            05_dumps \
            minigono \
            31 \
            '$tmpdir' \
            1"
)
