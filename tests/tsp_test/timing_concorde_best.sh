#! /usr/bin/env bash

set -e
set -o pipefail
set -u

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "$script_dir/../.." && pwd)"
d=$(mktemp -d)

(cd "$d"; python3 "$repo_root/bin/galitime.py" -l stderr -n timing_concorde_best "'$repo_root/bin/concorde' '$script_dir/instance_ngono_best/instance.tsp'")
