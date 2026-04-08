#! /usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

find "$(realpath "$SCRIPT_DIR/../data/minigono")" -name '*.fa' \
  | sort \
  > "$(realpath "$SCRIPT_DIR/01_datasets")/minigono.txt"
