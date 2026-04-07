#! /usr/bin/env bash

set -e
set -o pipefail
set -u

find "$(realpath ../data/minigono)" -name '*.fa' \
  | sort \
  > "$(realpath 01_datasets)/minigono.txt"
