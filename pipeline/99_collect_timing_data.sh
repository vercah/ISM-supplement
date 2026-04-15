#! /usr/bin/env bash

set -e
set -o pipefail
set -u

(cat */*.time | awk '!a[$0]++') \
  | {
    read -r header
    printf '%s\n' "$header"
    sort
  } \
    | tee timing_data.tsv
