#! /usr/bin/env bash

set -e
set -o pipefail
set -u

2>&1 echo
2>&1 echo "1. Collecting timing data"
2>&1 echo

(cat */*.time | awk '!a[$0]++') \
  | {
    read -r header
    printf '%s\n' "$header"
    sort
  } \
    > final_data.time

2>&1 echo
2>&1 echo "2. Collecting results"
2>&1 echo

./_aggregate_final_runs.py ./10_runs/*.runs > final_data.tsv
