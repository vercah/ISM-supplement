#!/bin/bash
echo -e "genome\tcount\tk\torder\tn_bit_changes\tn_row_runs"
for file in 09_runs/*.runs; do
  base=$(basename "$file" .runs)
  IFS='_' read -r genome order kstr ncount <<< "$base"
  k="${kstr#k}"
  count="${ncount#N}"
  n_bit_changes=$(awk '/n_bit_changes/ {print $2}' "$file")
  n_runs=$(awk '/n_runs/ {print $2}' "$file")
  echo -e "$genome\t$count\t$k\t$order\t$n_bit_changes\t$n_runs"
done

