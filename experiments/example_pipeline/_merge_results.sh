#!/bin/bash
#
indir="$1"
if [ -z "$indir" ]; then
  echo "Usage: $0 <directory_with_runs>"
  exit 1
fi

echo -e "genome\tcount\tk\torder\ttype\tn_bit_changes\tn_runs\tn_cylinder_runs"

for file in "$indir"/*.runs; do
  base=$(basename "$file" .runs)
  IFS='_' read -r genome order kstr ncount type <<< "$base"
  k="${kstr#k}"
  count="${ncount#N}"
  n_bit_changes=$(awk '/n_bit_changes/ {print $2}' "$file")
  n_runs=$(awk '/n_runs/ {print $2}' "$file")
  n_cylinder_runs=$(awk '/n_cylinder_runs/ {print $2}' "$file")

  echo -e "$genome\t$count\t$k\t$order\t$type\t$n_bit_changes\t$n_runs\t$n_cylinder_runs"
done

