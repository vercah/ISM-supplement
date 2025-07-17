#!/bin/bash

set -euo pipefail

cp -v 09_runs/* ../old_results_backup/28_pipeline_single_tree/09_runs/

./_merge_results.sh > 10_results/merged.tsv
echo "merged results in 09_runs"
cp -v 10_results/merged.tsv ../old_results_backup/28_pipeline_single_tree/merged.tsv
