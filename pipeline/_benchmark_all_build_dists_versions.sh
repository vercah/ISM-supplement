#!/usr/bin/env bash

# Defaults correspond to:
#   ./_build_dists_v1.py 05_dumps gisp 31 06_ham_dists 10

set -euo pipefail

PIPELINE_DIR="${1:-.}"
INPUT_DIR="${2:-05_dumps}"
DATASET="${3:-gisp}"
K="${4:-31}"
OUT_BASE="${5:-06_ham_dists}"
THREADS="${6:-10}"

cd "$PIPELINE_DIR"

mkdir -p __aux_timing_logs

make -C _build_dists_v4
make -C _build_dists_v5
make -C _build_dists_v6

run_one() {
    local version="$1"
    local cmd="$2"
    local outdir="${OUT_BASE}/v${version}"

    mkdir -p "$outdir"

    ../bin/galitime \
        -l stdout \
        -n "build_dists_v${version}_${DATASET}_k${K}" \
        "$cmd"
}

(
run_one 1 "./_build_dists_v1.py ${INPUT_DIR} ${DATASET} ${K} ${OUT_BASE}/v1 ${THREADS}"
run_one 2 "./_build_dists_v2.py ${INPUT_DIR} ${DATASET} ${K} ${OUT_BASE}/v2 ${THREADS}"
run_one 3 "./_build_dists_v3.py ${INPUT_DIR} ${DATASET} ${K} ${OUT_BASE}/v3 ${THREADS}"
run_one 4 "./_build_dists_v4/_build_dists_v4 ${INPUT_DIR} ${DATASET} ${K} ${OUT_BASE}/v4 ${THREADS}"
run_one 5 "./_build_dists_v5/_build_dists_v5 ${INPUT_DIR} ${DATASET} ${K} ${OUT_BASE}/v5 ${THREADS}"
run_one 6 "./_build_dists_v6/_build_dists_v6 ${INPUT_DIR} ${DATASET} ${K} ${OUT_BASE}/v6 ${THREADS}"
) \
    | awk '!a[$0]++' \
    | tee galitime_build_dists_all.tsv
