#!/usr/bin/env bash
set -euo pipefail

# Configuration
EPSILONS=(1 2 3)
QS=(4 5 6 7 8)
P=2
NUM_ITERS=10
# Algorithms to test: 0 (Naive), 2, 3, 4
ALGORITHMS=(0 2 3 4)
ALPHA0=0.05
ALPHA1=0.60
ALPHA2=0.35
# 11th argument: probability filtering (0=disabled, 1=enabled)
PROB_FILTER=0

# Datasets to test (paths relative to this script's directory)
DATASETS=(
  "librec-filmtrust-ratings"
  "csbwiki"
  "lrcwiki"
  "rmwiki"
)

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
BIN="$ROOT_DIR/biclique"
DATA_ROOT="${ROOT_DIR}/../bidata"
OUT_ROOT="${ROOT_DIR}/results_vary_eps_q_p2"

mkdir -p "$OUT_ROOT"

for ds in "${DATASETS[@]}"; do
  DATA_PATH="${DATA_ROOT}/${ds}"
  # Accept either a directory dataset or a single file dataset with .e extension
  if [[ -d "$DATA_PATH" ]]; then
    RESOLVED_PATH="$DATA_PATH"
  elif [[ -f "${DATA_PATH}.e" ]]; then
    RESOLVED_PATH="${DATA_PATH}"
  else
    echo "[WARN] Dataset not found: $DATA_PATH (looked for dir or ${DATA_PATH}.e) (skipping)" >&2
    continue
  fi

  for eps in "${EPSILONS[@]}"; do
    for q in "${QS[@]}"; do
      for algo in "${ALGORITHMS[@]}"; do
        OUT_DIR="${OUT_ROOT}/${ds}/eps${eps}"
        mkdir -p "$OUT_DIR"
        LOG_FILE="${OUT_DIR}/p${P}_q${q}_algo${algo}.log"

        echo "Running: eps=${eps} ds=${ds} p=${P} q=${q} algo=${algo} iters=${NUM_ITERS} prob_filter=${PROB_FILTER}"
        "${BIN}" "${eps}" "${RESOLVED_PATH}" "${NUM_ITERS}" "${algo}" "${P}" "${q}" \
          "${ALPHA0}" "${ALPHA1}" "${ALPHA2}" "${PROB_FILTER}" | tee "${LOG_FILE}"
      done
    done
  done
done

echo "All runs completed. Logs in ${OUT_ROOT}" 


