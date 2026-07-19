#!/usr/bin/env bash
#SBATCH --job-name=fighi-run
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=fighi-run-%j.log

set -euo pipefail

: "${FIGHI_ENV_PREFIX:?Set FIGHI_ENV_PREFIX}"
: "${FIGHI_INPUT:?Set FIGHI_INPUT to a FIGHI-ready CSV/TSV}"
: "${FIGHI_PHENOTYPE_COLUMN:?Set FIGHI_PHENOTYPE_COLUMN}"
: "${FIGHI_RESULTS_DIR:?Set FIGHI_RESULTS_DIR to a new output directory}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/fighi-matplotlib-${SLURM_JOB_ID:-local}"

args=(
  run "$FIGHI_INPUT"
  --phenotype "$FIGHI_PHENOTYPE_COLUMN"
  --id-column "${FIGHI_IID_COLUMN:-IID}"
  --trait "${FIGHI_TRAIT:-auto}"
  --top-m "${FIGHI_TOP_M:-100}"
  --max-order "${FIGHI_MAX_ORDER:-2}"
  --screen-method "${FIGHI_SCREEN_METHOD:-hybrid}"
  --discovery-fraction "${FIGHI_DISCOVERY_FRACTION:-0.40}"
  --correction "${FIGHI_CORRECTION:-fdr_bh}"
  --alpha "${FIGHI_ALPHA:-0.05}"
  --seed "${FIGHI_SEED:-17}"
  --outdir "$FIGHI_RESULTS_DIR"
)

if [[ -n "${FIGHI_COVARIATES:-}" ]]; then
  args+=(--covariates "$FIGHI_COVARIATES")
fi
if [[ -n "${FIGHI_FEATURE_FILE:-}" ]]; then
  args+=(--feature-file "$FIGHI_FEATURE_FILE")
fi

conda run --prefix "$FIGHI_ENV_PREFIX" fighi "${args[@]}"
