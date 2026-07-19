#!/usr/bin/env bash
#SBATCH --job-name=fighi-benchmark
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=fighi-benchmark-%j.log

set -euo pipefail

: "${FIGHI_ENV_PREFIX:?Set FIGHI_ENV_PREFIX}"
: "${FIGHI_BENCHMARK_MANIFEST:?Set FIGHI_BENCHMARK_MANIFEST to the JSON manifest}"
: "${FIGHI_BENCHMARK_DIR:?Set FIGHI_BENCHMARK_DIR to a new output directory}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/fighi-matplotlib-${SLURM_JOB_ID:-local}"

conda run --prefix "$FIGHI_ENV_PREFIX" \
  fighi benchmark-validate "$FIGHI_BENCHMARK_MANIFEST"

conda run --prefix "$FIGHI_ENV_PREFIX" \
  fighi benchmark-run "$FIGHI_BENCHMARK_MANIFEST" \
  --outdir "$FIGHI_BENCHMARK_DIR" \
  ${FIGHI_BENCHMARK_STRICT:+--strict}
