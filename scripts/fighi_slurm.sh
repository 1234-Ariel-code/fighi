#!/usr/bin/env bash
#SBATCH --job-name=fighi
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=fighi-%j.log

set -euo pipefail

: "${FIGHI_INPUT:?Set FIGHI_INPUT to a CSV or TSV path}"
: "${FIGHI_PHENOTYPE:?Set FIGHI_PHENOTYPE to the phenotype column}"
: "${FIGHI_OUTDIR:?Set FIGHI_OUTDIR to an output directory}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/fighi-matplotlib-${SLURM_JOB_ID:-local}"

fighi run "$FIGHI_INPUT" \
  --phenotype "$FIGHI_PHENOTYPE" \
  --id-column "${FIGHI_ID_COLUMN:-auto}" \
  --covariates "${FIGHI_COVARIATES:-}" \
  --trait "${FIGHI_TRAIT:-auto}" \
  --top-m "${FIGHI_TOP_M:-50}" \
  --max-order "${FIGHI_MAX_ORDER:-2}" \
  --outdir "$FIGHI_OUTDIR"

