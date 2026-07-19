#!/usr/bin/env bash
#SBATCH --job-name=fighi-prepare
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=fighi-prepare-%j.log

set -euo pipefail

: "${FIGHI_ENV_PREFIX:?Set FIGHI_ENV_PREFIX to the Conda environment prefix}"
: "${FIGHI_PLINK_INPUT:?Set FIGHI_PLINK_INPUT to a BED/PGEN prefix or VCF path}"
: "${FIGHI_CANDIDATES:?Set FIGHI_CANDIDATES to the pre-specified variant list}"
: "${FIGHI_PHENOTYPE_FILE:?Set FIGHI_PHENOTYPE_FILE to the IID-keyed phenotype table}"
: "${FIGHI_PREPARED_DIR:?Set FIGHI_PREPARED_DIR to a new output directory}"

CONDA_RUN=(conda run --prefix "$FIGHI_ENV_PREFIX")
THREADS="${SLURM_CPUS_PER_TASK:-1}"
MEMORY_MB="${SLURM_MEM_PER_NODE:-16000}"

args=(
  prepare-plink
  --input "$FIGHI_PLINK_INPUT"
  --input-type "${FIGHI_PLINK_INPUT_TYPE:-auto}"
  --candidates "$FIGHI_CANDIDATES"
  --phenotype-file "$FIGHI_PHENOTYPE_FILE"
  --phenotype-column "${FIGHI_PHENOTYPE_COLUMN:-phenotype}"
  --iid-column "${FIGHI_IID_COLUMN:-IID}"
  --threads "$THREADS"
  --memory-mb "$MEMORY_MB"
  --outdir "$FIGHI_PREPARED_DIR"
)

if [[ -n "${FIGHI_KEEP_FILE:-}" ]]; then
  args+=(--samples "$FIGHI_KEEP_FILE")
fi
if [[ -n "${FIGHI_COVARIATES:-}" ]]; then
  args+=(--covariates "$FIGHI_COVARIATES")
fi

"${CONDA_RUN[@]}" fighi "${args[@]}"
