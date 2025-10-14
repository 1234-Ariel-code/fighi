#!/bin/bash
#SBATCH --job-name=FIGHI_batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=logs/fighi_batch_%A_%a.out
#SBATCH --error=logs/fighi_batch_%A_%a.err

set -euo pipefail

CSV="$1"
PHENO="$2"
OUTDIR="$3"
BATCH_LIST="$4"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${BATCH_LIST}")
COLS_FILE="${LINE}"
BASENAME=$(basename "${COLS_FILE}" .cols)

mkdir -p "${OUTDIR}/batches_out/${BASENAME}"
SUBCSV="${OUTDIR}/batches_out/${BASENAME}/${BASENAME}.csv"

# 1) extract subset CSV
srun -c "${SLURM_CPUS_PER_TASK}" python subset_columns.py \
  --csv "${CSV}" \
  --pheno "${PHENO}" \
  --cols_file "${COLS_FILE}" \
  --out_csv "${SUBCSV}"

# 2) run FIGHI on the subset
srun -c "${SLURM_CPUS_PER_TASK}" python fighi_ext/run_cli.py \
  --csv "${SUBCSV}" \
  --pheno "${PHENO}" \
  --trait binary \
  --outdir "${OUTDIR}/batches_out/${BASENAME}" \
  --max_order 4
