#!/bin/bash
#SBATCH --job-name=FIGHI_reconcile
#SBATCH --cpus-per-task=22
#SBATCH --mem=180G
#SBATCH --time=12:00:00
#SBATCH --output=logs/fighi_reconcile_%j.out
#SBATCH --error=logs/fighi_reconcile_%j.err

set -euo pipefail

CSV="$1"         # original merged CSV
PHENO="$2"
OUTDIR="$3"      # same OUTDIR as before

UNION_LIST="${OUTDIR}/reconcile_snps.txt"
UNION_CSV="${OUTDIR}/reconcile.csv"

# 1) collect union of SNPs
python reconcile_collect_snps.py \
  --batches_out "${OUTDIR}/batches_out" \
  --out_snplist "${UNION_LIST}"

# 2) build union CSV
python subset_columns.py \
  --csv "${CSV}" \
  --pheno "${PHENO}" \
  --cols_file "${UNION_LIST}" \
  --out_csv "${UNION_CSV}"

# 3) final FIGHI pass on the union
srun -c "${SLURM_CPUS_PER_TASK}" python fighi_ext/run_cli.py \
  --csv "${UNION_CSV}" \
  --pheno "${PHENO}" \
  --trait binary \
  --outdir "${OUTDIR}/final_union" \
  --max_order 4
