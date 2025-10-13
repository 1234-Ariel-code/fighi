#!/bin/bash
#SBATCH --job-name=FIGHITPED
#SBATCH --cpus-per-task=22
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fighitped_%j.out
#SBATCH --error=logs/fighitped_%j.err
# #SBATCH --partition=compute          # uncomment & set if your cluster requires a partition
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1


set -euo pipefail

WORKDIR="/work/long_lab/for_Ariel/files"
TPED="${WORKDIR}/CD_origin.tped"
TFAM="${WORKDIR}/CD_origin.tfam"          # optional; comment if not available
PHEN_FILE="${WORKDIR}/CD_origin.phen"

GENO_NAMED="${WORKDIR}/CD_filtered_named.csv"
MERGED_CSV="${WORKDIR}/CD_merged.csv"
PHENO_NAME="case"
OUTDIR="${WORKDIR}/fighi_out"

FIGHI_SRC="${WORKDIR}/fighi_ext"          # leave blank if installed site-wide
USE_CONDA=1
CONDA_ENV="fighi"

module purge || true
module load python/3.10 || module load python || true

if [[ "${USE_CONDA}" -eq 1 ]]; then
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV}" || conda create -y -n "${CONDA_ENV}" python=3.10
    conda activate "${CONDA_ENV}"
  else
    echo "[WARN] conda not found; using system Python."
  fi
fi

python - <<'PY'
import importlib, sys, subprocess
for pkg in ["numpy","pandas","scipy"]:
    try: importlib.import_module(pkg)
    except Exception: subprocess.check_call([sys.executable,"-m","pip","install","-q",pkg])
PY

[[ -n "${FIGHI_SRC}" && -d "${FIGHI_SRC}" ]] && export PYTHONPATH="${FIGHI_SRC}:${PYTHONPATH:-}"
mkdir -p "${OUTDIR}"

#echo "== Step 0: TPED -> named numeric CSV =="
#if [[ -f "${TPED}" ]]; then
#  if [[ -f "${TFAM}" ]]; then
#    srun python tped_to_named_csv.py \
#      --tped "${TPED}" --tfam "${TFAM}" \
#      --phen_file "${PHEN_FILE}" \
#      --out_csv "${GENO_NAMED}" --ref_mode major
#  else
#    echo "[INFO] No TFAM; using IID order from .phen."
#    srun python tped_to_named_csv.py \
#      --tped "${TPED}" \
#      --phen_file "${PHEN_FILE}" \
#      --out_csv "${GENO_NAMED}" --ref_mode major
#  fi
#else
#  echo "[ERROR] Missing TPED: ${TPED}"; exit 2
#fi

#echo "== Step 1: Merge with .phen (no pandas) =="
#srun python merge_pheno_geno_nopandas.py \
#  --geno_csv "${GENO_NAMED}" \
#  --phen_file "${PHEN_FILE}" \
#  --out_csv "${MERGED_CSV}" \
#  --phen_name "${PHENO_NAME}" \
#  --impute_mean


echo "== Step 2: Run FIGHI =="
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
srun python -m fighi_ext.run_cli \
  --csv "${MERGED_CSV}" \
  --pheno "${PHENO_NAME}" \
  --trait binary \
  --outdir "${OUTDIR}" \
  --kmax 4

echo "== Done. Outputs in: ${OUTDIR} =="

