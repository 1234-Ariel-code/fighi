#!/bin/bash
#SBATCH --job-name=FIGHI_local
#SBATCH --cpus-per-task=22
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fighi_%j.out
#SBATCH --error=logs/fighi_%j.err
# SBATCH --partition=gpu-a100,mtst
# SBATCH --gres=gpu:0

set -euo pipefail

# --------------------------
# 1. Disease setup
# --------------------------
DISEASE_NAME="CD"         # set once
PHENO_NAME="case"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# --------------------------
# 2. Paths
# --------------------------
# NOTE: Data is still in /work/... but all code & outputs stay here
DATA_DIR="/work/long_lab/for_Ariel/files"
TPED="${DATA_DIR}/${DISEASE_NAME}_origin.tped"
TFAM="${DATA_DIR}/${DISEASE_NAME}_origin.tfam"
PHEN_FILE="${DATA_DIR}/${DISEASE_NAME}_origin.phen"
GENO_NAMED="${DATA_DIR}/${DISEASE_NAME}_filtered_named.csv"
MERGED_CSV="${DATA_DIR}/${DISEASE_NAME}_merged.csv"

# Everything else (code, logs, outputs) stays local
CODE_DIR="$PWD"
OUTDIR="${CODE_DIR}/fighi_out"
LOGDIR="${CODE_DIR}/logs"
mkdir -p "${OUTDIR}" "${LOGDIR}"

# --------------------------
# 3. Environment
# --------------------------
USE_CONDA=1
CONDA_ENV="fighi"

module purge || true
module load python/3.10 || module load python || true

if [[ "${USE_CONDA}" -eq 1 ]]; then
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV}" || conda create -y -n "${CONDA_ENV}" python=3.10
    conda activate "${CONDA_ENV}"
  fi
fi

python - <<'PY'
import importlib, sys, subprocess
for pkg in ["numpy","pandas","scipy","matplotlib","statsmodels"]:
    try: importlib.import_module(pkg)
    except Exception: subprocess.check_call([sys.executable,"-m","pip","install","-q",pkg])
PY

# Add current folder (fighi_ext) to sys.path so imports work
export PYTHONPATH="${CODE_DIR}:${PYTHONPATH:-}"

echo "[INFO] Running FIGHI inside: ${CODE_DIR}"
echo "[INFO] Using data from: ${DATA_DIR}"
echo "[INFO] Outputs: ${OUTDIR}"

# --------------------------
# 4. Pipeline
# --------------------------
echo "== Step 0: TPED -> named CSV =="
if [[ -s "${GENO_NAMED}" ]]; then
  echo "[INFO] Found existing ${GENO_NAMED}; skipping conversion."
else
  if [[ -f "${TPED}" ]]; then
    if [[ -f "${TFAM}" ]]; then
      srun -c "${SLURM_CPUS_PER_TASK}" python tped_to_named_csv.py \
        --tped "${TPED}" --tfam "${TFAM}" \
        --phen_file "${PHEN_FILE}" \
        --out_csv "${GENO_NAMED}" --ref_mode major
    else
      srun -c "${SLURM_CPUS_PER_TASK}" python tped_to_named_csv.py \
        --tped "${TPED}" --phen_file "${PHEN_FILE}" \
        --out_csv "${GENO_NAMED}" --ref_mode major
    fi
  else
    echo "[ERROR] Missing TPED: ${TPED}"; exit 2
  fi
fi

echo "== Step 1: Merge (no-pandas) =="
srun -c "${SLURM_CPUS_PER_TASK}" python merge_pheno_geno_nopandas.py \
  --geno_csv "${GENO_NAMED}" \
  --phen_file "${PHEN_FILE}" \
  --out_csv "${MERGED_CSV}" \
  --phen_name "${PHENO_NAME}" \
  --impute_mean

echo "== Step 2: Run FIGHI =="
srun -c "${SLURM_CPUS_PER_TASK}" python run_cli.py \
  --csv "${MERGED_CSV}" \
  --pheno "${PHENO_NAME}" \
  --trait binary \
  --outdir "${OUTDIR}" \
  --max_order 4

# (Optional) Step 3: annotate SNPs
# CS2G_DIR="${CODE_DIR}/cS2G_1000GEUR"
# srun -c "${SLURM_CPUS_PER_TASK}" python annotate_fighi_features.py \
#   --feature_csv "${OUTDIR}/fighi_feature_scores.csv" \
#   --cs2g_dir "${CS2G_DIR}" \
#   --out_csv "${OUTDIR}/fighi_feature_scores_annotated.csv"

echo "== Done. Outputs in: ${OUTDIR} =="
