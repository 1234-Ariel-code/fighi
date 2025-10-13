#!/bin/bash
#SBATCH --job-name=FIGHI_full
#SBATCH --cpus-per-task=22
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fighi_%j.out
#SBATCH --error=logs/fighi_%j.err
# #SBATCH --partition=compute
# SBATCH --partition=gpu-a100,theia
# SBATCH --gres=gpu:1

set -euo pipefail

# --------------------------
# 1) Global knobs (set once)
# --------------------------
DISEASE_NAME="CD"              # e.g., "UC", "T2D", etc.
PHENO_NAME="case"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# --------------------------
# 2) Paths
# --------------------------
# Data lives here (tped/tfam/phen, etc.)
WORKDIR="/work/long_lab/for_Ariel/files"
TPED="${WORKDIR}/${DISEASE_NAME}_origin.tped"
TFAM="${WORKDIR}/${DISEASE_NAME}_origin.tfam"         # optional
PHEN_FILE="${WORKDIR}/${DISEASE_NAME}_origin.phen"

# Intermediates in DATA area (ok) – or change to CODE_DIR if you prefer
GENO_NAMED="${WORKDIR}/${DISEASE_NAME}_filtered_named.csv"
MERGED_CSV="${WORKDIR}/${DISEASE_NAME}_merged.csv"

# Code directory = where you submit the job from (where your .py files are)
CODE_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
export PYTHONPATH="${CODE_DIR}:${PYTHONPATH:-}"

# Outputs/logs live next to your code
OUTDIR="${CODE_DIR}/fighi_out"
LOGDIR="${CODE_DIR}/logs"
mkdir -p "${OUTDIR}" "${LOGDIR}"

# --------------------------
# 3) Environment
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

# tee logs into ./logs
ts="$(date +%Y%m%d_%H%M%S)"
exec > >(tee -a "${LOGDIR}/fighi_${DISEASE_NAME}_${ts}.out")
exec 2> >(tee -a "${LOGDIR}/fighi_${DISEASE_NAME}_${ts}.err" >&2)

echo "[INFO] CODE_DIR=${CODE_DIR}"
echo "[INFO] OUTDIR=${OUTDIR}"
echo "[INFO] Using data from WORKDIR=${WORKDIR}"

# --------------------------
# 4) Pipeline
# --------------------------
echo "== Step 0: TPED -> named CSV =="
if [[ -s "${GENO_NAMED}" ]]; then
  echo "[INFO] Found existing ${GENO_NAMED}; skipping TPED conversion."
else
  if [[ -f "${TPED}" ]]; then
    if [[ -f "${TFAM}" ]]; then
      srun -c "${SLURM_CPUS_PER_TASK}" python "${CODE_DIR}/tped_to_named_csv.py" \
        --tped "${TPED}" --tfam "${TFAM}" \
        --phen_file "${PHEN_FILE}" \
        --out_csv "${GENO_NAMED}" --ref_mode major
    else
      srun -c "${SLURM_CPUS_PER_TASK}" python "${CODE_DIR}/tped_to_named_csv.py" \
        --tped "${TPED}" --phen_file "${PHEN_FILE}" \
        --out_csv "${GENO_NAMED}" --ref_mode major
    fi
  else
    echo "[ERROR] Missing TPED: ${TPED}"; exit 2
  fi
fi

echo "== Step 1: Merge (no-pandas) =="
srun -c "${SLURM_CPUS_PER_TASK}" python "${CODE_DIR}/merge_pheno_geno_nopandas.py" \
  --geno_csv "${GENO_NAMED}" \
  --phen_file "${PHEN_FILE}" \
  --out_csv "${MERGED_CSV}" \
  --phen_name "${PHENO_NAME}" \
  --impute_mean

echo "== Step 2: Run FIGHI =="
srun -c "${SLURM_CPUS_PER_TASK}" python -m fighi_ext.run_cli \
  --csv "${MERGED_CSV}" \
  --pheno "${PHENO_NAME}" \
  --trait binary \
  --outdir "${OUTDIR}" \
  --max_order 4

# (Optional) Step 3: annotate SNPs → genes/pathways
# CS2G_DIR="${CODE_DIR}/cS2G_1000GEUR"
# srun -c "${SLURM_CPUS_PER_TASK}" python "${CODE_DIR}/annotate_fighi_features.py" \
#   --feature_csv "${OUTDIR}/fighi_feature_scores.csv" \
#   --cs2g_dir "${CS2G_DIR}" \
#   --out_csv "${OUTDIR}/fighi_feature_scores_annotated.csv"
#   # Or add: --gmt_files ...  or  --use_gprofiler --organism hsapiens

echo "== Done. Outputs in: ${OUTDIR} =="

