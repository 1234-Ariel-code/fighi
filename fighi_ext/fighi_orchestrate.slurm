#!/bin/bash
#SBATCH --job-name=FIGHI_orchestrate
#SBATCH --cpus-per-task=22
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fighi_orch_%j.out
#SBATCH --error=logs/fighi_orch_%j.err

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


set -euo pipefail

# ---- USER INPUTS ----
WORKDIR="/work/long_lab/for_Ariel/files"
MERGED_CSV="${WORKDIR}/CD_merged.csv"
PHENO_NAME="case"
OUTDIR="${WORKDIR}/fighi_out_batched"

# prescreen settings
TOP_M=5000
COL_BLOCK=800
ROW_CHUNK=50000

# batching settings
BATCH_SIZE=1000
HALO=200
ANCHORS=200

mkdir -p "${OUTDIR}" logs

echo "[1/5] Prescreen streaming -> top pool"
srun -c "${SLURM_CPUS_PER_TASK}" python run_cli.py \
  --csv "${MERGED_CSV}" \
  --pheno "${PHENO_NAME}" \
  --trait binary \
  --outdir "${OUTDIR}/prescreen" \
  --prescreen_top_m "${TOP_M}" \
  --col_block "${COL_BLOCK}" \
  --row_chunksize "${ROW_CHUNK}" \
  --max_order 2 \
  --write_top_cols

echo "[2/5] Make overlapping batches from top pool"
python prep_batches.py \
  --top_cols "${OUTDIR}/prescreen/top_columns.txt" \
  --batch_size "${BATCH_SIZE}" \
  --halo "${HALO}" \
  --anchors "${ANCHORS}" \
  --outdir "${OUTDIR}/batches"

BATCH_LIST="${OUTDIR}/batches/batches.txt"
NUM_BATCHES=$(wc -l < "${BATCH_LIST}")
echo "[INFO] ${NUM_BATCHES} batches prepared"

echo "[3/5] Submit array for per-batch FIGHI"
BATCH_JOBID=$(sbatch --parsable --array=1-${NUM_BATCHES} fighi_batch.slurm "${MERGED_CSV}" "${PHENO_NAME}" "${OUTDIR}" "${BATCH_LIST}")
echo "[INFO] batch array job id: ${BATCH_JOBID}"

echo "[4/5] Submit reconcile dependent on array completion"
sbatch --dependency=afterok:${BATCH_JOBID} --parsable reconcile_run.slurm "${MERGED_CSV}" "${PHENO_NAME}" "${OUTDIR}"
echo "[INFO] submitted reconcile job (afterok dependency)"

echo "[5/5] Orchestration submitted. Check logs/ for progress."
