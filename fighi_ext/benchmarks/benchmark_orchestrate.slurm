#!/bin/bash
#SBATCH --job-name=FIGHI_vs_PLINK_MDR
#SBATCH --cpus-per-task=22
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=benchmarks/logs/bench_%j.out
#SBATCH --error=benchmarks/logs/bench_%j.err

set -euo pipefail

# -------- USER INPUTS --------
WORKDIR="/work/long_lab/for_Ariel/files"
CSV="${WORKDIR}/CD_merged.csv"       # IID,case,<rs...>
PHENO="case"
TOP_M=3000                           # common SNP pool size after prescreen
COL_BLOCK=1000
ROW_CHUNK=50000
OUTDIR="${WORKDIR}/benchmarks/out"   # results go here
FIGHI_DIR="${WORKDIR}/fighi_ext"     # package with run_cli.py
# -----------------------------

mkdir -p "${OUTDIR}" benchmarks/logs

# --- Environment (Conda + Python + tools) ---
USE_CONDA=1
CONDA_ENV="fighi-bench"

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
pkgs = ["numpy","pandas","scipy","matplotlib","statsmodels","scikit-learn","scikit-mdr"]
for p in pkgs:
    try: importlib.import_module(p)
    except Exception: subprocess.check_call([sys.executable,"-m","pip","install","-q",p])
PY

# Download PLINK v1.9 locally if not available
if ! command -v plink >/dev/null 2>&1; then
  echo "[INFO] Downloading PLINK v1.9 locally..."
  mkdir -p "${OUTDIR}/tools"
  cd "${OUTDIR}/tools"
  # Linux x86_64 binary
  wget -qO plink.zip https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip
  unzip -o plink.zip >/dev/null
  chmod +x plink
  export PATH="${OUTDIR}/tools:${PATH}"
  cd - >/dev/null
fi

# For Python to find fighi_ext
export PYTHONPATH="${FIGHI_DIR}:${PYTHONPATH:-}"

echo "========== [1/5] Prescreen with FIGHI to fix a common SNP set =========="
# We only prescreen here; we won't use FIGHI results from this step for runtime.
# We just need a fair list of TOP_M SNPs used by all methods.
srun -c "${SLURM_CPUS_PER_TASK}" python "${FIGHI_DIR}/run_cli.py" \
  --csv "${CSV}" \
  --pheno "${PHENO}" \
  --trait binary \
  --outdir "${OUTDIR}/prescreen" \
  --prescreen_top_m "${TOP_M}" \
  --col_block "${COL_BLOCK}" \
  --row_chunksize "${ROW_CHUNK}" \
  --max_order 2 \
  --write_top_cols \
  --no_plots

TOP_COLS="${OUTDIR}/prescreen/top_columns.txt"
echo "[OK] Wrote SNP list: ${TOP_COLS}"

echo "========== [2/5] Materialize filtered CSV with only TOP_M SNPs =========="
FILTERED_CSV="${OUTDIR}/cd_top${TOP_M}.csv"
python benchmarks/filter_csv_by_columns.py \
  --csv "${CSV}" \
  --pheno "${PHENO}" \
  --keep "${TOP_COLS}" \
  --out "${FILTERED_CSV}"

echo "========== [3/5] FIGHI timing =========="
FIGHI_OUT="${OUTDIR}/fighi_run"
mkdir -p "${FIGHI_OUT}"
/usr/bin/time -v python "${FIGHI_DIR}/run_cli.py" \
  --csv "${FILTERED_CSV}" \
  --pheno "${PHENO}" \
  --trait binary \
  --outdir "${FIGHI_OUT}" \
  --max_order 3 \
  --no_plots 2> "${OUTDIR}/fighi.time"

# Summarize FIGHI stats (pairs tested approximate from results)
python - <<PY
import pandas as pd, json, os
outdir = "${FIGHI_OUT}"
res = pd.read_csv(os.path.join(outdir,"fighi_results.csv"))
pair_n = (res['order']==2).sum() if 'order' in res.columns else 'NA'
print(f"[FIGHI] interactions total={len(res)}; pairs={pair_n}")
PY

echo "========== [4/5] PLINK --epistasis timing (same SNP subset) =========="
PLINK_DIR="${OUTDIR}/plink"
mkdir -p "${PLINK_DIR}"
/usr/bin/time -v python benchmarks/csv_to_plink.py \
  --csv "${FILTERED_CSV}" \
  --pheno "${PHENO}" \
  --out_prefix "${PLINK_DIR}/cd_top${TOP_M}" 2> "${OUTDIR}/plink_convert.time"

# Pairwise epistasis (all pairs among TOP_M SNPs)
# Use --fast-epistasis for speed (reports t-stat-like) or --epistasis (logistic)
# We'll use --fast-epistasis for a fair runtime baseline.
 /usr/bin/time -v plink \
  --bfile "${PLINK_DIR}/cd_top${TOP_M}" \
  --fast-epistasis \
  --epistasis-snp-set all \
  --allow-no-sex \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --out "${PLINK_DIR}/epi_fast" 2> "${OUTDIR}/plink.time"

echo "========== [5/5] MDR (pairwise) timing (same SNP subset) =========="
MDR_DIR="${OUTDIR}/mdr"
mkdir -p "${MDR_DIR}"
/usr/bin/time -v python benchmarks/mdr_benchmark.py \
  --csv "${FILTERED_CSV}" \
  --pheno "${PHENO}" \
  --out "${MDR_DIR}/mdr_pairs.csv" \
  --cv_folds 5 \
  --max_pairs 100000000 2> "${OUTDIR}/mdr.time"

echo "========== [ACCURACY] 80/20 held-out AUC =========="
python benchmarks/accuracy_eval.py \
  --csv "${FILTERED_CSV}" \
  --pheno "${PHENO}" \
  --fighi_results "${FIGHI_OUT}/fighi_results.csv" \
  --plink_epi "${PLINK_DIR}/epi_fast.epi.cc" \
  --mdr_pairs "${MDR_DIR}/mdr_pairs.csv" \
  --out_json "${OUTDIR}/accuracy.json" \
  --top_k_fighi 50

echo "========== [SUMMARY] =========="
echo "FIGHI timing:"
grep -E 'Elapsed|Maximum resident' "${OUTDIR}/fighi.time" || true
echo
echo "PLINK conversion timing:"
grep -E 'Elapsed|Maximum resident' "${OUTDIR}/plink_convert.time" || true
echo
echo "PLINK --epistasis timing:"
grep -E 'Elapsed|Maximum resident' "${OUTDIR}/plink.time" || true
echo
echo "MDR timing:"
grep -E 'Elapsed|Maximum resident' "${OUTDIR}/mdr.time" || true

echo
echo "[NOTE] Outputs:"
echo "  FIGHI:   ${FIGHI_OUT}/fighi_results.csv"
echo "  PLINK:   ${PLINK_DIR}/epi_fast.epi.cc"
echo "  MDR:     ${MDR_DIR}/mdr_pairs.csv"
