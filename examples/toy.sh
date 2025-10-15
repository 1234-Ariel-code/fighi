#!/usr/bin/env bash
# FIGHI_FULL (node-run)
# Former SBATCH settings for reference:
# --cpus-per-task=22, --mem=300G, --time=24:00:00, --partition=gpu-a100,theia, --gres=gpu:1
set -euo pipefail

# --------------------------
# 0. Logging
# --------------------------
mkdir -p logs
LOGSTAMP="$(date +%Y%m%d_%H%M%S)"
LOGFILE="logs/fighi_${LOGSTAMP}.out"
exec > >(tee -a "${LOGFILE}") 2>&1
echo "[INFO] Node-run started at ${LOGSTAMP}"

# --------------------------
# 1. Disease setup
# --------------------------
DISEASE_NAME="toy"         # set once
PHENO_NAME="case"

# Threads: use SLURM_CPUS_PER_TASK if present, else nproc, else 8
THREADS="${SLURM_CPUS_PER_TASK:-$(command -v nproc >/dev/null 2>&1 && nproc || echo 8)}"
export OMP_NUM_THREADS="${THREADS}"
export MKL_NUM_THREADS="${THREADS}"
export OPENBLAS_NUM_THREADS="${THREADS}"
export NUMEXPR_NUM_THREADS="${THREADS}"

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

# Try modules if available (safe no-ops on systems without environment modules)
if command -v module >/dev/null 2>&1; then
  module purge || true
  module load python/3.10 || module load python || true
fi

if [[ "${USE_CONDA}" -eq 1 ]]; then
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV}" || conda create -y -n "${CONDA_ENV}" python=3.10
    conda activate "${CONDA_ENV}"
  else
    echo "[WARN] conda not found; proceeding with system Python."
  fi
fi

# Ensure required Python packages are present
python - <<'PY'
import importlib, sys, subprocess
pkgs = ["numpy","pandas","scipy","matplotlib","statsmodels"]
for pkg in pkgs:
    try:
        importlib.import_module(pkg)
    except Exception:
        subprocess.check_call([sys.executable,"-m","pip","install","-q",pkg])
PY

# Add current folder (fighi_ext) to sys.path so imports work
export PYTHONPATH="${CODE_DIR}:${PYTHONPATH:-}"

echo "[INFO] Running FIGHI inside: ${CODE_DIR}"
echo "[INFO] Using data from:     ${DATA_DIR}"
echo "[INFO] Outputs to:          ${OUTDIR}"
echo "[INFO] Threads:             ${THREADS}"

# --------------------------
# 4. Pipeline
# --------------------------
# == Step 0: TPED -> named CSV ==
# (Uncomment if you want to regenerate from TPED/TFAM.)
#: <<'STEP0'
#echo "== Step 0: TPED -> named CSV =="
#if [[ -s "${GENO_NAMED}" ]]; then
#  echo "[INFO] Found existing ${GENO_NAMED}; skipping conversion."
#else
#  if [[ -f "${TPED}" ]]; then
#    if [[ -f "${TFAM}" ]]; then
#      python tped_to_named_csv.py \
#        --tped "${TPED}" --tfam "${TFAM}" \
#        --phen_file "${PHEN_FILE}" \
#        --out_csv "${GENO_NAMED}" --ref_mode major
#    else
#      python tped_to_named_csv.py \
#        --tped "${TPED}" --phen_file "${PHEN_FILE}" \
#        --out_csv "${GENO_NAMED}" --ref_mode major
#    fi
#  else
#    echo "[ERROR] Missing TPED: ${TPED}"; exit 2
#  fi
#fi
#STEP0

# == Step 1: Merge (no-pandas) ==
# (Uncomment if you want to recreate the merged file.)
#: <<'STEP1'
#echo "== Step 1: Merge (no-pandas) =="
#python merge_pheno_geno_nopandas.py \
#  --geno_csv "${GENO_NAMED}" \
#  --phen_file "${PHEN_FILE}" \
#  --out_csv "${MERGED_CSV}" \
#  --phen_name "${PHENO_NAME}" \
#  --impute_mean
#STEP1

# == Step 1.5: Fix-up merged CSV to include IID + phenotype (case) and proper order ==
# This block guarantees:
#   - Column "IID" exists (prefers TFAM IID; else synthetic I1..IN)
#   - Column "case" exists (prefers PHEN file; else renames first non-IID col)
#   - Column order: IID, case, rs*
TMP_MERGED="${MERGED_CSV%.csv}_with_iid.csv"
python - <<PY
import pandas as pd, sys
csv_path = r"""${MERGED_CSV}"""
tfam     = r"""${TFAM}"""
phen     = r"""${PHEN_FILE}"""

df = pd.read_csv(csv_path)

# Attach/ensure IID
if "IID" not in df.columns:
    iids = None
    try:
        tf = pd.read_csv(tfam, delim_whitespace=True, header=None)
        # TFAM columns: FID, IID, father, mother, sex, pheno
        if len(tf) == len(df):
            iids = tf.iloc[:,1].astype(str).tolist()
    except Exception:
        pass
    if iids is None:
        iids = [f"I{i}" for i in range(1, len(df)+1)]
    df.insert(0, "IID", iids)

# Ensure 'case' phenotype
if "case" not in df.columns:
    try:
        ph = pd.read_csv(phen, delim_whitespace=True, header=None, names=["FID","IID","PHENO"])
        if len(ph) == len(df):
            ph["IID"] = ph["IID"].astype(str)
            # unify types
            df["IID"] = df["IID"].astype(str)
            df = df.merge(ph[["IID","PHENO"]], on="IID", how="left")
            df.rename(columns={"PHENO":"case"}, inplace=True)
        else:
            raise ValueError("PHEN length mismatch")
    except Exception:
        # fallback: rename the first non-IID column to 'case'
        cand = [c for c in df.columns if c != "IID"][0]
        df.rename(columns={cand:"case"}, inplace=True)

# Reorder: IID, case, then the rest (keep existing SNP order)
front = ["IID","case"]
rest = [c for c in df.columns if c not in front]
df = df[front + rest]

out = r"""${TMP_MERGED}"""
df.to_csv(out, index=False)
print(f"[INFO] Using merged CSV with IID header -> {out}")
PY
MERGED_CSV="${TMP_MERGED}"

# == Step 2: Run FIGHI ==
echo "== Step 2: Run FIGHI =="
python run_cli.py \
  --csv    "${MERGED_CSV}" \
  --pheno  "${PHENO_NAME}" \
  --trait  binary \
  --outdir "${OUTDIR}" \
  --max_order 4 \
  --prescreen_top_m 200 \
  --col_block 1000 \
  --row_chunksize 20000

# (Optional) Step 3: annotate SNPs
# CS2G_DIR="${CODE_DIR}/cS2G_1000GEUR"
# echo "== Step 3: Annotate SNPs (optional) =="
# python annotate_fighi_features.py \
#   --feature_csv "${OUTDIR}/fighi_feature_scores.csv" \
#   --cs2g_dir "${CS2G_DIR}" \
#   --out_csv "${OUTDIR}/fighi_feature_scores_annotated.csv"

echo "== Done. Outputs in: ${OUTDIR} =="
echo "[INFO] Log saved to: ${LOGFILE}"

