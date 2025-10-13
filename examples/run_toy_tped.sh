#!/usr/bin/env bash
set -euo pipefail
THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_DIR="$( cd "${THIS_DIR}/.." && pwd )"
OUTDIR="${THIS_DIR}/toy_out_tped"
mkdir -p "${OUTDIR}"

python "${REPO_DIR}/fighi_ext/tped_to_named_csv.py"   --tped "${THIS_DIR}/toy.tped"   --tfam "${THIS_DIR}/toy.tfam"   --phen_file "${THIS_DIR}/toy.phen"   --out_csv "${THIS_DIR}/toy_named.csv"   --ref_mode major

python "${REPO_DIR}/fighi_ext/merge_pheno_geno_nopandas.py"   --geno_csv "${THIS_DIR}/toy_named.csv"   --phen_file "${THIS_DIR}/toy.phen"   --out_csv "${THIS_DIR}/toy_merged.csv"   --phen_name case   --impute_mean

python -m fighi_ext.run_cli   --csv "${THIS_DIR}/toy_merged.csv"   --pheno case   --trait binary   --outdir "${OUTDIR}"   --max_order 3
echo "Done. Outputs in: ${OUTDIR}"
