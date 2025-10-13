#!/usr/bin/env bash
set -euo pipefail
THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTDIR="${THIS_DIR}/toy_out_csv"
mkdir -p "${OUTDIR}"
python -m fighi_ext.run_cli   --csv "${THIS_DIR}/toy_merged.csv"   --pheno case   --trait binary   --outdir "${OUTDIR}"   --max_order 3
echo "Done. Outputs in: ${OUTDIR}"
