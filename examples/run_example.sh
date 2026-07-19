#!/usr/bin/env bash
set -euo pipefail

fighi demo --outdir fighi_demo

# Replace the demo with a real analysis:
# fighi run cohort.csv \
#   --phenotype case \
#   --id-column IID \
#   --covariates age,sex,PC1,PC2 \
#   --trait binary \
#   --outdir results

