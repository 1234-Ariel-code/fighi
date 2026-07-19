#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="${FIGHI_REPO_DIR:-$PWD}"
ENV_PREFIX="${FIGHI_ENV_PREFIX:-$REPO_DIR/.conda/fighi-arc}"

cd "$REPO_DIR"

if command -v mamba >/dev/null 2>&1; then
  SOLVER=mamba
elif command -v conda >/dev/null 2>&1; then
  SOLVER=conda
else
  echo "Conda or Mamba is required. Load your cluster's Conda module first." >&2
  exit 2
fi

"$SOLVER" env create --prefix "$ENV_PREFIX" --file environment-arc.yml
"$SOLVER" run --prefix "$ENV_PREFIX" fighi --version
"$SOLVER" run --prefix "$ENV_PREFIX" python -m unittest discover -s tests -v

echo "FIGHI environment ready at $ENV_PREFIX"
echo "Activate with: conda activate $ENV_PREFIX"
