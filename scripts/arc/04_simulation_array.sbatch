#!/usr/bin/env bash
#SBATCH --job-name=fighi-sim
#SBATCH --array=0-5
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=fighi-sim-%A_%a.log

set -euo pipefail

: "${FIGHI_ENV_PREFIX:?Set FIGHI_ENV_PREFIX}"
: "${FIGHI_SIMULATION_ROOT:?Set FIGHI_SIMULATION_ROOT}"

scenarios=(null main pairwise pairwise_main threeway structure)
scenario="${scenarios[${SLURM_ARRAY_TASK_ID:-0}]}"
seed="$(( ${FIGHI_BASE_SEED:-1700} + ${SLURM_ARRAY_TASK_ID:-0} ))"
run_dir="$FIGHI_SIMULATION_ROOT/${scenario}_seed${seed}"
trait="${FIGHI_SIM_TRAIT:-binary}"
phenotype="case"
if [[ "$trait" = "linear" ]]; then
  phenotype="trait"
fi
max_order="${FIGHI_SIM_MAX_ORDER:-2}"
if [[ "$scenario" = "threeway" && -z "${FIGHI_SIM_MAX_ORDER:-}" ]]; then
  max_order="3"
fi

conda run --prefix "$FIGHI_ENV_PREFIX" fighi simulate \
  --outdir "$run_dir" \
  --samples "${FIGHI_SIM_SAMPLES:-1000}" \
  --features "${FIGHI_SIM_FEATURES:-100}" \
  --trait "$trait" \
  --scenario "$scenario" \
  --effect "${FIGHI_SIM_EFFECT:-1.0}" \
  --ld-rho "${FIGHI_SIM_LD_RHO:-0.0}" \
  --seed "$seed"

conda run --prefix "$FIGHI_ENV_PREFIX" fighi benchmark-template \
  --output "$run_dir/benchmark.json" \
  --phenotype "$phenotype" \
  --trait "$trait" \
  --max-order "$max_order"

conda run --prefix "$FIGHI_ENV_PREFIX" fighi benchmark-run \
  "$run_dir/benchmark.json" \
  --outdir "$run_dir/benchmark_results"
