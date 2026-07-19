# ARC/Slurm quick start

The scripts in `scripts/arc/` use portable Slurm directives and environment variables. They deliberately do not hard-code an account, partition, module name or project-specific filesystem. Add your ARC allocation at submission time.

## Install

```bash
git clone https://github.com/1234-Ariel-code/fighi.git
cd fighi

# Load the Conda/Mamba module provided by your cluster, then:
export FIGHI_REPO_DIR="$PWD"
export FIGHI_ENV_PREFIX="$PWD/.conda/fighi-arc"
bash scripts/arc/00_setup_arc.sh
```

If the cluster does not permit Conda on compute nodes, build an Apptainer image from the login/build service:

```bash
apptainer build fighi-1.1.0.sif containers/Apptainer.def
apptainer run fighi-1.1.0.sif --version
```

PLINK 2 is included in `environment-arc.yml` for data preparation. Other comparators must be installed under their own licenses and referenced by absolute executable paths in the benchmark manifest.

## Prepare real data

Copy `configs/arc/real_data.example.env` outside the repository, fill in secure paths, and source it:

```bash
source /secure/project/fighi-real-data.env
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/arc/01_prepare_plink.sbatch
```

After reviewing the preparation manifest and validation output:

```bash
export FIGHI_INPUT="$FIGHI_PREPARED_DIR/fighi_input.tsv.gz"
export FIGHI_FEATURE_FILE="$FIGHI_PREPARED_DIR/fighi_candidates.txt"
export FIGHI_RESULTS_DIR="$FIGHI_PROJECT_ROOT/results/primary"
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/arc/02_run_fighi.sbatch
```

## Run the comparator manifest

```bash
export FIGHI_BENCHMARK_MANIFEST="$FIGHI_PROJECT_ROOT/design/benchmark.json"
export FIGHI_BENCHMARK_DIR="$FIGHI_PROJECT_ROOT/results/benchmark"
export FIGHI_BENCHMARK_STRICT=1
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/arc/03_run_benchmark.sbatch
```

`--strict` makes unavailable, timed-out, failed or unparseable enabled methods fail the Slurm job after all result artifacts have been written. Leave a method disabled while its installation or parser is being validated.

## Run the simulation matrix

```bash
export FIGHI_SIMULATION_ROOT="$FIGHI_PROJECT_ROOT/results/simulations"
export FIGHI_SIM_SAMPLES=2000
export FIGHI_SIM_FEATURES=200
export FIGHI_SIM_EFFECT=1.0
export FIGHI_SIM_LD_RHO=0.2
sbatch --account=YOUR_ACCOUNT --partition=YOUR_PARTITION \
  scripts/arc/04_simulation_array.sbatch
```

The six array tasks cover null, main-effect-only, pairwise, pairwise-plus-main, three-way and population-structure scenarios. For a publication, repeat seeds and parameter grids; record the submitted environment, Slurm job IDs, node type and allocated resources.

## Operational notes

- Keep the repository and environment on a project filesystem; keep temporary matrices and logs where the cluster policy permits.
- Export `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, and `OPENBLAS_NUM_THREADS` to avoid accidental oversubscription.
- Never write restricted cohort data to Git, a container image, or a public artifact directory.
- Output directories are protected against accidental overwrite. Use new, versioned directories for reruns.
- Benchmark commands execute serially inside one allocation so their resource measurements are not confounded by concurrent comparators. Use separate manifests/jobs when tools require materially different nodes.
