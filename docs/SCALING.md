# Scaling and HPC guidance

The search size grows combinatorially. With `M` screened features, pairwise evaluation requires `M(M-1)/2` candidate tests. Order 3 can be much larger without beam pruning.

## Recommended progression

1. Validate the input.
2. Run the demo to verify the environment.
3. Start with `max_order=2`, `top_m=30-50`, and no stability repeats.
4. Inspect QC, convergence, runtime, and candidate counts.
5. Increase `top_m` or order only with a sample-size and compute justification.

The candidate safety limit stops a run before an accidental combinatorial expansion. Raising it is an explicit decision.

## Genome-wide data

Outcome-based top-M screening is pragmatic but may miss pure interactions and changes the scientific target. For genome-wide analyses, consider an externally defined candidate set based on prior biology, LD pruning, functional annotation, or an independent cohort.

Independent chromosome or batch runs cannot guarantee detection of cross-batch interactions. Anchors and overlaps improve coverage but do not create a mathematical guarantee. Any distributed batching strategy must document which cross-block pairs were never evaluated.

## Slurm

`scripts/fighi_slurm.sh` is a conservative single-analysis template. Set its variables, copy it alongside the project, and submit with `sbatch`. The script fails on unset variables and records the exact FIGHI outputs.

## Threading

FIGHI itself is deterministic and currently evaluates candidates serially. Numerical libraries may use BLAS threads. On shared clusters, set `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, and `OPENBLAS_NUM_THREADS` to match allocated CPUs.

