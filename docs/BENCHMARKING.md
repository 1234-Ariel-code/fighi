# Reproducible comparison framework

FIGHI 1.1 includes a manifest-driven benchmark runner. It executes each enabled method as an argument array without a shell, captures logs and software versions, measures each run in a fresh resource-monitor process, parses results into a common interaction table, reapplies one declared multiple-testing correction to raw p-values, and produces truth-recovery or cross-method overlap summaries. On Linux/ARC the monitor records wall time, child CPU time and peak resident memory; unsupported platforms retain wall time and mark unavailable fields.

Third-party executables and restricted datasets are not bundled. Install and cite each comparator separately, then declare its exact command, output schema and statistical target in the JSON manifest.

## Two valid benchmark tracks

Do not mix these tracks in one headline result.

### Fixed-universe search benchmark

Give every method the same samples and the same pre-specified variants. For pairwise FIGHI runs, use `--screen-method all --discovery-fraction 0`; the all-pairs universe is fixed before observing the outcome, so no outcome-driven pair selection occurs inside FIGHI. This track compares ranking, calibration, recovery, time and memory. If a method adaptively searches higher orders on the same samples, label those results exploratory.

### Protected-inference benchmark

Freeze candidate construction using prior knowledge or a discovery partition, then run final tests on the same untouched inference samples for every method. A method that cannot accept the common inference sample set or the required covariates belongs in a separate capability stratum. This is the preferred track for inferential claims.

## Start a truth-known benchmark

```bash
fighi simulate \
  --outdir work/pairwise_seed17 \
  --samples 1000 \
  --features 100 \
  --scenario pairwise \
  --effect 1.0 \
  --ld-rho 0.2 \
  --seed 17

fighi benchmark-template \
  --output work/pairwise_seed17/benchmark.json

fighi benchmark-validate work/pairwise_seed17/benchmark.json
fighi benchmark-run work/pairwise_seed17/benchmark.json \
  --outdir work/pairwise_seed17/comparison \
  --strict
```

Simulation scenarios are `null`, `main`, `pairwise`, `pairwise_main`, `threeway`, and `structure`; both binary and continuous traits are supported. The generated `truth.json` records configuration, causal hyperedges, MAFs, file hashes and sample counts.

## Manifest contract

The required top-level values are `analysis_id`, `methods`, and valid method objects. `shared_inputs` paths are resolved relative to the manifest and must exist during validation.

Each enabled method declares:

- `name`: filesystem-safe unique identifier;
- `command`: argument array, never a shell string;
- `version_command`: optional short version query;
- `cwd` and `env`: optional execution context;
- `timeout_seconds`: optional safety limit;
- `target` and `notes`: scientific capability disclosure;
- `result.path`: output file, with placeholders allowed;
- `result.format`: `fighi`, `plink_epistasis`, or `generic`.

Available placeholders are `{manifest_dir}`, `{outdir}`, `{run_dir}`, `{python_executable}` (the interpreter running the benchmark), and every key declared in `shared_inputs` such as `{candidate_file}`.

The generic parser accepts `feature_columns`, `p_column`, `effect_column`, `separator`, and `comment`. Ranking-only methods may set `p_column` to `null` and provide `score_column`; set `score_ascending` to `true` only when smaller scores rank higher. Ranking scores are used for top-N recovery and overlap but never converted into invented p-values or q-values.

## Comparator coverage

The framework supports, but does not redistribute, the following comparator families.

| Family | Typical scope | Adapter | Required disclosure |
|---|---|---|---|
| FIGHI | Binary/continuous; pairwise to order 6; covariates; all within-candidate lower terms | Native | Split mode and candidate construction |
| PLINK 1.9 | Widely used pairwise epistasis workflows | `plink_epistasis` | Exact command, phenotype mode and PLINK statistic |
| MDR/GMDR | Non-parametric combinations and cross-validation ranking | `generic` | Cross-validation, permutation and covariate procedure |
| BOOST | Fast pairwise case-control screening | `generic` | Binary phenotype coding, filtering and second-stage test |
| BEAM/Bayesian mapping | Bayesian case-control interaction modeling | `generic` | Priors, MCMC/search settings and reported posterior quantity |
| WTCCC-derived exhaustive tools | Usually method- and dataset-specific exhaustive pair scans | `generic` | Executable/version, test statistic and filters |
| Pathway-aware methods | Candidate construction informed by sets, networks or pathways | `generic` | Knowledge source version and whether outcome information entered selection |
| Information-theoretic tools | Entropy, mutual information or information gain ranking | `generic` | Estimator, discretization, bias correction and significance procedure |

This matrix is not a claim that the methods are interchangeable. Check the documentation for the exact installed version. PLINK's official association documentation describes `--epistasis` and `--fast-epistasis`, and its data-export documentation defines the additive export used by `fighi prepare-plink`: [PLINK 1.9 association](https://www.cog-genomics.org/plink/1.9/assoc) and [PLINK 2 data export](https://www.cog-genomics.org/plink/2.0/data). A command-line GMDR implementation is described in the [GMDR paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5320543/).

## PLINK example

The generated template contains a disabled PLINK 1.9 adapter. Replace `PLINK_PREFIX`, verify case/control coding and the output suffix for the installed version, then enable it. Use the same `candidate_file` and `sample_file` hashes recorded by the runner.

For quantitative traits, do not reuse the case-control result suffix. Configure the command and `result.path` for the actual PLINK output and document the target in `notes`.

## Outputs

| Output | Purpose |
|---|---|
| `validated_inputs.json` | Resolved shared inputs and SHA-256 hashes |
| `benchmark_runs.csv` | Status, exact argv, version, runtime, CPU, memory and logs |
| `benchmark_interactions.csv` | Normalized hyperedges, raw p-values/scores and harmonized q-values |
| `benchmark_method_summary.csv` | Test counts, discoveries, truth recovery and average precision |
| `benchmark_overlap.csv` | Pairwise top-N Jaccard overlap |
| `benchmark_summary.json` | Machine-readable design and status summary |
| `benchmark_report.html` | Human-readable comparison report |

## What a publishable comparison must report

Report dataset accession and release, consent/access constraints, sample and variant QC, ancestry/relatedness policy, phenotype coding, covariates, allele coding, candidate-universe construction, train/test boundaries, every software version and command, failed or unavailable methods, hardware allocation, seeds, raw tested counts, method-native correction, harmonized correction, and whether each result is confirmatory or exploratory.

Do not interpret faster runtime, greater result count, or cross-method agreement as evidence of biological truth. Calibration, recovery under truth-known simulations, independent replication, and method-specific assumptions remain essential.
