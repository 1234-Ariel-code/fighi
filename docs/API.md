# Python API

## Configuration

`AnalysisConfig` is a dataclass containing every analysis decision. Call `validate()` to check a configuration and `to_dict()` for serialization.

## Model

`FIGHI(config).fit(features, phenotype, covariates=None)` accepts pandas objects and returns `SearchResult`. Genotype validation is enabled by default; pass `allow_non_genotype=True` only for intentional generic numeric analyses.

`FIGHI(config).run(prepared_data)` accepts the validated object returned by `prepare_dataframe` or `prepare_file`.

## Result

`SearchResult.to_frame()` returns every inference-stage candidate. `SearchResult.significant` returns interaction objects passing the configured global threshold.

The public API follows semantic versioning. Fields may be added in minor versions; removals or changed meanings require a major version.

## Preparation, simulation and benchmarking

- `fighi.plink.build_plink2_export_command(...)` validates PLINK inputs and returns a shell-free argument array.
- `fighi.plink.prepare_plink_dataset(...)` runs the curated export and writes a FIGHI-ready table with provenance.
- `fighi.simulation.SimulationConfig` and `simulate_dataset(...)` create truth-known in-memory benchmarks.
- `fighi.simulation.write_simulation(...)` writes data, candidates, samples and truth metadata.
- `fighi.benchmark.validate_benchmark_manifest(...)` resolves and hashes shared inputs without running methods.
- `fighi.benchmark.run_benchmark(...)` executes enabled adapters and writes normalized comparisons.
- `fighi.benchmark.write_benchmark_template(...)` creates an editable JSON manifest.

These modules are public in 1.1 but intentionally remain explicit imports; the compact top-level namespace continues to expose the statistical model classes.
