# Python API

## Configuration

`AnalysisConfig` is a dataclass containing every analysis decision. Call `validate()` to check a configuration and `to_dict()` for serialization.

## Model

`FIGHI(config).fit(features, phenotype, covariates=None)` accepts pandas objects and returns `SearchResult`. Genotype validation is enabled by default; pass `allow_non_genotype=True` only for intentional generic numeric analyses.

`FIGHI(config).run(prepared_data)` accepts the validated object returned by `prepare_dataframe` or `prepare_file`.

## Result

`SearchResult.to_frame()` returns every inference-stage candidate. `SearchResult.significant` returns interaction objects passing the configured global threshold.

The public API follows semantic versioning. Fields may be added in minor versions; removals or changed meanings require a major version.

