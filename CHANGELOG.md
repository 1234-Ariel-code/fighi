# Changelog

## 1.1.0 — 2026-07-18

Cluster and comparative-validation release.

### Added

- `prepare-plink` for curated BED/PGEN/VCF exports using PLINK 2, IID-safe phenotype merging, logs, allele convention and SHA-256 provenance.
- `simulate` for truth-known binary and continuous null, main-effect, pairwise, three-way, structure and LD scenarios.
- Manifest-driven `benchmark-validate` and `benchmark-run` commands with safe argument-array execution.
- FIGHI, PLINK epistasis and generic result parsers, including ranking-only comparator support without invented p-values.
- Harmonized multiple-testing correction, truth recovery, average precision, top-N overlap, runtime, CPU and memory reports.
- ARC/Slurm setup, preparation, analysis, benchmarking and simulation-array scripts.
- Apptainer recipe, ARC environment, real-data protocol and comparator capability documentation.
- Unit and integration coverage for simulations, PLINK preparation and multi-method benchmarking.

### Changed

- Automatic delimiter detection now recognizes compressed `.tsv.gz` and `.tab.gz` input.
- Cytoscape provenance now uses the package version dynamically.

### Scientific boundary

- Comparator outputs are harmonized only where quantities are validly available. Different null models, search strategies and phenotype/covariate support remain explicitly disclosed rather than treated as equivalent.

## 1.0.0 — 2026-07-17

First production-oriented release.

### Added

- Installable `fighi` package and multi-command CLI.
- Independent discovery/inference split for final statistical testing.
- Covariate- and hierarchy-adjusted efficient score tests.
- Global FDR or Bonferroni correction.
- Data validation, QC, categorical covariates, and explicit exploratory mode.
- Complete reports, network exports, provenance, and portable JSON manifests.
- Disk-backed TPED conversion and offline gene/pathway annotation.
- Demo generator, Python API, tests, CI, container, and release documentation.

### Changed

- Replaced ambiguous positive-FI retention with explicit adjusted-p-value significance.
- Replaced pickle result storage with versioned JSON.
- Removed bundled external annotation datasets and misleading batching guarantees.

### Fixed

- All advertised graph formats are now created.
- Version, phenotype, runtime, convergence, and output counts are now accurate.
- Interaction tests no longer confuse lower-order effects with the highest-order term.
