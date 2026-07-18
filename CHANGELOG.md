# Changelog

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

