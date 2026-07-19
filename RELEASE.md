# FIGHI 1.1.0 release summary

Release date: 2026-07-18

FIGHI 1.1.0 is the cluster and comparative-validation release of Fisher-Information-Guided Hyperinteraction Inference.

## Release contents

- Everything in the production-oriented 1.0 statistical engine and reporting package.
- Curated PLINK 2 BED/PGEN/VCF preparation with explicit alternate-allele dosage convention.
- Truth-known binary and continuous interaction simulations with reproducible provenance.
- A safe, versioned comparator runner for FIGHI, PLINK and generic MDR/GMDR, BOOST, BEAM, pathway-aware or information-theoretic outputs.
- Harmonized correction, truth recovery, ranking overlap, resource measurement, logs and HTML reports.
- ARC/Slurm jobs, an Apptainer recipe, a cluster Conda environment and a full real-data protocol.

## Verification target

- All 17+ automated tests pass.
- `scripts/verify_release.sh` completes the tests and a truth-known simulate-to-benchmark workflow.
- Source and wheel distributions build and pass package metadata checks.
- A clean installation reports version `1.1.0` and runs the CLI demo.

Exact verification results and artifact hashes are regenerated immediately before the release commit.

## Real-data status

No controlled-access real dataset is bundled, and this repository does not claim that a real-cohort comparison has already run. The ARC jobs and protocol are ready to submit once authorized genotype/phenotype paths, the allocation/partition and separately installed comparator executables are provided. Results must be reported as study-specific validation, not as a property of the package alone.

## Scientific status

The implementation is research software, not a clinical diagnostic product. A substantive publication still requires study-specific calibration, power, population-structure and relatedness sensitivity, comparator-specific assumptions, independent replication and transparent failed-method reporting.

## Release procedure

Follow `docs/RELEASING.md`. The release tag must be `v1.1.0` and match `pyproject.toml`, `src/fighi/__init__.py`, `CITATION.cff`, the container label and built artifacts.
