# Validation strategy

The automated suite covers:

- known Benjamini-Hochberg calculations;
- planted logistic and linear interactions;
- recovery of a planted genotype interaction;
- binary phenotype normalization and QC;
- candidate-limit safety behavior;
- output-directory protection;
- complete CLI demo and artifact parsing;
- GraphML well-formedness;
- TPED/TFAM conversion;
- local SNP-gene annotation and enrichment;
- reproducible truth-known binary and continuous simulations;
- PLINK 2 command construction, dosage-column normalization, IID-safe merging and hashes;
- generic multi-method execution, parsing, harmonized correction, truth recovery and overlap.

Before a scientific publication, add study-specific validation beyond unit tests:

1. Null type-I error across sample sizes, MAFs, case ratios, LD, and covariate structures.
2. Power across main-effect and interaction strengths.
3. Calibration under population structure and relatedness.
4. Comparison with full likelihood-ratio tests for manageable candidate sets.
5. Stability across random discovery/inference splits.
6. Runtime and memory benchmarks with versioned hardware information.
7. External biological replication.

The `simulate`, `benchmark-template`, and `benchmark-run` commands provide the executable foundation for items 1–6. The ARC simulation-array script covers scenario classes, but a manuscript should repeat multiple seeds and a pre-registered parameter grid rather than treating one seeded example as validation.

The repository's release verification uses a small seeded pairwise scenario to test the entire software path. It is a smoke/integration test, not an estimate of general power.

Software tests demonstrate implementation behavior; they do not establish universal scientific validity.
