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
- local SNP-gene annotation and enrichment.

Before a scientific publication, add study-specific validation beyond unit tests:

1. Null type-I error across sample sizes, MAFs, case ratios, LD, and covariate structures.
2. Power across main-effect and interaction strengths.
3. Calibration under population structure and relatedness.
4. Comparison with full likelihood-ratio tests for manageable candidate sets.
5. Stability across random discovery/inference splits.
6. Runtime and memory benchmarks with versioned hardware information.
7. External biological replication.

Software tests demonstrate implementation behavior; they do not establish universal scientific validity.

