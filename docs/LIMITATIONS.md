# Limitations and responsible interpretation

FIGHI is a focused interaction-testing engine, not a complete GWAS study-design system.

- Holdout inference reduces selection bias but uses fewer samples. Independent cohort replication remains stronger evidence.
- The score-test approximation can be inaccurate in very small samples, rare variants, separation, or extremely unbalanced phenotypes.
- Relatedness and population structure are not modeled as random effects. Include suitable principal components and use an appropriate study design; mixed-model extensions remain future work.
- LD creates correlated candidates and can complicate FDR interpretation and biological localization.
- Product coding captures a specific interaction parameterization. Dominant, recessive, haplotype, and nonmultiplicative encodings require explicit preprocessing or future extensions.
- Screening and hierarchical expansion can miss interactions without lower-order discovery evidence.
- Stability is conditional on a fixed finalist set and is not external replication.
- Pathway enrichment depends on the mapping, gene universe, database version, and multiple-testing assumptions.
- The software does not make clinical decisions and has not been cleared for diagnostic use.

Report null results, preprocessing, candidate universe, split seed, correction method, ancestry, covariates, and all deviations from defaults.

