# Project audit and remediation record

The supplied prototype contained a promising high-level idea and useful prototype scripts, but it was not publishable software. This release preserves the Fisher-information-guided ranking concept while replacing unsafe or scientifically ambiguous parts.

| Prototype issue | Release remediation |
|---|---|
| Installation produced no `fighi` command | Standards-based package with console entry point |
| Main workflow crashed on a missing GraphML import | Rebuilt and tested exporters for every advertised format |
| GML, JSON, and Cytoscape files were returned as paths without being written | All artifact creation is validated end to end |
| Default threshold retained every positive FI gain | Formal score-test p-values and adjusted-p-value decisions |
| `pval` was always missing | Real asymptotic p-values with convergence diagnostics |
| `n_perm` and stability claims were not connected to the workflow | Removed false permutation claims; implemented clearly labeled conditional stability |
| Product terms were tested mostly against an intercept-only null | Covariates and all within-candidate lower-order terms are nuisance-adjusted |
| Screening and testing reused outcomes without warning | Independent discovery/inference holdout is the default |
| “Significant interactions” meant retained candidates | Significance now has one explicit adjusted-p-value definition |
| FI main-effect score column was always zero | Replaced by honest per-feature interaction contributions |
| Pickle output exposed portability and trust risks | Portable versioned JSON result manifest |
| Runtime, phenotype name, and version were hard-coded or missing | Accurate runtime, QC, trait, version, configuration, and provenance |
| Version values disagreed across package and report | Single v1.0.0 release version |
| Wide CSV “chunking” still concatenated everything | Removed misleading scalability promise; explicit limits and guidance |
| Batching claimed guaranteed cross-block interactions | Documentation now states coverage limitations |
| TPED conversion repeatedly rewrote the full CSV | Disk-backed two-pass conversion with structural validation |
| Annotation installed packages during execution | Offline local mapping/GMT flow; no runtime package installation |
| Large reference datasets were bundled in source | External resources are user-supplied and versionable |
| Documentation contained malformed fences/placeholders | Rewritten installation, methods, I/O, limitations, API, and release docs |
| No tests or functioning CI | Statistical, QC, CLI, conversion, annotation, export, and safety tests plus CI |

## Deliberate scope decisions

The release does not claim native PLINK BED/PGEN input, mixed-model relatedness correction, exhaustive genome-wide pair enumeration, distributed cross-block completeness, or clinical validity. These require separate validated engineering and methodological work and should not be implied by marketing text.

