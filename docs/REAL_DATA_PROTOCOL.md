# Real-data analysis protocol

This protocol turns an authorized PLINK cohort into a reproducible FIGHI analysis without placing genotype or phenotype data in the repository.

## 1. Freeze the design

Before looking at interaction results, record:

- primary phenotype and binary/continuous coding;
- inclusion/exclusion rules and one row per independent analysis unit;
- covariates, ancestry PCs, study site/batch terms and any matching variables;
- candidate-variant construction and external knowledge versions;
- discovery/inference split or an external replication cohort;
- primary interaction order, alpha and multiplicity correction;
- planned comparator methods and the capability differences that prevent direct pooling.

Use a versioned, outcome-independent `candidates.txt` whenever possible. Pathway-aware construction must record the pathway database and release. If outcome data influence candidate selection, preserve an untouched inference set.

## 2. Complete cohort QC upstream

FIGHI is not a full GWAS QC suite. Use established tooling to address sample call rate, sex checks where scientifically appropriate, duplicates/relatedness, heterozygosity anomalies, ancestry structure, variant call rate, MAF, Hardy-Weinberg filters appropriate to the design, allele alignment, imputation quality, batch effects and genome build. Preserve the exact commands and counts.

Create:

- a BED/BIM/FAM or PGEN/PVAR/PSAM prefix, or a VCF/BCF;
- `candidates.txt`, one stable variant ID per line;
- `keep.txt`, normally PLINK `FID IID` columns;
- an IID-keyed phenotype/covariate TSV with unique IDs.

Never commit controlled-access or identifiable data. Follow the cohort's data-use agreement and cluster storage policy.

## 3. Prepare the curated matrix

First inspect the command without executing PLINK:

```bash
fighi prepare-plink \
  --input /secure/cohort/qc/cohort \
  --input-type pgen \
  --candidates design/candidates.txt \
  --samples design/keep.txt \
  --phenotype-file /secure/cohort/phenotypes.tsv \
  --phenotype-column case \
  --covariates age,sex,PC1,PC2,PC3,PC4 \
  --threads 8 \
  --memory-mb 32000 \
  --outdir work/prepared_dry_run \
  --dry-run
```

Review `plink2_command.txt`, then run into a new directory without `--dry-run`. The command uses PLINK 2 additive export with `include-alt`, so prepared genotype columns count alternate alleles. Outputs include the compressed FIGHI table, exact retained candidate list, sample list, command, logs, missing-candidate list and SHA-256 provenance.

`prepare-plink` is designed for curated or already filtered candidate sets. A dense genome-wide sample-by-variant text matrix is usually inappropriate; use a scientifically justified reduction strategy and record it.

## 4. Validate before searching

```bash
fighi validate work/prepared/fighi_input.tsv.gz \
  --phenotype case \
  --id-column IID \
  --covariates age,sex,PC1,PC2,PC3,PC4 \
  --feature-file work/prepared/fighi_candidates.txt \
  --trait binary \
  --screen-method all \
  --max-order 2
```

Check sample count, class count, missingness, MAF filtering, constant features, encoded covariates and warnings. Reconcile every count with the upstream cohort flow diagram.

## 5. Run the primary analysis

```bash
fighi run work/prepared/fighi_input.tsv.gz \
  --phenotype case \
  --id-column IID \
  --covariates age,sex,PC1,PC2,PC3,PC4 \
  --feature-file work/prepared/fighi_candidates.txt \
  --trait binary \
  --screen-method all \
  --max-order 2 \
  --discovery-fraction 0.40 \
  --correction fdr_bh \
  --alpha 0.05 \
  --seed 2026 \
  --outdir results/primary
```

For a fixed, fully pre-specified pairwise universe, setting `--discovery-fraction 0` does not introduce outcome-driven pair selection inside FIGHI. It still changes the estimand's data reuse and must be declared. For adaptive screening or higher-order beam expansion, retain the default holdout or use an externally frozen candidate set.

## 6. Compare methods fairly

Generate and edit a manifest:

```bash
fighi benchmark-template \
  --output design/benchmark.json \
  --phenotype case \
  --trait binary \
  --covariates age,sex,PC1,PC2,PC3,PC4
```

Point `shared_inputs` to the common sample and variant files. Enable only installed comparators whose commands and output schemas have been verified on a small test. Read [Benchmarking](BENCHMARKING.md) before interpreting the harmonized tables.

## 7. Sensitivity and replication

At minimum, evaluate null calibration and power under MAF/LD/prevalence settings resembling the cohort, alternate discovery splits, population-structure sensitivity, genotype uncertainty policy, stricter multiplicity correction, top-result stability, and an independent replication dataset where permitted.

For replication, test the frozen discovery hyperedges in the replication cohort; do not rerun outcome-driven screening and call it replication. Align alleles/builds and preserve the original interaction direction convention.

## 8. Release an auditable analysis, not restricted data

Publish code, environment lock or container recipe, redacted design manifest, commands, hashes of non-public inputs where policy permits, aggregate flow counts, software logs, benchmark summaries and figures. Do not publish raw individual-level genotypes or phenotypes unless the governing consent and repository policy explicitly allow it.
