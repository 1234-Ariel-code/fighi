```yaml

---

```markdown
# Data Formats

## Genotype input
- **CSV**: rows=samples, columns=rsIDs (+1 phenotype column). Values can be 0/1/2 encodings.
- **TPED/TFAM**: use `tped_to_named_csv.py` to produce a named matrix CSV; rsIDs become column names.

## Phenotype
- PLINK `.phen` with columns: `FID IID PHENO`. Supported encodings:
  - Binary {1,2} → mapped to {0,1}
  - Already {0,1}
  - Missing coded `-9` is dropped
- Name of the phenotype column in the merged CSV is given by `--phen_name` (default `case`).

## Merging
`merge_pheno_geno_nopandas.py`:
- If genotype CSV lacks an `IID` column but has the **same number of rows** as `.phen`, IIDs are injected from `.phen`.
- Inner-joins on `IID`, coerces numerics safely, optional mean imputation (`--impute_mean`).
```
