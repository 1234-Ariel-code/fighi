```markdown
# Quick Start

## 0) Environment
```bash
conda env create -f environment.yml
conda activate fighi
export PYTHONPATH="$(pwd):${PYTHONPATH:-}"
```

### 1) From merged CSV

python fighi_ext/run_cli.py \
  --csv path/to/merged.csv \
  --pheno case \
  --trait binary \
  --outdir fighi_ext/fighi_out \
  --max_order 4

### 2) From TPED/TFAM + .phen

python fighi_ext/tped_to_named_csv.py \
  --tped data/CD_origin.tped --tfam data/CD_origin.tfam \
  --phen_file data/CD_origin.phen \
  --out_csv data/CD_filtered_named.csv

python fighi_ext/merge_pheno_geno_nopandas.py \
  --geno_csv data/CD_filtered_named.csv \
  --phen_file data/CD_origin.phen \
  --out_csv data/CD_merged.csv \
  --phen_name case --impute_mean

python fighi_ext/run_cli.py \
  --csv data/CD_merged.csv --pheno case \
  --trait binary --outdir fighi_ext/fighi_out \
  --max_order 4

### 3) Annotation & enrichment (optional)

SNP→Gene via cS2G or g:Profiler; add Gene/Pathway columns
python fighi_ext/annotate_fighi_features.py \
  --feature_csv fighi_ext/fighi_out/fighi_feature_scores.csv \
  --cs2g_dir fighi_ext/cS2G_1000GEUR \
  --gmt_files path/KEGG.gmt path/REACTOME.gmt path/GO_BP.gmt

### 4) Slurm example

See fighi_ext/fighi_job.slurm One line to set disease:

DISEASE_NAME="CD"

Saves outputs to fighi_ext/fighi_out.


```yaml

---

## 3) `docs/USAGE.md`

```markdown
# Usage & Outputs

## CLI
```bash
python fighi_ext/run_cli.py \
  --csv <merged.csv> \
  --pheno <phen_col> \
  --trait {binary|linear} \
  --outdir <outdir> \
  --max_order 4 \
  [--target_or 1.3] [--read_chunksize 0]
```

### Important args

--max_order: hard cap on interaction order K (the search may stop earlier via information ratio).

--target_or: planner’s detectable odds ratio (guides feasible K under sample size & MAF).

--read_chunksize: if your merged CSV is huge, read in chunks (reduces peak RAM).

### Output files

#### fighi_results.csv

column	meaning

- hyperedge	tuple of rsIDs separated by `
order	interaction order K

- fi_gain	Fisher-Information gain (score-test)

- pval	optional, if permutation is enabled

- beta_hat	one-step estimate(s) for coefficients

- info	information (observed Fisher) contributing to cumulative ratio

#### fighi_feature_scores.csv

- Per-SNP totals aggregated across all retained edges:

- FI_total, FI_main, FI_interact, Rank, MAF

- Optional: Gene, Pathway (via annotate_fighi_features.py)

#### Hypergraphs

- fighi_hypergraph.gml — Gephi/Cytoscape import

- fighi_cytoscape.cyjs — direct Cytoscape session element JSON

- fighi_hypergraph.hyper — simple JSON with nodes + hyperedges

#### Summary & diagnostics

- fighi_summary.json — n_samples, n_snps, feasible K, runtime, etc.

- fighi_log.txt — run log with convergence notes & pruning summaries

- plots/ — FI distribution, interaction heatmap, FI vs order, convergence traces

#### Memory knobs

- Use merge_pheno_geno_nopandas.py to avoid pandas during the heaviest join.

- Use --read_chunksize in run_cli.py for very wide CSVs.

- Limit --max_order and/or lower the atom screening quota (see pipeline.py:screen).

- Control BLAS threads: OMP_NUM_THREADS, MKL_NUM_THREADS, OPENBLAS_NUM_THREADS.

```yaml

---

## 4) `docs/DATAFORMATS.md`

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
