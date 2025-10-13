# FIGHI
FIGHI: Fisher-Information–Guided Hyper-interaction Inference for genome-wide epistasis discovery, summaries, and hypergraph outputs.

<p align="center">
  <img src="docs/IMAGES/fighi_cover.png" width="400" alt="FIGHI protocol sketch"/>
</p>

**FIGHI** is a fast, memory-aware method to discover **SNP–SNP (and higher-order) interactions** from genotype × phenotype data. It uses score-test approximations to **Fisher Information gain** to rank multi-locus combinations without fitting an enormous number of full models.

- **Goal:** find biologically meaningful epistasis and higher-order interactions
- **Fast & light:** vectorized score tests; memory-aware joins; chunked IO
- **Genomics-friendly IO:** TPED/TFAM/PLINK `.phen` → tidy CSV, rsIDs preserved
- **Network outputs:** Cytoscape/Gephi-ready hypergraphs
- **Annotation & enrichment:** optional SNP→gene mapping (cS2G/g:Profiler) & pathway/disease enrichment
- **HPC-friendly:** drop-in Slurm job (with disease tag), thread controls, no heavy dependencies

<p align="center">
  <img src="docs/IMAGES/fighi_protocol.png" width="600" alt="FIGHI protocol sketch"/>
</p>

---

## Why FIGHI?

Combinatorial interaction searches explode in cost. FIGHI sidesteps this by:
1. Screening variants (MAF-aware, standardized) to a feasible atom set.
2. Growing candidate tuples using an Apriori-style generator with variance pruning.
3. Ranking candidates via **score-test Fisher Information gain** (logistic/linear).
4. Adapting the maximum order `K` using information sufficiency and a power planner.
5. Emitting a **hypergraph** view and **SNP-level feature scores**, with optional annotation.

> TL;DR: You get interpretable interactions and network artifacts without OOMs.

---

## Install

### Option A: local (recommended while developing)
```bash
git clone https://github.com/<you>/FIGHI.git
cd FIGHI
# (optional) conda
conda env create -f environment.yml
conda activate fighi
# ensure parent of fighi_ext/ is visible
export PYTHONPATH="$(pwd):${PYTHONPATH:-}"
```


## Option B: editable install (lets you call python -m fighi_ext.run_cli)
```bash
pip install -e .
```

## Quick start

Minimal end-to-end on a merged CSV with phenotype column:

```bash
# ensure repo root on PYTHONPATH if running scripts in-place
export PYTHONPATH="$(pwd):${PYTHONPATH:-}"

python fighi_ext/run_cli.py \
  --csv path/to/merged.csv \
  --pheno case \
  --trait binary \
  --outdir fighi_ext/fighi_out \
  --max_order 4
```

From TPED/TFAM + PLINK .phen:

```bash
# 1) TPED → named numeric matrix (rsID columns, IID preserved)
python fighi_ext/tped_to_named_csv.py \
  --tped data/CD_origin.tped --tfam data/CD_origin.tfam \
  --phen_file data/CD_origin.phen \
  --out_csv data/CD_filtered_named.csv

# 2) Merge with phenotype (pandas-free, memory-lean)
python fighi_ext/merge_pheno_geno_nopandas.py \
  --geno_csv data/CD_filtered_named.csv \
  --phen_file data/CD_origin.phen \
  --out_csv data/CD_merged.csv \
  --phen_name case --impute_mean

# 3) Run FIGHI
python fighi_ext/run_cli.py \
  --csv data/CD_merged.csv --pheno case \
  --trait binary --outdir fighi_ext/fighi_out \
  --max_order 4
```

## Outputs (in fighi_ext/fighi_out/ by default):

- fighi_results.csv — interaction rows with order, fi_gain, beta_hat, info

- fighi_feature_scores.csv — per-SNP FI main/interaction totals (+ optional Gene/Pathway)

- fighi_hypergraph.gml, fighi_cytoscape.cyjs, fighi_hypergraph.hyper

- fighi_summary.json, fighi_model.pkl, fighi_log.txt

- plots/ — FI gain distribution, interaction heatmap, FI vs order, convergence traces


## Slurm (HPC)

A ready job is in fighi_ext/examples/demo_slurm.job. It supports a single disease tag used everywhere, saves outputs in fighi_ext/fighi_out, and exports thread envs. See docs/QUICKSTART.md for the complete script and notes.

## Method in 30 seconds

- Score-test FI gain (logistic/linear) estimates the information contributed by a new product feature (tuple of SNPs) without full refits.

- Apriori-style candidate growth + variance/MAF pruning keep the search sparse.

- Information ratio and planner stop at a sensible max order K under power constraints.

- Optional Westfall–Young controls family-wise error for permutations.

- Full derivations and references: docs/ALGORITHM.md.

## Cite

If you use FIGHI, please cite us (see also CITATION.cff):

```bash
@article{fighi2025,
  title={FIGHI: Fisher-Information–Guided Hyperinteraction Inference for Epistasis Discovery},
  author={Ariel G. Kemogne K.},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/xxxxx}
}
```

## License

MIT (see LICENSE).


## Contributing

Issues and PRs welcome! Please see CONTRIBUTING.md & CODE_OF_CONDUCT.md.

```yaml

---

## 2) `docs/QUICKSTART.md`

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

% SNP→Gene via cS2G or g:Profiler; add Gene/Pathway columns
python fighi_ext/annotate_fighi_features.py \
  --feature_csv fighi_ext/fighi_out/fighi_feature_scores.csv \
  --cs2g_dir fighi_ext/cS2G_1000GEUR \
  --gmt_files path/KEGG.gmt path/REACTOME.gmt path/GO_BP.gmt

### 4) Slurm example

See fighi_ext/examples/demo_slurm.job. One line to set disease:

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



For full derivations and equations, see [docs/THEORY.md](docs/THEORY.md)
