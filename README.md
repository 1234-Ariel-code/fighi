# FIGHI
FIGHI: Fisher-Information–Guided Hyper-interaction Inference for genome-wide epistasis discovery, summaries, and hypergraph outputs.

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

```bash
fighi/
├─ fighi_ext/                     # Python package + scripts
│  ├─ __init__.py
│  ├─ adaptive.py
│  ├─ apriori.py
│  ├─ fisher.py
│  ├─ glm.py
│  ├─ hypergraph.py
│  ├─ knockoff_perm.py
│  ├─ pipeline.py
│  ├─ report.py
│  ├─ run_cli.py                  # can run as `python -m fighi_ext.run_cli` or `python fighi_ext/run_cli.py`
│  ├─ tped_to_named_csv.py
│  ├─ merge_pheno_geno_nopandas.py
│  ├─ annotate_fighi_features.py
│  └─ utils.py
├─ examples/
│  └─ (tiny toy data / minimal commands)
├─ docs/
│  ├─ README_theory.md            # full math/theory (we drafted earlier)
│  └─ assets/
│     ├─ fighi_protocol.png       # flowchart
│     └─ fighi_cover.png          # cover art
├─ fighi_job.slurm                # SLURM job (disease-var enabled)
├─ README.md                      # user-facing quickstart
├─ LICENSE                        # MIT (if chosen)
└─ .gitignore
```

## Install

### Option A: local (recommended while developing)
```bash
git clone https://github.com/1234-Ariel-code/FIGHI.git
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

## Quick start (see [docs/QUICKSTART.md](docs/QUICKSTART.md))

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

A ready job is in fighi_ext/fighi_job.slurm It supports a single disease tag used everywhere, saves outputs in fighi_ext/fighi_out, and exports thread envs. See [docs/ALGORITHM.md](docs/ALGORITHM.md) for the complete script and notes.

## Method in 30 seconds

- Score-test FI gain (logistic/linear) estimates the information contributed by a new product feature (tuple of SNPs) without full refits.

- Apriori-style candidate growth + variance/MAF pruning keep the search sparse.

- Information ratio and planner stop at a sensible max order K under power constraints.

- Optional Westfall–Young controls family-wise error for permutations.

- Full derivations and references: [docs/THEORY.md](docs/THEORY.md).

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

Issues and PRs welcome! Please see [docs/CONTRIBUTING.md](docs/CONTRIBUTING.md) & [docs/CODE_OF_CONDUCT.md](docs/CODE_OF_CONDUCT.md).


For full derivations and equations, see [docs/THEORY.md](docs/THEORY.md)



<p align="center">
  <img src="docs/IMAGES/fighi_cover.png" width="200" alt="FIGHI protocol sketch"/>
</p>

