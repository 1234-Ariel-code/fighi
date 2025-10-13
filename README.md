# FIGHI
FIGHI: Fisher-Information–Guided Hyper-interaction Inference for genome-wide epistasis discovery, summaries, and hypergraph outputs.

**FIGHI** is a fast, memory-aware method to discover **SNP–SNP (and higher-order) interactions** from genotype × phenotype data. It uses score-test approximations to **Fisher Information gain** to rank multi-locus combinations without fitting an enormous number of full models.

- 🔎 **Goal:** find biologically meaningful epistasis and higher-order interactions
- ⚡ **Fast & light:** vectorized score tests; memory-aware joins; chunked IO
- 🧬 **Genomics-friendly IO:** TPED/TFAM/PLINK `.phen` → tidy CSV, rsIDs preserved
- 🕸️ **Network outputs:** Cytoscape/Gephi-ready hypergraphs
- 🧪 **Annotation & enrichment:** optional SNP→gene mapping (cS2G/g:Profiler) & pathway/disease enrichment
- 🖥️ **HPC-friendly:** drop-in Slurm job (with disease tag), thread controls, no heavy dependencies

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
