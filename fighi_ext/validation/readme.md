# FIGHI – Theory Validation 

This folder contains **reproducible scripts** (headless-friendly) that validate the theoretical backbone of FIGHI.

## Contents
- `validate_score_vs_lrt.py` – Shows **rank equivalence** between FIGHI’s score-based Fisher Information gain and **likelihood-ratio test (LRT)** on simulated logistic data with true epistatic pairs.
- `validate_planner.py` – Validates the **K-planner detectability rule** (based on $N$, MAF, OR) against **empirical power** from simulations.
- `realdata_yeast_sanity.py` – Template to check recovery of **known epistasis pairs** on yeast (supply your dataset).

All scripts write outputs to `notebooks/out/`.

## How to run
```bash
cd notebooks
python validate_score_vs_lrt.py
python validate_planner.py
# optional (set paths inside first):
python realdata_yeast_sanity.py
```

## Expected outputs

### 1) Score vs LRT

```text
- score_vs_lrt_pairs.csv – per-pair metrics (ΔI, score-χ², LRT, p-values, ground truth flag)
- score_vs_lrt_scatter.png – scatter of ΔI vs LRT (with Spearman ρ on title)
- score_vs_lrt_qq.png – QQ plot of score χ² vs theoretical χ²(1)
- score_vs_lrt_summary.json – summary with rank correlations
```
```text
Interpretation: high Spearman ρ confirms ranking equivalence between FIGHI’s score-based ΔI and LRT when effects are modest, validating our use of the score test.
```

### 2) Planner (Detectability)

```text
- planner_grid.csv – rows across $(N,\text{MAF},\text{OR})$ with columns feasible_by_rule and empirical_power
- planner_heatmap.png – power heatmap for a fixed MAF slice
- planner_curves.png – empirical power vs $N$ for several OR
```
```text
Interpretation: the planner’s inequality correctly tracks the 80% power boundary; increases in $N$, MAF, and OR move the system into detectable regime.
```

### 3) Yeast sanity (optional)

```text
- yeast_pairs_score.csv – per known pair score and p-value
- yeast_overlap_stats.json – count of significant overlaps
```
```text
Interpretation: recovery of known pairs, or enrichment of low p-values, provides a quick real-data sanity check for the pipeline.
```

### Reproducibility

```text
- Deterministic RNG (NumPy Generator with fixed seeds)
- Pure Python + statsmodels, scipy, numpy, pandas, matplotlib (Agg backend)
- No internet / downloads required

### Requirements

```bash
pip install numpy pandas scipy statsmodels matplotlib
```
