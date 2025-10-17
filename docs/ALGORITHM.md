# Algorithm (FIGHI)

We seek multi-locus interactions (hyperedges) that increase Fisher Information with respect to the phenotype model.


## Inputs

```text
   - Genotype–phenotype CSV file with:

       - Phenotype column y

       - Genotype columns X₁ … Xₚ

       - Optional ID column

   - Trait type ∈ {binary, linear}

   - Maximum interaction order K_max

   - Target effect size proxy (e.g. odds ratio OR₀)

   - Prescreen parameters: (M, row_chunksize, col_block)

   - Pruning thresholds: (ε_var, ε_maf, τ_gain)

   - Optional permutation/stability settings
```

## Outputs

```text
    - Hypergraph H = (V, E) where:

         - nodes V = SNPs

         - edges E = discovered interactions

         - weights 𝑤𝑒=ΔI(e)

    - Per-SNP Fisher Information (FI) summaries

    - Tables, logs, and plots
```

## Algorithm 1 — Main FIGHI Procedure

```text
FIGHI(csv, y_name, trait, K_max, OR₀, prescreen, prune):
```
<p align="center">
  <img src="docs/IMAGES/alg1.png" width="600" alt="FIGHI algo sketch"/>
</p>

<p align="center">
  <img src="docs/IMAGES/alg2.png" width="600" alt="FIGHI algo sketch"/>
</p>

<p align="center">
  <img src="docs/IMAGES/alg3.png" width="600" alt="FIGHI algo sketch"/>
</p>

## Algorithm 2 — Prescreen: Streaming Correlation-based SNP Filtering

```text
Prescreen(csv, y_name, M, row_chunksize, col_block):
```
<p align="center">
  <img src="docs/IMAGES/alg3.png" width="600" alt="FIGHI algo sketch"/>
</p>

## Algorithm 3 — Score-Test Statistics and Fisher Information Gain

```text
ScoreTestStats(z, X, y, trait):
```
<p align="center">
  <img src="docs/IMAGES/alg3.png" width="600" alt="FIGHI algo sketch"/>
</p>

## Algorithm 4 — Apriori-style Candidate Expansion with Pruning

```text
ExpandCandidates(E, K):
```
<p align="center">
  <img src="docs/IMAGES/alg3.png" width="600" alt="FIGHI algo sketch"/>
</p>


## Algorithm 5 — Practical Planner for K

<p align="center">
  <img src="docs/IMAGES/alg3.png" width="600" alt="FIGHI algo sketch"/>
</p>

<p align="center">
  <img src="docs/IMAGES/alg3.png" width="600" alt="FIGHI algo sketch"/>
</p>












## 1. Screening
Standardize each SNP; compute a simple correlation (or point-biserial). Keep top-M SNPs passing MAF threshold.

## 2. Candidate growth
- For K=2, enumerate all pairs from atoms.
- For K>2, use an Apriori join of frequent (K−1) tuples, prune by subset presence and variance of the product feature.

## 3. Score-test Fisher Information gain
For binary traits (logistic):
$$\Delta I = \tfrac{1}{2} \beta_{\text{1step}}^2 \cdot I_{xx}$$
where $\beta_{\text{1step}} = U / I_{xx}$, $U = x^\top(y - p)$, $I_{xx} = x^\top W x$, $W = p(1-p)$, $p = \sigma(X\beta)$.  
For linear traits (OLS) analogous forms hold with $\sigma^2$.

## 4. Adaptive K
- `planner_max_K`: uses sample size N, MAF range and target OR to estimate feasible maximum K.
- `information_ratio`: stops early if cumulative information at K−1 is ≥ `eps_ratio` (default 0.95).

## 5. Multiple testing (optional)
Westfall–Young max-T permutation (family-wise error control).

## 6. Artifacts
- Hypergraph exports: GML and Cytoscape JSON nodes + hyperedge nodes.
- SNP-level aggregation builds FI totals and ranks.

See in-code references in `fighi_ext/fisher.py`, `fighi_ext/adaptive.py`, `fighi_ext/apriori.py`.
