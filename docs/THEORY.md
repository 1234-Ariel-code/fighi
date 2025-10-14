
# Theoretical Foundations of FIGHI  
**Fisher-Information–Guided Hyperinteraction Inference**

---

## 1. Motivation

Genome-wide interaction analysis (epistasis) explores joint effects of multiple loci on a phenotype.  
Traditional regression or likelihood-ratio tests become computationally prohibitive for high-order interactions because the number of candidate tuples grows combinatorially with the number of SNPs.

**FIGHI** re-formulates the search as an *incremental Fisher Information Gain (FI-gain)* problem:

> “Which new combination of variables most increases the expected Fisher Information of the model with respect to the phenotype?”

This reframing enables interaction discovery without refitting full models at each step.

---

## 2. Model Setting

Let:
- $X \in \mathbb{R}^{n \times p}$ be the standardized genotype matrix  
  (rows = individuals, columns = SNPs encoded as 0, 1, 2)
- $y \in \{0,1\}^n$ (binary trait) or $y \in \mathbb{R}^n$ (continuous trait)
- $f_\beta(x)$ be the model mean (logistic or linear)

### 2.1 Logistic Trait

$$
\Pr(y_i = 1 \mid x_i) = \sigma(x_i^\top \beta)
= \frac{1}{1 + e^{-x_i^\top \beta}}.
$$

### 2.2 Linear Trait

$$
y_i = x_i^\top \beta + \varepsilon_i, \quad \varepsilon_i \sim \mathcal{N}(0, \sigma^2).
$$

---

## 3. Fisher Information (Score-Test Formulation)

Let $\ell(\beta)$ be the log-likelihood and $U(\beta) = \frac{\partial \ell}{\partial \beta}$ the score.

For a small perturbation in a new feature $z = x_{j_1} * x_{j_2} * \dots * x_{j_K}$,  
the one-step update under the **score test** is:

$$
\hat{\beta}_{\text{new}}^{(1)} = I^{-1}_{zz} U_z,
$$

where  
$I_{zz} = \mathbb{E}\!\left[-\frac{\partial^2\ell}{\partial\beta_z^2}\right]$ and  
$U_z = \frac{\partial\ell}{\partial\beta_z}$.

The **Fisher Information Gain** for adding this feature is approximated as:

$$
\Delta \mathcal{I}(z) = \tfrac{1}{2} (\hat{\beta}_{\text{new}}^{(1)})^2 I_{zz}
= \tfrac{1}{2} \frac{U_z^2}{I_{zz}}.
$$

---

## 4. Computing $U_z$ and $I_{zz}$

### Binary Trait (Logistic)

Let $p = \sigma(X\beta)$, $W = \mathrm{diag}(p(1-p))$:

$$
U_z = z^\top (y - p), \qquad I_{zz} = z^\top W z.
$$

Hence:

$$
\Delta \mathcal{I}(z)
= \tfrac{1}{2}\frac{(z^\top(y-p))^2}{z^\top W z}.
$$

### Linear Trait

$$
U_z = z^\top (y - X\beta), \qquad
I_{zz} = \frac{z^\top z}{\sigma^2}, \quad
\Delta \mathcal{I}(z)
= \tfrac{1}{2}\frac{(z^\top(y - X\beta))^2}{z^\top z}.
$$

These require only vector operations—no full model refits.

---

## 5. Interaction Generation

### 5.1 Apriori Growth

For $K$-way interactions:

$$
z_{(j_1,\dots,j_K)} = \prod_{k=1}^K x_{j_k}.
$$

Candidates at order $K+1$ are formed only from frequent or high-information subsets at order $K$.

### 5.2 Pruning Criteria

- **Subset pruning:** Drop any tuple whose all $(K−1)$ subsets are not retained.  
- **Variance pruning:** Discard tuples with low variance $\mathrm{Var}(z) < \epsilon_v$.  
- **MAF filter:** Remove SNPs with minor allele frequency below threshold $\epsilon_{maf}$.

---

## 6. Adaptive Order Selection

FIGHI automatically determines a practical upper limit $K_{\max}$ via an *information ratio* rule:

$$
r_K = \dfrac{\sum_{k'=1}^{K} \Delta \mathcal{I}_{k'}}{\sum_{k'=1}^{K_{\max}^{\text{theor}}} \Delta \mathcal{I}_{k'}}
$$

Stop increasing $K$ if $r_K > \tau$ (default 0.95).  

A theoretical upper bound $K_{\max}^{\mathrm{theor}}$ can be estimated from sample size $N$, minor allele frequency (MAF), and a target detectable odds ratio $\mathrm{OR}_0$:

$$
K_{\max}^{\mathrm{theor}} \approx
\max \{ K : N \cdot \mathrm{MAF}_{\min}^K \cdot (\log \mathrm{OR}_0)^2 \ge z_{1-\alpha/2}^2 \}.
$$

This is implemented in `adaptive.py:planner_max_K`.

---

## 7. Aggregation: SNP-Level Feature Information

After exploring all edges,  
for each SNP $s$ appearing in interactions $\mathcal{E}(s)$:

$$
\mathrm{FI}_{\text{main}}(s)
 = \sum_{e \in \mathcal{E}(s),\,|e|=1} \Delta \mathcal{I}(e),
\qquad
\mathrm{FI}_{\text{interact}}(s)
 = \sum_{e \in \mathcal{E}(s),\,|e|>1} \Delta \mathcal{I}(e)
$$

and the total contribution:

$$
\mathrm{FI}_{\text{total}}(s)
 = \mathrm{FI}_{\text{main}}(s) + \mathrm{FI}_{\text{interact}}(s)
$$

Ranks are computed by descending $\mathrm{FI}_{\text{total}}$.

---

## 8. Connection to Classical Tests

| Method | Statistic | Needs Model Refit? | Notes |
|---------|------------|--------------------|-------|
| Likelihood-ratio | $2(\ell_1-\ell_0)$ | ✅ Yes | Full logistic regression |
| Wald | $\hat{\beta}^2 / \mathrm{Var}(\hat{\beta})$ | ✅ Yes |  |
| Score (FIGHI) | $U_z^2 / I_{zz}$ | ❌ No | One-step Fisher-information test |

Thus, **FIGHI** provides the *score-test analogue* of the LRT, yielding comparable ranking under small-effect assumptions but at much lower computational cost.

---

## 9. Output as Hypergraph

Each retained tuple $e = \{s_1,\dots,s_K\}$ becomes a **hyperedge** in $H=(V,E)$, where $V$ are SNPs.  
Weights correspond to FI gain.

Adjacency tensors or incidence matrices can be constructed as:

$$
A_{ij} = \sum_{e\ni i,j} w_e, \quad w_e = \Delta\mathcal{I}(e).
$$

**Exports:**
- `.gml` for Gephi/NetworkX  
- `.cyjs` for Cytoscape  
- JSON hyperedge list  

---

## 10. Multiple Testing Correction (Optional)

Under the permutation framework (Westfall–Young):

1. Shuffle the phenotype y^(π)  
2. Re-run the FI pipeline  
3. Record the maximum FI per order K  
4. Estimate the empirical p-value as:  
   p_e = Probability(FI_π ≥ FI_obs)

---

## 11. Complexity Analysis

Let:  
- $N$: samples  
- $P$: retained SNPs  
- $K$: maximum interaction order  
- $M_K$: number of surviving tuples at order $K$

Rough complexity:

$$
O\!\left(\sum_{K=1}^{K_{\max}} M_K N\right),
\qquad
M_K \ll \binom{P}{K} \text{ due to pruning.}
$$

Memory: $O(N)$ if each feature is computed on-the-fly (chunked vector product).

---

## 12. Interpretation

FI-gain correlates with potential **predictive stability** and **causal relevance**:

- High $\Delta \mathcal{I}$ — strong evidence that a combination explains phenotype variance beyond marginals.  
- Comparing $FI_{main}$ vs. $FI_{interact}$ distinguishes additive vs. epistatic signal.  
- Aggregated FI profiles can feed downstream enrichment or polygenic risk estimation.

---

## 13. References

Kemogne, A.G. et al. (2025). *FIGHI: Fisher-Information–Guided Hyperinteraction Inference for Epistasis Discovery.*

---

*This document serves as the mathematical appendix to the FIGHI repository.*
