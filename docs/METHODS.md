# Methods

## Objective

FIGHI tests whether an order-`K` product term adds phenotype information after conditioning on covariates and every lower-order term formed from the same candidate features.

For candidate features `x1, ..., xK`, the tested feature is:

`z = x1 × ... × xK`.

The nuisance design contains an intercept, user covariates, all main effects, and all products of orders 2 through `K-1`. Genotype columns are standardized before products are constructed.

## Efficient score statistic

Under the fitted null model, let `U = z'(y - μ)` be the score. The information for `z` is residualized against the nuisance design so that correlated lower-order terms do not contribute to the highest-order test. The one-degree-of-freedom statistic is:

`T = U² / I_effective`.

Under regularity conditions and the null hypothesis, `T` is asymptotically chi-squared with one degree of freedom. FIGHI reports:

- `statistic = T`
- `fi_gain = T / 2`
- `beta_hat = U / I_effective`, a one-step effect estimate
- `p_value = P(ChiSq(1) >= T)`

Binary outcomes use logistic nuisance models fitted with damped Newton updates. Continuous outcomes use ordinary least squares and residual variance scaling.

## Independent discovery and inference

The default stratified holdout divides samples into discovery and inference partitions using the configured seed.

1. Discovery partition: quality-controlled features are screened and higher-order candidates are expanded.
2. Inference partition: the frozen candidate set is scored independently.
3. All inference p-values across all evaluated orders are adjusted together.

This design avoids testing on the same outcomes used to select candidates. `discovery_fraction=0` disables splitting and is labeled exploratory in every output.

## Screening

- `marginal`: absolute association with the covariate-residualized phenotype.
- `variance`: highest feature variance after QC.
- `hybrid`: 75% marginal and 25% high-variance features.
- `all`: no outcome-driven feature screen.

Hybrid screening diversifies the selected set but cannot guarantee recovery of pure epistasis. A biologically curated feature list is often preferable.

## Candidate expansion

Pairs are enumerated from screened atoms. For order 3 and above, Apriori joins require every order-`K-1` subset to be present among the discovery beam. The beam contains statistically strongest discovery candidates up to `beam_width`.

This is a computational search rule, not a biological law. It encodes hierarchical heredity and may miss interactions whose lower-order subsets have weak discovery evidence.

## Multiple testing

Final inference p-values are adjusted globally across every evaluated candidate and order using one of:

- Benjamini-Hochberg FDR (`fdr_bh`, default)
- Bonferroni family-wise error control (`bonferroni`)
- no adjustment (`none`, exploratory)

Dependence among variants can affect FDR assumptions. Bonferroni is more conservative and is the safer option when strict family-wise control is required.

## Stability

Optional stability repeats subsample only the inference partition and re-score the strongest final candidates. The stability value is the fraction of repeats in which a candidate passes multiplicity correction within that fixed candidate subset. It is a conditional robustness diagnostic, not a replacement for replication.

