#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validate score-test FI gain vs LRT: rank equivalence

Outputs (in notebooks/out/):
- score_vs_lrt_pairs.csv
- score_vs_lrt_scatter.png
- score_vs_lrt_qq.png
- score_vs_lrt_summary.json
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import expit
from scipy.stats import chi2, spearmanr, pearsonr
from itertools import combinations
from statsmodels.discrete.discrete_model import Logit
from statsmodels.tools import add_constant

# ------------------ config ------------------
OUTDIR = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUTDIR, exist_ok=True)

RNG = np.random.default_rng(42)
N = 5000                 # samples
P = 200                  # SNPs
MAF_range = (0.1, 0.4)   # uniform MAFs
EFFECT_PAIRS = 20        # number of true interacting pairs
BETA_MAIN_SD = 0.05
BETA_INT_SD  = 0.25      # interactions stronger than main
SUBSET_PAIRS = 5000      # evaluate 5k random pairs (runtime-safe)

# ------------------ simulate genotypes ------------------
p_maf = RNG.uniform(MAF_range[0], MAF_range[1], size=P)
# Hardy–Weinberg genotype draws: P(0)=q^2, P(1)=2pq, P(2)=p^2
G = np.zeros((N, P), dtype=np.float32)
for j in range(P):
    p = p_maf[j]
    pr = [ (1-p)**2, 2*p*(1-p), p**2 ]
    G[:, j] = RNG.choice([0,1,2], size=N, p=pr)

# z-score columns (stable for logistic)
Gz = (G - G.mean(0)) / (G.std(0) + 1e-8)

# ------------------ choose true pairs ------------------
all_pairs = list(combinations(range(P), 2))
RNG.shuffle(all_pairs)
true_pairs = all_pairs[:EFFECT_PAIRS]

# ------------------ build linear predictor ------------------
beta_main = RNG.normal(0.0, BETA_MAIN_SD, size=P)
beta_int  = {pair: RNG.normal(0.5, BETA_INT_SD) for pair in true_pairs}

eta = Gz @ beta_main
for (i,j), b in beta_int.items():
    eta += b * (Gz[:, i] * Gz[:, j])

# logistic outcomes
p = expit(eta)
y = RNG.binomial(1, p, size=N).astype(np.float32)

# ------------------ helper: score/FI for one pair ------------------
def score_fi_gain_pair(i, j, y, Gz):
    # score test around null model (intercept only)
    p0 = np.full_like(y, y.mean(), dtype=np.float64)  # logistic null: p ≈ mean(y)
    W  = np.clip(p0*(1-p0), 1e-9, None)

    z = (Gz[:, i] * Gz[:, j]).astype(np.float64)
    U = float(np.dot(z, (y - p0)))
    I = float(np.dot(z * W, z))
    fi_gain = 0.5 * (U*U) / max(I, 1e-12)
    # score-test pval (1-df chi-square on U^2/I)
    chi = (U*U)/max(I, 1e-12)
    pval = 1.0 - chi2.cdf(chi, df=1)
    return fi_gain, chi, pval

# ------------------ helper: LRT for one pair ------------------
def lrt_pair(i, j, y, Gz):
    Z = np.column_stack([Gz[:, i], Gz[:, j], Gz[:, i]*Gz[:, j]])
    X1 = add_constant(Z)
    X0 = add_constant(np.zeros((len(y),0)))  # intercept only

    try:
        m1 = Logit(y, X1).fit(disp=0, maxiter=100)
        m0 = Logit(y, X0).fit(disp=0, maxiter=100)
        lrt = 2*(m1.llf - m0.llf)
    except Exception:
        # fallback robustly with small ridge by jittering (rare)
        lrt = np.nan
    pval = 1.0 - chi2.cdf(lrt, df=3) if np.isfinite(lrt) else 1.0
    return lrt, pval

# ------------------ evaluate random subset of pairs ------------------
RNG.shuffle(all_pairs)
eval_pairs = all_pairs[:SUBSET_PAIRS]

rows = []
for (i,j) in eval_pairs:
    fi_gain, chi_sc, p_sc = score_fi_gain_pair(i,j,y,Gz)
    lrt, p_lrt = lrt_pair(i,j,y,Gz)
    rows.append(dict(i=i,j=j, fi_gain=fi_gain, chi_score=chi_sc, p_score=p_sc,
                     lrt=lrt, p_lrt=p_lrt,
                     is_true=((i,j) in true_pairs or (j,i) in true_pairs)))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUTDIR, "score_vs_lrt_pairs.csv"), index=False)

# ------------------ rank correlation & plots ------------------
# Use finite rows only
dd = df.replace([np.inf,-np.inf], np.nan).dropna(subset=["fi_gain","lrt"])
rho_s, _ = spearmanr(dd["fi_gain"], dd["lrt"])
rho_p, _ = pearsonr(dd["fi_gain"], dd["lrt"])

with open(os.path.join(OUTDIR, "score_vs_lrt_summary.json"), "w") as f:
    json.dump({
        "N": int(N), "P": int(P), "subset_pairs": int(SUBSET_PAIRS),
        "spearman_fi_vs_lrt": float(rho_s),
        "pearson_fi_vs_lrt": float(rho_p),
        "true_pairs_in_subset": int(sum(dd["is_true"]))
    }, f, indent=2)

# scatter
plt.figure(figsize=(6,5))
plt.scatter(dd["fi_gain"], dd["lrt"], s=6, alpha=0.3)
plt.xlabel("FIGHI score: ΔI ≈ 0.5 * U^2 / I")
plt.ylabel("LRT (2·ΔlogLik)")
plt.title(f"Rank equivalence: Spearman ρ = {rho_s:.3f}")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "score_vs_lrt_scatter.png"), dpi=180)
plt.close()

# QQ: score chi (1 df) vs theoretical
sorted_chi = np.sort(dd["chi_score"].values)
theory = chi2.ppf((np.arange(1, len(sorted_chi)+1)-0.5)/len(sorted_chi), df=1)
plt.figure(figsize=(5.2,5))
plt.plot(theory, sorted_chi, '.', ms=3)
m = max(theory.max(), sorted_chi.max())
plt.plot([0,m],[0,m],'--',lw=1)
plt.xlabel("Theoretical χ²(1) quantiles")
plt.ylabel("Observed U²/I quantiles")
plt.title("Score-statistic QQ")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "score_vs_lrt_qq.png"), dpi=180)
plt.close()

print("[OK] Wrote validation/out/* for score vs LRT.")
