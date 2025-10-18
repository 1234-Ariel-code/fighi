#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validate the planner's detectability inequality against empirical power.

We vary N, MAF, OR (odds-ratio for a K-way interaction term) and simulate
binary outcomes under a logistic model with ONLY a K-way product term.

Outputs:
- planner_grid.csv: grid with feasible_by_rule, empirical_power
- planner_heatmap.png: heatmap of power vs (N, OR) at fixed MAF and K
- planner_curves.png: curves of power across N for multiple OR
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import expit
from scipy.stats import norm
from itertools import product

OUTDIR = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUTDIR, exist_ok=True)

RNG = np.random.default_rng(7)

# ------------------ grids ------------------
K = 2                 # test pairwise detectability; you can change to 3
N_grid = [1000, 2000, 4000, 8000, 16000]
OR_grid = [1.1, 1.2, 1.3, 1.5, 1.7, 2.0]
MAF_grid = [0.1, 0.2, 0.3, 0.4]   # we’ll loop, and plot one slice
alpha = 0.05
power_target = 0.8

B = 200   # permutations/replicates for empirical power

# ------------------ planner inequality ------------------
def planner_feasible(N, maf, K, OR0, alpha=0.05, power=0.8):
    z = norm.ppf(1 - alpha/2.0) + norm.ppf(power)
    beta = np.log(max(OR0, 1.0))
    if beta <= 0: 
        return False
    I_req = (z/beta)**2
    # crude per-tuple information ~ (p*(1-p))^K
    p = maf
    per_tuple_I = (p*(1-p))**K
    return (N * per_tuple_I) >= I_req

# ------------------ simulate power ------------------
def simulate_power(N, maf, K, OR0, B=200):
    # simulate K SNPs, others irrelevant
    q = 1 - maf
    pr = [q*q, 2*maf*q, maf*maf]
    # K product term
    beta = np.log(OR0)
    hits = 0
    for b in range(B):
        X = []
        for j in range(K):
            X.append(RNG.choice([0,1,2], size=N, p=pr))
        X = np.column_stack(X).astype(np.float32)
        Xz = (X - X.mean(0)) / (X.std(0)+1e-8)
        zprod = np.prod(Xz, axis=1)
        # logistic y with only interaction effect
        eta = beta * zprod
        p = expit(eta)
        y = RNG.binomial(1, p, size=N)

        # score test (null: intercept only)
        p0 = np.full_like(y, y.mean(), dtype=np.float64)
        W  = np.clip(p0*(1-p0), 1e-9, None)
        z  = zprod.astype(np.float64)
        U = float(np.dot(z, (y - p0)))
        I = float(np.dot(z*W, z))
        chi = (U*U)/max(I,1e-12)
        # reject at alpha?
        from scipy.stats import chi2
        if 1.0 - chi2.cdf(chi, df=1) < alpha:
            hits += 1
    return hits / B

rows = []
for maf in MAF_grid:
    for N in N_grid:
        for OR0 in OR_grid:
            feasible = planner_feasible(N, maf, K, OR0, alpha=alpha, power=power_target)
            pow_emp = simulate_power(N, maf, K, OR0, B=B)
            rows.append(dict(K=K, N=N, MAF=maf, OR=OR0,
                             feasible_by_rule=feasible,
                             empirical_power=pow_emp))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUTDIR, "planner_grid.csv"), index=False)

# --------------- Plot a heatmap at fixed MAF ---------------
maf_show = 0.2
dd = df[df["MAF"]==maf_show].pivot(index="N", columns="OR", values="empirical_power")
plt.figure(figsize=(7,4.5))
im = plt.imshow(dd.values, aspect="auto", origin="lower")
plt.xticks(range(len(dd.columns)), [f"{c:.2f}" for c in dd.columns], rotation=0)
plt.yticks(range(len(dd.index)), [str(i) for i in dd.index])
plt.colorbar(im, label="Empirical power")
plt.title(f"Empirical power heatmap (K={K}, MAF={maf_show})")
plt.xlabel("Odds ratio (OR)")
plt.ylabel("N")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "planner_heatmap.png"), dpi=180)
plt.close()

# --------------- Curves across N for multiple OR ---------------
plt.figure(figsize=(7,4.2))
for OR0 in OR_grid:
    yy = [df[(df["OR"]==OR0) & (df["MAF"]==maf_show) & (df["K"]==K) & (df["N"]==n)]["empirical_power"].iloc[0]
          for n in N_grid]
    plt.plot(N_grid, yy, marker='o', label=f"OR={OR0:.2f}")
plt.axhline(0.8, ls="--", lw=1, label="0.80 target")
plt.xlabel("N")
plt.ylabel("Empirical power")
plt.title(f"Power vs N at MAF={maf_show}, K={K}")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "planner_curves.png"), dpi=180)
plt.close()

print("[OK] Wrote notebooks/out/* for planner validation.")
