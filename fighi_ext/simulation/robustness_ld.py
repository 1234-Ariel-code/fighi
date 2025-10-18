#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Robustness to LD (correlated SNPs) and model misspecification.

Scenarios:
  A) No LD (baseline)
  B) Moderate LD within blocks (rho=0.6)
  C) Strong LD (rho=0.85)
  D) Misspecified: outcome uses a nonlinear link (probit), we test with logistic score

Outputs:
- robust_summary.csv
- robust_bar.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import norm
from common_sim import make_rng, simulate_genotypes, zscore_cols, inject_ld, score_fi_gain_vector, bh_fdr

OUTDIR = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUTDIR, exist_ok=True)

RNG = make_rng(11)
REPS = 200
N = 6000
P = 400
M_SCAN = 3000
OR_TRUE = 1.4
MAF = 0.2

def power_scenario(kind):
    hits = 0
    for _ in range(REPS):
        G, _ = simulate_genotypes(N, P, maf_low=MAF, maf_high=MAF, rng=RNG)
        Gz = zscore_cols(G)

        # inject LD
        if kind == "LD-moderate":
            Gz = inject_ld(Gz, blocks=8, rho=0.6, rng=RNG)
        elif kind == "LD-strong":
            Gz = inject_ld(Gz, blocks=8, rho=0.85, rng=RNG)

        # choose true pair
        i, j = np.random.default_rng().choice(P, size=2, replace=False)
        zprod = Gz[:, i]*Gz[:, j]
        beta = np.log(OR_TRUE)

        # outcome generation
        if kind == "Misspec-probit":
            # latent normal -> probit-like link
            eta = beta * zprod
            y = (norm.cdf(eta) > np.random.rand(len(eta))).astype(np.float32)
        else:
            from scipy.special import expit
            p = expit(beta * zprod)
            y = np.random.binomial(1, p, size=len(p)).astype(np.float32)

        # scan pairs (including true)
        all_pairs = list(combinations(range(P), 2))
        np.random.default_rng().shuffle(all_pairs)
        pairs = set(all_pairs[:M_SCAN-1])
        pairs.add(tuple(sorted((i,j))))
        pairs = list(pairs)

        pvals = []
        true_index = None
        for k,(a,b) in enumerate(pairs):
            z = (Gz[:, a] * Gz[:, b])
            _, _, p = score_fi_gain_vector(y, z)
            pvals.append(p)
            if (a,b)==tuple(sorted((i,j))):
                true_index = k
        q = bh_fdr(pvals)
        if q[true_index] < 0.05:
            hits += 1
    return hits/REPS

scenarios = ["No-LD", "LD-moderate", "LD-strong", "Misspec-probit"]
rows = []
for s in scenarios:
    pow_est = power_scenario(s)
    rows.append(dict(scenario=s, power=pow_est))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUTDIR, "robust_summary.csv"), index=False)

plt.figure(figsize=(6,4))
plt.bar(df["scenario"], df["power"])
plt.ylim(0,1)
plt.ylabel("Power (BH<0.05)")
plt.title("Robustness across LD / misspecification")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "robust_bar.png"), dpi=180)
plt.close()

print("[OK] Robustness results written to simulation/out/")
