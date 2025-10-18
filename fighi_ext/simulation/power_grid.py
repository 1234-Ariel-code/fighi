#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Power for detecting a known interacting pair with BH correction
among M scanned pairs (includes distractors).

Outputs:
- power_grid.csv
- power_heatmap.png
- power_curves.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import combinations
from common_sim import make_rng, simulate_genotypes, zscore_cols, score_fi_gain_vector, bh_fdr

OUTDIR = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUTDIR, exist_ok=True)

RNG = make_rng(7)
REPS = 200

N_grid = [1000, 2000, 4000, 8000, 16000]
OR_grid = [1.1, 1.2, 1.3, 1.5, 2.0]
MAF_grid= [0.1, 0.2, 0.3]

P = 400           # total SNPs
M_SCAN = 3000     # scanned pairs total

def one_power(N, maf, OR):
    # simulate
    G, _ = simulate_genotypes(N, P, maf_low=maf, maf_high=maf, rng=RNG)
    Gz   = zscore_cols(G)
    # plant a true pair
    i, j = RNG.choice(P, size=2, replace=False)
    beta = np.log(OR)
    eta  = beta * (Gz[:, i] * Gz[:, j])
    from scipy.special import expit
    p    = expit(eta)
    y    = RNG.binomial(1, p, size=N).astype(np.float32)

    # choose pairs to scan (ensure true pair included)
    all_pairs = list(combinations(range(P), 2))
    RNG.shuffle(all_pairs)
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
    pvals = np.asarray(pvals)
    q = bh_fdr(pvals)
    hit = (q[true_index] < 0.05)
    return hit

rows = []
for maf in MAF_grid:
    for N in N_grid:
        for OR in OR_grid:
            hits = 0
            for _ in range(REPS):
                hits += one_power(N, maf, OR)
            pow_est = hits/REPS
            rows.append(dict(N=N, MAF=maf, OR=OR, power=pow_est))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUTDIR, "power_grid.csv"), index=False)

# heatmap at MAF=0.2
maf_show = 0.2
dd = df[df["MAF"]==maf_show].pivot(index="N", columns="OR", values="power")
plt.figure(figsize=(7,4.5))
im = plt.imshow(dd.values, aspect="auto", origin="lower")
plt.xticks(range(len(dd.columns)), [f"{c:.2f}" for c in dd.columns])
plt.yticks(range(len(dd.index)), [str(i) for i in dd.index])
plt.colorbar(im, label="Power (BH<0.05)")
plt.title(f"Power heatmap (MAF={maf_show})")
plt.xlabel("OR")
plt.ylabel("N")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "power_heatmap.png"), dpi=180)
plt.close()

# curves across N for each OR (MAF=0.2)
plt.figure(figsize=(7,4.2))
for OR in OR_grid:
    yy = [df[(df["MAF"]==maf_show)&(df["OR"]==OR)&(df["N"]==n)]["power"].iloc[0] for n in N_grid]
    plt.plot(N_grid, yy, marker='o', label=f"OR={OR:.2f}")
plt.axhline(0.8, ls='--', lw=1, label='0.80')
plt.xlabel("N")
plt.ylabel("Power")
plt.title(f"Power vs N (MAF={maf_show})")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "power_curves.png"), dpi=180)
plt.close()

print("[OK] Power grid written to simulation/out/")
