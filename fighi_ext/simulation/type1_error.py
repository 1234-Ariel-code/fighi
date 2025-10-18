#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Type I error under the null (no genetic effects).
Measures empirical rejection rates at nominal alpha after BH correction
when scanning M random pairs.

Outputs (sim/out/):
- t1e_summary.csv
- t1e_reject_rates.png
- t1e_qq.png
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

# ------------- settings -------------
RNG = make_rng(202)
REPS = 200               # repetitions
N = 4000                 # samples
P = 500                  # SNPs
M_PAIRS = 2000           # random tested pairs/repetition
ALPHAS = [0.10, 0.05, 0.01]

def one_rep():
    G, _ = simulate_genotypes(N, P, rng=RNG)
    Gz   = zscore_cols(G)
    # pure null: Bernoulli with base rate ~ 0.5 (slight jitter)
    y = (RNG.random(N) < 0.5).astype(np.float32)

    # random pairs
    all_pairs = list(combinations(range(P), 2))
    RNG.shuffle(all_pairs)
    pairs = all_pairs[:M_PAIRS]

    pvals = []
    for (i,j) in pairs:
        z = (Gz[:, i] * Gz[:, j])
        _, _, p = score_fi_gain_vector(y, z)
        pvals.append(p)
    pvals = np.asarray(pvals)
    q = bh_fdr(pvals)
    return pvals, q

rows = []
qq_all = []
for r in range(REPS):
    p, q = one_rep()
    qq_all.extend(p.tolist())
    for a in ALPHAS:
        rej = float((q < a).mean())
        rows.append(dict(rep=r, alpha=a, reject_rate=rej))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUTDIR, "t1e_summary.csv"), index=False)

# reject-rate plot
plt.figure(figsize=(6,4))
for a in ALPHAS:
    dd = df[df["alpha"]==a]["reject_rate"].values
    plt.plot(np.arange(1, REPS+1), dd, alpha=0.4, label=f"α={a}")
plt.axhline(0.10, ls='--', lw=1, color='C0')
plt.axhline(0.05, ls='--', lw=1, color='C1')
plt.axhline(0.01, ls='--', lw=1, color='C2')
plt.xlabel("Repetition")
plt.ylabel("Empirical rejection rate (BH)")
plt.title("Type I error under null (pairwise scan)")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "t1e_reject_rates.png"), dpi=180)
plt.close()

# QQ plot for raw p-values
qq = np.sort(np.asarray(qq_all))
u  = (np.arange(1, len(qq)+1) - 0.5)/len(qq)
plt.figure(figsize=(4.8,4.8))
plt.plot(-np.log10(u), -np.log10(qq), '.', ms=2)
m = max((-np.log10(u)).max(), (-np.log10(qq)).max())
plt.plot([0,m],[0,m],'--',lw=1)
plt.xlabel("Expected -log10(p)")
plt.ylabel("Observed -log10(p)")
plt.title("Null calibration QQ")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "t1e_qq.png"), dpi=180)
plt.close()

print("[OK] Type I error results written to simulation/out/")
