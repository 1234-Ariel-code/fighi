#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scaling experiment: how runtime grows with P (number of SNPs) when scanning M pairs.

Outputs:
- scaling_pairs.csv
- scaling_pairs.png
"""

import os, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import combinations
from common_sim import make_rng, simulate_genotypes, zscore_cols, score_fi_gain_vector

OUTDIR = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUTDIR, exist_ok=True)

RNG = make_rng(9)
N = 4000
P_grid = [200, 400, 800, 1200, 1600]
M_SCAN = 4000  # fixed scanned pairs

rows = []
for P in P_grid:
    G, _ = simulate_genotypes(N, P, rng=RNG)
    Gz = zscore_cols(G)
    y  = (RNG.random(N) < 0.5).astype(np.float32)

    all_pairs = list(combinations(range(P), 2))
    RNG.shuffle(all_pairs)
    pairs = all_pairs[:M_SCAN]

    t0 = time.time()
    for (i,j) in pairs:
        _ = score_fi_gain_vector(y, Gz[:, i]*Gz[:, j])
    t1 = time.time()
    rows.append(dict(P=P, pairs_tested=M_SCAN, seconds=t1-t0))

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUTDIR, "scaling_pairs.csv"), index=False)

plt.figure(figsize=(6,4))
plt.plot(df["P"], df["seconds"], marker='o')
plt.xlabel("P (SNPs)")
plt.ylabel("Seconds (scan M pairs)")
plt.title("FIGHI score-scan scaling (pairwise)")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "scaling_pairs.png"), dpi=180)
plt.close()

print("[OK] Scaling results written to simulation/out/")
