#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Real-data sanity template: yeast epistasis pairs.

Expected inputs:
- yeast.csv with columns: IID, phenotype (binary or continuous), SNP_... columns
- known_pairs.csv with columns: snp1, snp2 (matching yeast.csv column names)

Outputs (notebooks/out/):
- yeast_overlap_stats.json
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2

OUTDIR = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUTDIR, exist_ok=True)

YEAST_CSV   = "/path/to/yeast.csv"         # TODO: set
PHENO_NAME  = "pheno"                       # TODO: set
KNOWN_PAIRS = "/path/to/known_pairs.csv"   # TODO: set

if not (os.path.exists(YEAST_CSV) and os.path.exists(KNOWN_PAIRS)):
    raise FileNotFoundError("Please set YEAST_CSV and KNOWN_PAIRS to real files.")

df = pd.read_csv(YEAST_CSV)
y  = df[PHENO_NAME].to_numpy()
X  = df.drop(columns=[PHENO_NAME, "IID"], errors="ignore")
Xz = (X - X.mean(0)) / (X.std(0) + 1e-8)

pairs = pd.read_csv(KNOWN_PAIRS)
pairs = [tuple(row) for row in pairs[["snp1","snp2"]].itertuples(index=False, name=None)]

def score_pair(c1, c2):
    z = (Xz[c1].values * Xz[c2].values).astype(np.float64)
    if set(np.unique(y)) <= {0,1}:
        p0 = np.full_like(y, y.mean(), dtype=np.float64)
        W  = np.clip(p0*(1-p0), 1e-9, None)
        U = float(np.dot(z, (y - p0)))
        I = float(np.dot(z*W, z))
    else:
        r = y - y.mean()
        U = float(np.dot(z, r))
        I = float(np.dot(z, z) / (np.var(r) + 1e-9))
    chi = (U*U)/max(I, 1e-12)
    pval = 1.0 - chi2.cdf(chi, df=1)
    return chi, pval

obs_stats = []
hits = 0
for s1, s2 in pairs:
    if s1 not in Xz.columns or s2 not in Xz.columns:
        continue
    chi, pval = score_pair(s1, s2)
    obs_stats.append(dict(snp1=s1, snp2=s2, chi=chi, pval=pval))
    if pval < 0.05 / max(1,len(pairs)):  # Bonferroni for quick sanity
        hits += 1

pd.DataFrame(obs_stats).to_csv(os.path.join(OUTDIR, "yeast_pairs_score.csv"), index=False)
with open(os.path.join(OUTDIR, "yeast_overlap_stats.json"), "w") as f:
    json.dump({"tested_pairs": len(obs_stats), "significant_hits_bonf": hits}, f, indent=2)

print("[OK] Yeast sanity finished. See validate/out/")
