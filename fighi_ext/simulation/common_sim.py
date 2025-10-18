#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common simulation utilities for FIGHI simulation studies.
"""

import os
import time
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import chi2

# ---------- RNG ----------
def make_rng(seed=123):
    return np.random.default_rng(seed)

# ---------- Genotype simulation ----------
def simulate_genotypes(N, P, maf_low=0.1, maf_high=0.4, rng=None):
    """
    Hardy–Weinberg draws per SNP with MAF~Uniform(maf_low, maf_high).
    Returns:
      G  (N,P) int8  in {0,1,2}
      p_maf (P,)     MAF for each SNP
    """
    rng = make_rng() if rng is None else rng
    p_maf = rng.uniform(maf_low, maf_high, size=P)
    G = np.empty((N, P), dtype=np.int8)
    for j, p in enumerate(p_maf):
        q = 1 - p
        pr = [q*q, 2*p*q, p*p]
        G[:, j] = rng.choice([0,1,2], size=N, p=pr)
    return G, p_maf

def zscore_cols(M):
    M = M.astype(np.float32)
    mu = M.mean(axis=0)
    sd = M.std(axis=0) + 1e-8
    return (M - mu)/sd

# ---------- Outcome simulation ----------
def simulate_logistic_y(Gz, beta_main=None, beta_pairs=None, beta_triples=None, rng=None):
    """
    Logistic outcome with optional main/pair/triple terms.
      beta_main: (P,) or None
      beta_pairs: dict[(i,j)] -> beta
      beta_triples: dict[(i,j,k)] -> beta
    """
    rng = make_rng() if rng is None else rng
    N, P = Gz.shape
    eta = np.zeros(N, dtype=np.float64)
    if beta_main is not None:
        eta += Gz @ beta_main
    if beta_pairs:
        for (i,j), b in beta_pairs.items():
            eta += b * (Gz[:, i] * Gz[:, j])
    if beta_triples:
        for (i,j,k), b in beta_triples.items():
            eta += b * (Gz[:, i] * Gz[:, j] * Gz[:, k])
    p = expit(eta)
    y = rng.binomial(1, p, size=N).astype(np.float32)
    return y

# ---------- FIGHI score/FI-gain (null = intercept only) ----------
def score_fi_gain_vector(y, z):
    """
    Score-based FI gain for a single candidate feature z (vector).
    Logistic (binary) with intercept-only null.
    Returns: (fi_gain, chi, pval)
    """
    y = y.astype(np.float64)
    z = z.astype(np.float64)
    p0 = np.full_like(y, y.mean())
    W  = np.clip(p0*(1-p0), 1e-9, None)
    U = float(np.dot(z, (y - p0)))
    I = float(np.dot(z*W, z))
    chi = (U*U)/max(I, 1e-12)
    pval = 1.0 - chi2.cdf(chi, df=1)
    fi_gain = 0.5 * chi  # 0.5 * U^2/I
    return fi_gain, chi, pval

def scan_pairs(y, Gz, pairs):
    out = []
    for (i,j) in pairs:
        z = (Gz[:, i] * Gz[:, j])
        fi, chi, p = score_fi_gain_vector(y, z)
        out.append((i,j,fi,chi,p))
    df = pd.DataFrame(out, columns=["i","j","fi_gain","chi","pval"])
    return df

# ---------- BH FDR ----------
def bh_fdr(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty_like(ranked)
    minq = 1.0
    for i in range(m-1, -1, -1):
        minq = min(minq, ranked[i] * m / (i+1))
        q[i] = minq
    qfull = np.empty_like(q)
    qfull[order] = q
    return qfull

# ---------- LD injection (simple) ----------
def inject_ld(G, blocks=10, rho=0.8, rng=None):
    """
    Crude LD: within each block, linearly mix columns to induce correlation.
    """
    rng = make_rng() if rng is None else rng
    N, P = G.shape
    block_size = P // blocks
    G = G.astype(np.float64)
    for b in range(blocks):
        lo = b*block_size
        hi = P if b==blocks-1 else (b+1)*block_size
        k  = hi - lo
        if k <= 1: continue
        # AR(1)-like mixing
        for j in range(lo+1, hi):
            G[:, j] = rho*G[:, j-1] + (1-rho)*G[:, j]
    # re-quantize and zscore
    return zscore_cols(G)
