import numpy as np
import pandas as pd
from itertools import combinations
from typing import Iterable, List, Tuple

def standardize(X: pd.DataFrame) -> pd.DataFrame:
    Z = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0).replace(0, 1.0)
    return Z

def maf(series: pd.Series) -> float:
    x = series.dropna().values
    if x.size == 0:
        return 0.0
    m = np.clip(x, 0, 2)
    p_alt = m.mean() / 2.0
    return float(min(p_alt, 1.0 - p_alt))

def product_feature(X: pd.DataFrame, cols: Tuple[str, ...]) -> np.ndarray:
    z = np.ones(len(X), dtype=float)
    for c in cols:
        z = z * X[c].values
    return z

def variance_of_product(X: pd.DataFrame, cols: Tuple[str, ...]) -> float:
    z = product_feature(X, cols)
    return float(np.var(z, ddof=0))

def jaccard(set_a: Iterable, set_b: Iterable) -> float:
    A, B = set(set_a), set(set_b)
    if not A and not B:
        return 1.0
    return len(A & B) / max(1, len(A | B))

def bh_fdr(pvals, alpha=0.05):
    """Benjamini-Hochberg q-values and threshold index."""
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty_like(ranked)
    min_q = 1.0
    for i in range(m-1, -1, -1):
        min_q = min(min_q, ranked[i] * m / (i+1))
        q[i] = min_q
    q_full = np.empty_like(q)
    q_full[order] = q
    return q_full
