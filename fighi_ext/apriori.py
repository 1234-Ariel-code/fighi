from itertools import combinations
from typing import List, Tuple, Set

def apriori_generate(prev_level: List[Tuple[str, ...]]) -> List[Tuple[str, ...]]:
    k = len(prev_level[0]) if prev_level else 0
    C = set()
    for a in prev_level:
        for b in prev_level:
            if a[:-1] == b[:-1] and a[-1] < b[-1]:
                C.add(tuple(sorted(a + (b[-1],))))
    return sorted(C)

def prune_by_subsets(cands: List[Tuple[str, ...]], prev_level: Set[Tuple[str, ...]]) -> List[Tuple[str, ...]]:
    from itertools import combinations as comb
    pruned = []
    for c in cands:
        ok = True
        for s in comb(c, len(c)-1):
            if tuple(sorted(s)) not in prev_level:
                ok = False
                break
        if ok:
            pruned.append(c)
    return pruned

def prune_by_variance(X, cands, min_var=1e-6):
    kept = []
    import numpy as np
    for c in cands:
        z = np.ones(len(X), dtype=float)
        for col in c:
            z *= X[col].values
        if np.var(z) >= min_var:
            kept.append(c)
    return kept
