import numpy as np

# --- robust import for jaccard ---
try:
    from .utils import jaccard
except Exception:
    import os, sys
    sys.path.append(os.path.dirname(__file__))
    from utils import jaccard
# ----------------------------------


#from .utils import jaccard

def stability_selection(run_fn, n_boot=50, subsample=0.7, rng=0):
    rng = np.random.RandomState(rng)
    supports = []
    for b in range(n_boot):
        mask = rng.rand(run_fn.N) < subsample
        res = run_fn(mask)  # expected to return set/list of selected tuples
        supports.append(set(res))
    # compute selection probabilities
    all_items = sorted(set().union(*supports))
    probs = {}
    for item in all_items:
        probs[item] = np.mean([item in s for s in supports])
    return probs
