import numpy as np, pandas as pd
from fighi_ext.pipeline import FIGHIPipeline

rng = np.random.RandomState(0)
N = 1000
X = pd.DataFrame({
    "g1": rng.binomial(2, 0.3, size=N),
    "g2": rng.binomial(2, 0.3, size=N),
    "g3": rng.binomial(2, 0.3, size=N),
    "noise": rng.normal(size=N)
})
z = (X - X.mean())/X.std()
y = (z["g1"]*z["g2"]*z["g3"]) + 0.1*rng.normal(size=N)

pipe = FIGHIPipeline(trait_type="linear", n_perm=50)
res = pipe.adaptive_search(X, y, max_order=4, target_OR=1.3)
paths = pipe.export("/mnt/data/fighi_ext_output_linear", res)
print(paths)
