#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np
from itertools import combinations
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from mdr import MDR

ap = argparse.ArgumentParser()
ap.add_argument("--csv", required=True)
ap.add_argument("--pheno", default="case")
ap.add_argument("--out", required=True)
ap.add_argument("--cv_folds", type=int, default=5)
ap.add_argument("--max_pairs", type=int, default=10_000_000)
args = ap.parse_args()

df = pd.read_csv(args.csv)
y = df[args.pheno].astype(int).values
X = df.drop(columns=["IID", args.pheno]).values
names = [c for c in df.columns if c not in ("IID", args.pheno)]
p = X.shape[1]

skf = StratifiedKFold(n_splits=args.cv_folds, shuffle=True, random_state=42)

rows = []
n_scanned = 0
for i,j in combinations(range(p), 2):
    if n_scanned >= args.max_pairs: break
    Xi = X[:, [i,j]]

    aucs = []
    for train_idx, test_idx in skf.split(Xi, y):
        mdl = MDR()     # default: balanced accuracy scoring
        mdl.fit(Xi[train_idx], y[train_idx])
        yhat = mdl.predict(Xi[test_idx])
        # AUC on predicted class; not ideal for MDR, but gives a numeric comparator
        aucs.append(roc_auc_score(y[test_idx], yhat))
    rows.append((names[i], names[j], float(np.mean(aucs))))
    n_scanned += 1

res = pd.DataFrame(rows, columns=["SNP1","SNP2","CV_AUC"])
res.sort_values("CV_AUC", ascending=False, inplace=True)
res.to_csv(args.out, index=False)
print(f"[OK] MDR wrote {args.out}; pairs={len(res)}")
