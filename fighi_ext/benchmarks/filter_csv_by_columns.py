#!/usr/bin/env python3
import argparse, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--csv", required=True)
ap.add_argument("--pheno", required=True)
ap.add_argument("--keep", required=True, help="file with one column name per line")
ap.add_argument("--out", required=True)
args = ap.parse_args()

with open(args.keep) as f:
    cols = [c.strip().strip('"').strip("'") for c in f if c.strip()]

usecols = ["IID", args.pheno] + cols
df = pd.read_csv(args.csv, usecols=lambda c: c in set(usecols))
df.to_csv(args.out, index=False)
print(f"[OK] wrote {args.out} with {df.shape[0]} rows, {df.shape[1]} cols")
