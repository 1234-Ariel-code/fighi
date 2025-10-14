#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--cols_file", required=True, help="one SNP per line")
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    cols = [ln.strip() for ln in open(args.cols_file) if ln.strip()]
    usecols = [args.pheno] + cols
    # tolerate mixed dtypes, coerce numerics; keep as float32 to constrain memory
    df = pd.read_csv(args.csv, usecols=usecols)
    for c in df.columns:
        if c == args.pheno:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("float32")
        else:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("float32")
    df.to_csv(args.out_csv, index=False)
    print(f"[subset] wrote {args.out_csv} with {len(cols)} SNPs")

if __name__ == "__main__":
    main()
