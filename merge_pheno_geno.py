#!/usr/bin/env python3
import argparse, sys
import numpy as np
import pandas as pd

def main():
    ap = argparse.ArgumentParser(description="Merge genotype CSV and PLINK .phen by IID")
    ap.add_argument("--geno_csv", required=True, help="Genotype matrix CSV (rows=samples, cols=markers; must contain IID or row order matches .phen)")
    ap.add_argument("--phen_file", required=True, help="PLINK .phen with FID IID PHENO (whitespace-separated, no header)")
    ap.add_argument("--out_csv", required=True, help="Output merged CSV")
    ap.add_argument("--iid_col", default="IID", help="IID column name in genotype CSV (if absent, script will inject one from .phen)")
    ap.add_argument("--phen_name", default="case", help="Name for phenotype column in output CSV")
    args = ap.parse_args()

    # Phenotype
    ph = pd.read_csv(args.phen_file, delim_whitespace=True, header=None, names=["FID","IID","PHENO"])
    # map phenotype {1,2} -> {0,1}; keep 0/1 if already coded; drop missing (-9)
    ph = ph.replace({"PHENO": {-9: np.nan}})
    # If coded 1/2, convert to 0/1; if already 0/1, leave as is.
    uniq = set(pd.unique(ph["PHENO"].dropna()))
    if uniq.issubset({0,1}):
        ph[args.phen_name] = ph["PHENO"].astype("Int64")
    elif uniq.issubset({1,2}):
        ph[args.phen_name] = ph["PHENO"].map({1:0, 2:1}).astype("Int64")
    else:
        sys.stderr.write(f"[merge] Unexpected PHENO values {uniq}. Expect {0,1} or {1,2}. Aborting.\n")
        sys.exit(2)
    ph = ph.dropna(subset=[args.phen_name]).astype({args.phen_name:"int64"})
    ph = ph[["IID", args.phen_name]]

    # Genotypes
    geno = pd.read_csv(args.geno_csv)
    if args.iid_col not in geno.columns:
        # Inject IID from phenotype order if shapes match
        if len(geno) != len(ph):
            sys.stderr.write("[merge] Genotype CSV has no IID column and sample counts differ from .phen.\n")
            sys.exit(2)
        geno.insert(0, args.iid_col, ph["IID"].values)

    # Merge inner join on IID to ensure alignment
    merged = geno.merge(ph, left_on=args.iid_col, right_on="IID", how="inner")
    # Drop duplicate IID column if present twice
    if args.iid_col != "IID" and "IID" in merged.columns:
        merged = merged.drop(columns=["IID"])

    # Optional: ensure numeric genotypes (coerce NAs)
    gcols = [c for c in merged.columns if c not in (args.iid_col, args.phen_name)]
    merged[gcols] = merged[gcols].apply(pd.to_numeric, errors="coerce")

    # Drop rows with all-NA genotypes or missing phenotype (shouldn’t remain)
    merged = merged.dropna(subset=[args.phen_name])

    merged.to_csv(args.out_csv, index=False)
    print(f"[merge] Wrote {args.out_csv} with shape {merged.shape}")

if __name__ == "__main__":
    main()

