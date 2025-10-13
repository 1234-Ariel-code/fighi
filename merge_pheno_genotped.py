#!/usr/bin/env python3
import argparse, sys
import numpy as np
import pandas as pd

def _clean_and_map_pheno(ph: pd.DataFrame, phen_name: str) -> pd.DataFrame:
    ph = ph.replace({"PHENO": {-9: np.nan}})
    uniq = set(pd.unique(ph["PHENO"].dropna()))
    if not uniq:
        sys.stderr.write("[merge] No non-missing phenotypes after dropping -9.\n"); sys.exit(2)
    if uniq.issubset({0,1}):
        ph[phen_name] = ph["PHENO"].astype("Int64")
    elif uniq.issubset({1,2}):
        ph[phen_name] = ph["PHENO"].map({1:0, 2:1}).astype("Int64")
    else:
        sys.stderr.write(f"[merge] Unexpected PHENO values {sorted(uniq)}; expected {{0,1}} or {{1,2}}.\n"); sys.exit(2)
    ph = ph.dropna(subset=[phen_name]).astype({phen_name:"int64"})[["IID", phen_name]].drop_duplicates("IID")
    return ph

def main():
    ap = argparse.ArgumentParser(description="Merge genotype CSV and PLINK .phen by IID.")
    ap.add_argument("--geno_csv", required=True, help="Genotype CSV with IID + rsIDs")
    ap.add_argument("--phen_file", required=True, help=".phen (FID IID PHENO)")
    ap.add_argument("--out_csv", required=True, help="Output merged CSV")
    ap.add_argument("--iid_col", default="IID", help="IID column name in genotype CSV")
    ap.add_argument("--phen_name", default="case", help="Output phenotype column name")
    ap.add_argument("--impute_mean", action="store_true", help="Mean-impute missing genotypes")
    args = ap.parse_args()

    ph = pd.read_csv(args.phen_file, delim_whitespace=True, header=None, names=["FID","IID","PHENO"])
    ph = _clean_and_map_pheno(ph, args.phen_name)

    geno = pd.read_csv(args.geno_csv)
    if args.iid_col not in geno.columns:
        sys.stderr.write("[merge] Genotype CSV must contain an IID column.\n"); sys.exit(2)
    geno = geno.drop_duplicates(subset=[args.iid_col])

    merged = geno.merge(ph, left_on=args.iid_col, right_on="IID", how="inner")
    if args.iid_col != "IID" and "IID" in merged.columns:
        merged = merged.drop(columns=["IID"])

    non_geno = {args.iid_col, args.phen_name}
    gcols = [c for c in merged.columns if c not in non_geno]
    merged[gcols] = merged[gcols].apply(pd.to_numeric, errors="coerce")

    if args.impute_mean and gcols:
        means = merged[gcols].mean()
        merged[gcols] = merged[gcols].fillna(value=means.to_dict())

    merged = merged[[args.iid_col, args.phen_name] + [c for c in merged.columns if c not in (args.iid_col, args.phen_name)]]
    merged.to_csv(args.out_csv, index=False)
    print(f"[merge] Wrote {args.out_csv} with shape {merged.shape}")

if __name__ == "__main__":
    main()

