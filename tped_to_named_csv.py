#!/usr/bin/env python3
import argparse, sys
import numpy as np
import pandas as pd
from collections import Counter

def choose_ref(alleles_flat, mode="major"):
    vals = [a for a in alleles_flat if a in ("A","C","G","T")]
    if not vals:
        return "A"
    if mode == "first":
        return vals[0]
    return Counter(vals).most_common(1)[0][0]  # major

def encode_row_to_dosage(alleles_row, n, ref):
    g = alleles_row.reshape(n, 2)
    out = np.empty(n, dtype=float)
    for i, (a1,a2) in enumerate(g):
        if a1 == "0" or a2 == "0":
            out[i] = np.nan
        elif a1 == a2:
            out[i] = 0.0 if a1 == ref else 2.0
        else:
            out[i] = 1.0
    return out

def main():
    ap = argparse.ArgumentParser(
        description="Convert PLINK .tped to numeric 0/1/2 CSV with rsIDs and IID column."
    )
    ap.add_argument("--tped", required=True, help=".tped file")
    ap.add_argument("--tfam", default=None,
                    help=".tfam to define IID order (preferred)")
    ap.add_argument("--phen_file", required=True,
                    help=".phen (FID IID PHENO) for fallback IID order")
    ap.add_argument("--out_csv", required=True, help="Output genotype CSV with IID + rsIDs")
    ap.add_argument("--ref_mode", choices=["major","first"], default="major",
                    help="Reference allele rule per SNP (default: major)")
    args = ap.parse_args()

    # Infer IID order
    if args.tfam:
        tfam = pd.read_csv(args.tfam, delim_whitespace=True, header=None,
                           names=["FID","IID","PAT","MAT","SEX","PHENO"])
        iids = tfam["IID"].astype(str).tolist()
    else:
        ph = pd.read_csv(args.phen_file, delim_whitespace=True, header=None,
                         names=["FID","IID","PHENO"])
        iids = ph["IID"].astype(str).tolist()
    n = len(iids)

    snp_names, dosage_cols = [], []
    with open(args.tped, "r") as fh:
        for line_no, line in enumerate(fh, 1):
            parts = line.rstrip("\n").split()
            if len(parts) < 4 + 2*n:
                sys.stderr.write(f"[tped] Line {line_no}: expected >= {4+2*n}, got {len(parts)}.\n")
                sys.exit(2)
            chrom, snp, gd, bp = parts[:4]
            alleles_np = np.array(parts[4:], dtype=object)
            ref = choose_ref(alleles_np, mode=args.ref_mode)
            dosage = encode_row_to_dosage(alleles_np, n=n, ref=ref)
            snp_names.append(snp)
            dosage_cols.append(dosage)

    geno = pd.DataFrame(np.column_stack(dosage_cols), columns=snp_names)
    geno.insert(0, "IID", iids)
    geno.to_csv(args.out_csv, index=False)
    print(f"[tped] Wrote {args.out_csv} with shape {geno.shape} (samples={n}, snps={len(snp_names)})")

if __name__ == "__main__":
    main()

