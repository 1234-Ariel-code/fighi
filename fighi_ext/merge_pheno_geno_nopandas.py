#!/usr/bin/env python3
import argparse, csv, sys, math

def mean_impute_col(col):
    vals = [float(x) for x in col if x != ""]
    if not vals: return None
    m = sum(vals)/len(vals)
    return [ (m if x=="" else float(x)) for x in col ]

def main():
    ap = argparse.ArgumentParser(description="Merge a large genotype CSV (IID, rsIDs...) with PLINK .phen (FID IID PHENO) without pandas.")
    ap.add_argument("--geno_csv", required=True)
    ap.add_argument("--phen_file", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--iid_col", default="IID")
    ap.add_argument("--phen_name", default="case")
    ap.add_argument("--impute_mean", action="store_true", help="Mean-impute missing genotype cells per SNP")
    ap.add_argument("--decimal", default=".", help="Decimal point")
    ap.add_argument("--delimiter", default=",")
    args = ap.parse_args()

    # Load phenotype into dict (IID -> 0/1)
    pheno = {}
    with open(args.phen_file, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3: continue
            iid, ph = parts[1], parts[2]
            if ph in ("-9","NA","nan"): continue
            if ph in ("1","2"):   # PLINK 1/2 -> 0/1
                ph = "0" if ph=="1" else "1"
            pheno[iid] = ph

    # Stream genotype CSV: write output on the fly to avoid large memory
    with open(args.geno_csv, "r", newline="") as fin, open(args.out_csv, "w", newline="") as fout:
        rin, wout = csv.reader(fin, delimiter=args.delimiter), csv.writer(fout)
        header = next(rin)
        try:
            iid_idx = header.index(args.iid_col)
        except ValueError:
            sys.exit(f"[merge] IID column '{args.iid_col}' not found in genotype header.")
        out_header = header[:] + [args.phen_name]
        wout.writerow(out_header)

        # If impute requested, we need one pass to compute per-column means.
        means = None
        if args.impute_mean:
            col_sums = [0.0]*len(header)
            col_counts = [0]*len(header)
            # skip header – we will re-open the file
            for row in rin:
                for j in range(len(header)):
                    if j == iid_idx: continue
                    x = row[j].strip()
                    if x=="":
                        continue
                    try:
                        col_sums[j] += float(x.replace(",", "." if args.decimal=="." else args.decimal))
                        col_counts[j] += 1
                    except Exception:
                        pass
            means = [ (col_sums[j]/col_counts[j] if col_counts[j]>0 else None) for j in range(len(header)) ]
            fin.seek(0); rin = csv.reader(fin, delimiter=args.delimiter)
            header = next(rin)

        for row in rin:
            iid = row[iid_idx]
            if iid not in pheno: 
                continue
            # impute on the fly
            if means is not None:
                for j in range(len(header)):
                    if j == iid_idx: continue
                    if row[j].strip()=="" and means[j] is not None:
                        row[j] = f"{means[j]:.6f}"
            wout.writerow(row + [pheno[iid]])
    print(f"[merge] Wrote: {args.out_csv}")

if __name__ == "__main__":
    main()

