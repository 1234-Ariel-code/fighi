#!/usr/bin/env python3
import argparse, os, glob, csv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batches_out", required=True)
    ap.add_argument("--out_snplist", required=True)
    args = ap.parse_args()

    snps = set()
    for csvp in glob.glob(os.path.join(args.batches_out, "batch_*", "fighi_results.csv")):
        with open(csvp) as f:
            r = csv.reader(f)
            header = next(r, None)
            if not header:
                continue
            try:
                hidx = header.index("hyperedge")
            except ValueError:
                # fallback to first column
                hidx = 0
            for row in r:
                for s in row[hidx].split("|"):
                    s = s.strip()
                    if s:
                        snps.add(s)

    with open(args.out_snplist, "w") as f:
        for s in sorted(snps):
            f.write(s + "\n")

    print(f"[reconcile] wrote union of {len(snps)} SNPs to {args.out_snplist}")

if __name__ == "__main__":
    main()
