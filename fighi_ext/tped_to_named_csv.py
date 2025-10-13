#!/usr/bin/env python3
import argparse, csv, sys, gzip, os
from collections import defaultdict

def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode, newline="")

def main():
    ap = argparse.ArgumentParser(description="Convert PLINK .tped (+optional .tfam) to sample x SNP CSV with rsID cols")
    ap.add_argument("--tped", required=True)
    ap.add_argument("--tfam", help="Optional .tfam: if given, we take IID order from here; else we will ask for .phen")
    ap.add_argument("--phen_file", help="If no TFAM, we derive IID order from .phen (FID IID PHENO)")
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--ref_mode", choices=["major","first"], default="major",
                    help="How to choose reference allele for genotype encoding (AA=0, AB=1, BB=2)")
    ap.add_argument("--chunk", type=int, default=200, help="Write columns in blocks to avoid large RAM")
    args = ap.parse_args()

    # Determine IIDs (sample order)
    iids = []
    if args.tfam and os.path.exists(args.tfam):
        with open_maybe_gz(args.tfam) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    iids.append(parts[1])
    elif args.phen_file:
        with open_maybe_gz(args.phen_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    iids.append(parts[1])
    else:
        sys.exit("[tped_to_named_csv] Need either --tfam or --phen_file to define IID order.")

    n = len(iids)
    # We'll write header first (IID + later SNPs in blocks)
    out = open(args.out_csv, "w", newline="")
    w = csv.writer(out)
    w.writerow(["IID"])  # start with only IID; we’ll append columns by rewriting temp files per block

    # Prepare a temp file to accumulate columns, block-wise
    temp_rows = [[iid] for iid in iids]
    block_names = []
    written_any = False

    def flush_block():
        nonlocal temp_rows, block_names, written_any
        if not block_names:
            return
        # If first block, overwrite file with full header; else append columns.
        if not written_any:
            with open(args.out_csv, "w", newline="") as fout:
                ww = csv.writer(fout)
                ww.writerow(["IID"] + block_names)
                for row in temp_rows:
                    ww.writerow(row)
            written_any = True
        else:
            # read current CSV, append block names as header and the block values to each row
            import pandas as pd
            df = pd.read_csv(args.out_csv)
            for j, name in enumerate(block_names):
                df[name] = [r[len(r)-len(block_names)+j] for r in temp_rows]
            df.to_csv(args.out_csv, index=False)
        # reset block buffer
        temp_rows = [[iid] for iid in iids]
        block_names = []

    # Parse TPED streaming
    with open_maybe_gz(args.tped) as f:
        for ln, line in enumerate(f, 1):
            parts = line.strip().split()
            if len(parts) < 6 or (len(parts)-4) % 2 != 0:
                continue
            chrom, rsid, cm, bp = parts[:4]
            alleles = parts[4:]
            # genotype coding per sample
            a_counts = defaultdict(int)
            geno = [0]*n
            for i in range(n):
                a1, a2 = alleles[2*i], alleles[2*i+1]
                if a1 == "0" or a2 == "0":
                    geno[i] = ""  # NA
                    continue
                if args.ref_mode == "major":
                    a_counts[a1] += 1; a_counts[a2] += 1
                geno[i] = (0 if a1 == a2 else 1)  # temporary; will finalize after ref allele set
            ref = None
            if args.ref_mode == "major":
                if a_counts:
                    ref = max(a_counts.items(), key=lambda kv: kv[1])[0]
            else:
                ref = alleles[0] if alleles[0] != "0" else alleles[1]
            # finalize 0/1/2 relative to ref (ref/ref=0, ref/alt=1, alt/alt=2)
            for i in range(n):
                a1, a2 = alleles[2*i], alleles[2*i+1]
                if a1 == "0" or a2 == "0":
                    val = ""
                else:
                    if a1 == ref and a2 == ref:
                        val = 0
                    elif a1 == ref or a2 == ref:
                        val = 1
                    else:
                        val = 2
                temp_rows[i].append(val)
            block_names.append(rsid)
            if len(block_names) >= args.chunk:
                flush_block()

    flush_block()
    out.close()
    print(f"[tped_to_named_csv] Wrote: {args.out_csv}")

if __name__ == "__main__":
    main()

