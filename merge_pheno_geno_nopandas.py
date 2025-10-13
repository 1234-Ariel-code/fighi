#!/usr/bin/env python3
import argparse, csv, gzip, sys, os

MISSING_STRS = {"", "NA", "NaN", "nan", "N/A", "null", "NULL"}

def open_maybe_gz(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", newline="")
    return open(path, mode, newline="")

def parse_args():
    ap = argparse.ArgumentParser(
        description="Merge genotype CSV (with IID column) and PLINK .phen by IID, no pandas."
    )
    ap.add_argument("--geno_csv", required=True, help="Genotype CSV or CSV.GZ; must have IID column")
    ap.add_argument("--phen_file", required=True, help="PLINK .phen (FID IID PHENO), whitespace-sep")
    ap.add_argument("--out_csv", required=True, help="Output merged CSV")
    ap.add_argument("--iid_col", default="IID", help="IID column name in genotype CSV")
    ap.add_argument("--phen_name", default="case", help="Output phenotype column name")
    ap.add_argument("--impute_mean", action="store_true", help="Two-pass mean imputation for missing genotypes")
    ap.add_argument("--decimal", default=".", help="Decimal separator in input CSV (default '.')")
    ap.add_argument("--delimiter", default=",", help="Field delimiter in genotype CSV (default ',')")
    return ap.parse_args()

def load_phen(path, phen_name):
    iid2y = {}
    with open(path, "r") as f:
        for line in f:
            if not line.strip(): 
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            _fid, iid, ph = parts[0], parts[1], parts[2]
            if ph == "-9":
                continue
            try:
                v = int(ph)
            except ValueError:
                continue
            # normalize to 0/1
            if v in (0, 1):
                y = v
            elif v in (1, 2):
                y = 0 if v == 1 else 1
            else:
                continue
            iid2y[iid] = y
    if not iid2y:
        sys.stderr.write("[merge] No valid phenotypes found.\n")
        sys.exit(2)
    return iid2y

def try_float(x, decimal="."):
    if x in MISSING_STRS:
        return None
    if decimal != ".":
        x = x.replace(decimal, ".")
    try:
        return float(x)
    except Exception:
        return None

def compute_means(geno_csv, iid_col, delimiter, decimal):
    """First pass: per-column mean for numeric genotype columns."""
    with open_maybe_gz(geno_csv, "r") as f:
        rdr = csv.reader(f, delimiter=delimiter)
        header = next(rdr)
        if iid_col not in header:
            sys.stderr.write(f"[merge] '{iid_col}' not found in header.\n")
            sys.exit(2)
        idx_iid = header.index(iid_col)
        # columns to average = all except IID; keep order
        cols = [i for i, c in enumerate(header) if i != idx_iid]
        sums = [0.0] * len(cols)
        cnts = [0] * len(cols)
        for row in rdr:
            if not row:
                continue
            for j, ci in enumerate(cols):
                v = try_float(row[ci], decimal=decimal)
                if v is not None:
                    sums[j] += v
                    cnts[j] += 1
        means = {}
        for j, ci in enumerate(cols):
            means[ci] = (sums[j] / cnts[j]) if cnts[j] > 0 else 0.0
        return header, idx_iid, means

def merge(geno_csv, out_csv, iid_col, iid2y, phen_name, delimiter, decimal, impute_mean, means):
    n_out = 0
    with open_maybe_gz(geno_csv, "r") as fin, open(out_csv, "w", newline="") as fout:
        rdr = csv.reader(fin, delimiter=delimiter)
        wtr = csv.writer(fout)
        header = next(rdr)
        idx_iid = header.index(iid_col)
        out_header = header + [phen_name]
        wtr.writerow(out_header)
        for row in rdr:
            if not row:
                continue
            iid = row[idx_iid]
            if iid not in iid2y:
                continue  # drop samples without phenotype
            # optional mean impute on numeric genotype columns (all except IID)
            if impute_mean:
                for ci, colname in enumerate(header):
                    if ci == idx_iid:
                        continue
                    if row[ci] in MISSING_STRS:
                        row[ci] = f"{means.get(ci, 0.0):.6g}"
                    else:
                        # validate numeric; if not numeric, try to coerce, else keep as-is
                        v = try_float(row[ci], decimal=decimal)
                        if v is None:
                            row[ci] = f"{means.get(ci, 0.0):.6g}"
                        else:
                            # keep original formatting (already fine)
                            pass
            wtr.writerow(row + [iid2y[iid]])
            n_out += 1
    return n_out

def main():
    a = parse_args()
    iid2y = load_phen(a.phen_file, a.phen_name)
    # Validate header, maybe compute means
    with open_maybe_gz(a.geno_csv, "r") as f:
        rdr = csv.reader(f, delimiter=a.delimiter)
        header = next(rdr)
    if a.iid_col not in header:
        sys.stderr.write(f"[merge] '{a.iid_col}' not found in header.\n")
        sys.exit(2)

    means = {}
    if a.impute_mean:
        header, idx_iid, means = compute_means(a.geno_csv, a.iid_col, a.delimiter, a.decimal)

    n = merge(
        a.geno_csv, a.out_csv, a.iid_col, iid2y, a.phen_name,
        a.delimiter, a.decimal, a.impute_mean, means
    )
    print(f"[merge] Wrote {a.out_csv} with {n} rows.")

if __name__ == "__main__":
    main()

