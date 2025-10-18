#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np, subprocess, os, tempfile, shutil

"""
Create a minimal PLINK dataset from a FIGHI CSV (IID, pheno, rs...).
We fabricate alleles as A/G and map dosage:
  0 -> A A
  1 -> A G
  2 -> G G
We also create a dummy MAP with chr=1 and bp=index.
"""

ap = argparse.ArgumentParser()
ap.add_argument("--csv", required=True)
ap.add_argument("--pheno", required=True)
ap.add_argument("--out_prefix", required=True)
args = ap.parse_args()

df = pd.read_csv(args.csv)
assert "IID" in df.columns and args.pheno in df.columns, "CSV must contain IID and phenotype"

# Build MAP
snp_cols = [c for c in df.columns if c not in ("IID", args.pheno)]
map_lines = []
for i, snp in enumerate(snp_cols, start=1):
    map_lines.append(f"1\t{snp}\t0\t{i}\n")

tmpdir = tempfile.mkdtemp()
ped_path = os.path.join(tmpdir, "cd.ped")
map_path = os.path.join(tmpdir, "cd.map")
with open(map_path, "w") as f:
    f.writelines(map_lines)

# PLINK PED: FID IID PAT MAT SEX PHENO followed by genotypes (pairs)
# PHENO must be 1/2 (1=control, 2=case); here pheno is 0/1, so +1
def encode_pair(v):
    if np.isnan(v): return ("0","0")
    v = int(v)
    if v==0: return ("A","A")
    if v==1: return ("A","G")
    if v==2: return ("G","G")
    return ("0","0")

with open(ped_path, "w") as f:
    for _, row in df.iterrows():
        iid = row["IID"]
        ph = int(row[args.pheno])
        ph = 2 if ph==1 else 1   # PLINK 1/2 coding
        header = [iid, iid, "0", "0", "1", str(ph)]
        gt_pairs = []
        for c in snp_cols:
            a1,a2 = encode_pair(row[c])
            gt_pairs.extend([a1,a2])
        f.write(" ".join(header + gt_pairs) + "\n")

outp = args.out_prefix
# Convert PED/MAP -> BED
subprocess.check_call(["plink", "--file", os.path.splitext(ped_path)[0],
                       "--make-bed", "--allow-no-sex", "--out", outp],
                      stdout=subprocess.DEVNULL)

shutil.rmtree(tmpdir)
print(f"[OK] Wrote PLINK: {outp}.bed/.bim/.fam")
