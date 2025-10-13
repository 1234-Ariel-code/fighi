#!/usr/bin/env python3
import argparse, os, json
import numpy as np, pandas as pd

# --- robust imports: works as module or script ---
try:
    # module mode: python -m fighi_ext.run_cli
    from .pipeline import FIGHIPipeline
    from .report import write_feature_scores, write_summary, save_model, write_log, make_plots
except Exception:
    # script mode: python run_cli.py from inside fighi_ext/
    import sys
    sys.path.append(os.path.dirname(__file__))
    from pipeline import FIGHIPipeline
    from report import write_feature_scores, write_summary, save_model, write_log, make_plots


def main():
    ap = argparse.ArgumentParser(description="FIGHI CLI")
    ap.add_argument("--csv", required=True, help="CSV with phenotype column and features")
    ap.add_argument("--pheno", required=True, help="Phenotype column name")
    ap.add_argument("--trait", default="binary", choices=["binary","linear"])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--max_order", type=int, default=4)
    ap.add_argument("--target_or", type=float, default=1.3)
    ap.add_argument("--read_chunksize", type=int, default=0, help="If >0, read CSV in chunks & assemble (for huge files)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Memory-aware read
    if args.read_chunksize and args.read_chunksize > 0:
        parts = []
        for chunk in pd.read_csv(args.csv, chunksize=args.read_chunksize):
            parts.append(chunk)
        df = pd.concat(parts, axis=0, copy=False)
        del parts
    else:
        df = pd.read_csv(args.csv)

    y = df[args.pheno].values.astype(float)
    X = df.drop(columns=[args.pheno])

    pipe = FIGHIPipeline(trait_type=args.trait, n_perm=0)  # keep light; can raise if you enable perms
    res = pipe.adaptive_search(X, y, max_order=args.max_order, target_OR=args.target_or)
    paths = pipe.export(args.outdir, res)

    # “Nice” deliverables
    write_feature_scores(args.outdir, res["results"], X)
    write_summary(args.outdir, res, X, y, args.max_order)
    save_model(args.outdir, res)
    write_log(args.outdir, X, y, res)
    make_plots(args.outdir, res, X)

    print(json.dumps(paths, indent=2))

if __name__ == "__main__":
    main()




















#import argparse, pandas as pd, numpy as np, os, json
#from fighi_ext.pipeline import FIGHIPipeline

#def main():
#    ap = argparse.ArgumentParser(description="FIGHI CLI")
#    ap.add_argument("--csv", required=True, help="CSV with phenotype column and features")
#    ap.add_argument("--pheno", required=True, help="Phenotype column name")
#    ap.add_argument("--trait", default="binary", choices=["binary","linear"])
#    ap.add_argument("--outdir", required=True)
#    #ap.add_argument("--max_order", type=int, default=4)
#    ap.add_argument("--kmax", type=int, help=argparse.SUPPRESS)
#    ap.add_argument("--target_or", type=float, default=1.3)
#    args = ap.parse_args()

#    df = pd.read_csv(args.csv)
#    y = df[args.pheno].values.astype(float)
#    X = df.drop(columns=[args.pheno])
#    pipe = FIGHIPipeline(trait_type=args.trait, n_perm=100)
#    res = pipe.adaptive_search(X, y, max_order=args.max_order, target_OR=args.target_or)
#    paths = pipe.export(args.outdir, res)
#    print(json.dumps(paths, indent=2))

#if __name__ == "__main__":
#    main()
