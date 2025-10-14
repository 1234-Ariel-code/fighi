#!/usr/bin/env python3
import argparse, os, json, sys
import numpy as np, pandas as pd

# --- robust imports: works as module or script ---
try:
    from .pipeline import FIGHIPipeline
    from .report import write_feature_scores, write_summary, save_model, write_log, make_plots
except Exception:
    sys.path.append(os.path.dirname(__file__))
    from pipeline import FIGHIPipeline
    from report import write_feature_scores, write_summary, save_model, write_log, make_plots


# ----------------- helpers -----------------
def _clean_name(name: str) -> str:
    """Trim whitespace and any surrounding single/double quotes."""
    if name is None:
        return name
    return name.strip().strip('"').strip("'")

def _read_header(csv_path: str) -> list[str]:
    """Read raw header safely and normalize each column name."""
    with open(csv_path, "r", encoding="utf-8") as f:
        first = f.readline().rstrip("\n")
    # Let pandas parse just the header line to be safe with commas/quotes
    # but avoid reading whole file: use a one-row DataFrame trick.
    raw_cols = [c for c in first.split(",")]
    return [_clean_name(c) for c in raw_cols]

def _downcast_df(df: pd.DataFrame, pheno: str, id_col: str | None) -> pd.DataFrame:
    for c in df.columns:
        if id_col and c == id_col:
            continue  # keep IDs/strings unchanged
        # convert to numeric, NaNs on failure, then float32
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("float32")
    return df

def _probe_numeric_cols(csv_path: str, candidates: list[str], nrows: int = 256) -> list[str]:
    """
    Probe a tiny sample to decide which columns are numeric-like.
    If 'usecols' fails due to header quoting/formatting, fall back to reading all
    columns for nrows and intersect afterward.
    """
    cand_set = set(candidates)
    try:
        # Primary path: read only the candidate columns for a tiny sample
        df = pd.read_csv(csv_path, usecols=candidates, nrows=nrows)
        keep = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        return keep
    except Exception:
        # Fallback: read small sample of ALL columns, then intersect
        sample = pd.read_csv(csv_path, nrows=nrows)
        sample.columns = [_clean_name(c) for c in sample.columns]
        keep = [c for c in sample.columns
                if (c in cand_set) and pd.api.types.is_numeric_dtype(sample[c])]
        return keep

def _compute_y(csv_path: str, pheno: str, chunksize: int) -> np.ndarray:
    parts = []
    for chunk in pd.read_csv(csv_path, usecols=[pheno], dtype={pheno: "float32"}, chunksize=chunksize):
        parts.append(chunk[pheno].to_numpy(dtype=np.float32, copy=False))
    return np.concatenate(parts, axis=0)

def _prescreen_columns(csv_path: str,
                       pheno: str,
                       id_col: str | None,
                       top_m: int = 200,
                       col_block: int = 500,
                       row_chunksize: int = 50000):
    """
    Streaming prescreen using absolute correlation with y.
    Robust to non-numeric columns (filters them) and optional ID column.
    """
    header = _read_header(csv_path)
    if pheno not in header:
        raise ValueError(f"Phenotype column '{pheno}' not found in CSV header. "
                         f"First 10 columns: {header[:10]}")

    # candidates = all columns minus phenotype and optional ID, cleaned
    snp_cols = [c for c in header if c != pheno and (not id_col or c != id_col)]

    # Keep only numeric-looking columns by sampling a few rows
    snp_cols = _probe_numeric_cols(csv_path, snp_cols, nrows=256)
    if not snp_cols:
        raise ValueError("No numeric SNP columns found after probing. "
                         "Check that your CSV uses numeric encodings (0/1/2) and the header isn’t malformed.")

    y = _compute_y(csv_path, pheno, row_chunksize).astype(np.float32)
    n = float(y.size)
    y_mean = float(np.nanmean(y))
    y_var  = max(1e-12, float(np.nanvar(y)))
    y_center = y - y_mean

    scores: dict[str, float] = {}
    for i in range(0, len(snp_cols), col_block):
        cols = snp_cols[i:i + col_block]

        # pass A: means/vars
        sum_x  = {c: 0.0 for c in cols}
        sum_x2 = {c: 0.0 for c in cols}
        for chunk in pd.read_csv(csv_path,
                                 usecols=[pheno] + cols,
                                 dtype={pheno: "float32", **{c: "float32" for c in cols}},
                                 chunksize=row_chunksize):
            arr = chunk[cols].to_numpy(dtype=np.float32, copy=False)
            arr = np.where(np.isfinite(arr), arr, np.nan)
            sum_x_block  = np.nansum(arr, axis=0)
            sum_x2_block = np.nansum(arr * arr, axis=0)
            for j, c in enumerate(cols):
                sum_x[c]  += float(sum_x_block[j])
                sum_x2[c] += float(sum_x2_block[j])
        mean_x = {c: sum_x[c] / max(1.0, n) for c in cols}

        # pass B: covariance with centered y (mean-impute per col)
        cov_xy = {c: 0.0 for c in cols}
        y_pos = 0
        for chunk in pd.read_csv(csv_path,
                                 usecols=[pheno] + cols,
                                 dtype={pheno: "float32", **{c: "float32" for c in cols}},
                                 chunksize=row_chunksize):
            m = len(chunk)
            yc = y_center[y_pos:y_pos + m]
            arr = chunk[cols].to_numpy(dtype=np.float32, copy=False)
            arr = np.where(np.isfinite(arr), arr, np.nan)
            for j, c in enumerate(cols):
                xj = arr[:, j]
                mu = mean_x[c]
                xj = np.where(np.isfinite(xj), xj, mu)
                cov_xy[c] += float(np.dot(xj - mu, yc))
            y_pos += m

        # finalize correlations
        for c in cols:
            mu = mean_x[c]
            var_x = max(1e-12, (sum_x2[c] / max(1.0, n)) - mu * mu)
            corr = cov_xy[c] / max(1e-12, np.sqrt(var_x * y_var) * (n - 1.0))
            scores[c] = abs(corr)

    top_cols = [c for c, _ in sorted(scores.items(), key=lambda kv: -kv[1])[:top_m]]
    return top_cols


# ----------------- main -----------------
def main():
    ap = argparse.ArgumentParser(description="FIGHI CLI (memory-safe)")
    ap.add_argument("--csv", required=True, help="CSV with phenotype column and features")
    ap.add_argument("--pheno", required=True, help="Phenotype column name (e.g., case)")
    ap.add_argument("--id_col", default="IID", help="Optional ID column to ignore (default: IID)")
    ap.add_argument("--trait", default="binary", choices=["binary","linear"])
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--max_order", type=int, default=4)
    ap.add_argument("--target_or", type=float, default=1.3)

    # memory controls
    ap.add_argument("--read_chunksize", type=int, default=0,
                    help="If >0, read CSV in row chunks then concat (still loads all columns).")
    ap.add_argument("--prescreen_top_m", type=int, default=0,
                    help="If >0, prescreen columns in column blocks and only load the top-M SNPs.")
    ap.add_argument("--col_block", type=int, default=500,
                    help="Number of SNP columns to scan per block during prescreen.")
    ap.add_argument("--row_chunksize", type=int, default=50000,
                    help="Row chunk size while streaming during prescreen and safe reads.")
    ap.add_argument("--write_top_cols", action="store_true",
                    help="If set with prescreen, write ranked SNP list to outdir/top_columns.txt")
    ap.add_argument("--no_plots", action="store_true",
                    help="Disable plotting (useful on headless/HPC nodes).")

    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Normalize names passed by user (defensive)
    args.pheno  = _clean_name(args.pheno)
    args.id_col = _clean_name(args.id_col) if args.id_col else args.id_col

    # --------- Mode A: prescreen then load only top-M ----------
    if args.prescreen_top_m and args.prescreen_top_m > 0:
        print(f"[INFO] Prescreening to top {args.prescreen_top_m} columns "
              f"(block={args.col_block}, rowsz={args.row_chunksize})")
        top_cols = _prescreen_columns(args.csv, args.pheno, args.id_col,
                                      top_m=args.prescreen_top_m,
                                      col_block=args.col_block,
                                      row_chunksize=args.row_chunksize)
        if args.write_top_cols:
            with open(os.path.join(args.outdir, "top_columns.txt"), "w") as f:
                f.write("\n".join(top_cols) + "\n")
        usecols = [args.pheno] + ([args.id_col] if args.id_col else []) + top_cols
        # Let pandas infer then we normalize types
        df = pd.read_csv(args.csv, usecols=usecols)
        df.columns = [_clean_name(c) for c in df.columns]
        df = _downcast_df(df, args.pheno, args.id_col)

    # --------- Mode B: row-chunked read of full CSV -----------
    elif args.read_chunksize and args.read_chunksize > 0:
        parts = []
        for chunk in pd.read_csv(args.csv, chunksize=args.read_chunksize):
            chunk.columns = [_clean_name(c) for c in chunk.columns]
            parts.append(_downcast_df(chunk, args.pheno, args.id_col))
        df = pd.concat(parts, axis=0, copy=False)
        del parts

    # --------- Mode C: simple read (small data) ---------------
    else:
        df = pd.read_csv(args.csv)
        df.columns = [_clean_name(c) for c in df.columns]
        df = _downcast_df(df, args.pheno, args.id_col)

    # sanity checks & drop ID
    if args.pheno not in df.columns:
        raise ValueError(f"Phenotype column '{args.pheno}' not found. "
                         f"Columns: {list(df.columns)[:20]} ...")
    if args.id_col and args.id_col in df.columns:
        df = df.drop(columns=[args.id_col])

    # Build X, y and run
    y = df[args.pheno].to_numpy(dtype=np.float32, copy=False)
    X = df.drop(columns=[args.pheno])

    pipe = FIGHIPipeline(trait_type=args.trait, n_perm=0)
    res = pipe.adaptive_search(X, y, max_order=args.max_order, target_OR=args.target_or)
    paths = pipe.export(args.outdir, res)

    write_feature_scores(args.outdir, res["results"], X)
    write_summary(args.outdir, res, X, y, args.max_order)
    save_model(args.outdir, res)
    write_log(args.outdir, X, y, res)

    if not args.no_plots:
        os.environ.setdefault("MPLBACKEND", "Agg")
        make_plots(args.outdir, res, X)

    print(json.dumps(paths, indent=2))


if __name__ == "__main__":
    main()
