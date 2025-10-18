#!/usr/bin/env python3
import argparse, json, os, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from mdr import MDR

def build_feature_from_tuple(df, cols):
    """Return vector = product of columns (handles NaN as column mean)."""
    z = np.ones(len(df), dtype=np.float32)
    for c in cols:
        x = pd.to_numeric(df[c], errors="coerce").astype(np.float32)
        mu = np.nanmean(x)
        x = np.where(np.isfinite(x), x, mu)
        z = z * x
    return z

def design_from_fighi(df, fighi_csv, top_k=50):
    """Top-K by ΔI; build matrix with singles and interactions as columns."""
    if not os.path.exists(fighi_csv):
        return None, []
    R = pd.read_csv(fighi_csv)
    if "fi_gain" not in R.columns and "FI_gain" in R.columns:
        R = R.rename(columns={"FI_gain":"fi_gain"})
    if "fi_gain" not in R.columns:
        # try weight column
        if "weight" in R.columns: R = R.rename(columns={"weight":"fi_gain"})
    R = R.sort_values("fi_gain", ascending=False).head(top_k)
    feats = []
    X_list = []
    for _, row in R.iterrows():
        # expect columns like: order, tuple (e.g., "rs1|rs2|rs3") or a list
        if "order" not in row or "tuple" not in R.columns:
            # try columns SNP1,SNP2,...
            cols = [c for c in row.index if str(c).lower().startswith("snp")]
            cols = [row[c] for c in cols if pd.notnull(row[c])]
        else:
            cols = str(row["tuple"]).split("|")
        cols = [c.strip() for c in cols if isinstance(c,str) and len(c.strip())>0]
        if not cols: 
            # fallback: if column 'snp' exists for main effect
            if "snp" in row.index and isinstance(row["snp"], str):
                cols = [row["snp"]]
            else:
                continue
        z = build_feature_from_tuple(df, cols)
        X_list.append(z)
        feats.append(tuple(cols))
    if not X_list:
        return None, []
    X = np.vstack(X_list).T  # n x k
    return X, feats

def design_from_plink_pair(df, plink_epi_file):
    """Use lowest p-value pair -> columns A,B,A*B."""
    if not os.path.exists(plink_epi_file):
        return None, ()
    epi = pd.read_csv(plink_epi_file, delim_whitespace=True)
    # Try to find p-value column (FAST epistasis has 'P', logistic has 'P')
    pcol = "P" if "P" in epi.columns else ( "CHISQ" if "CHISQ" in epi.columns else None )
    epi = epi.sort_values(pcol, ascending=True)
    snp_a = str(epi.iloc[0]["SNP_A"])
    snp_b = str(epi.iloc[0]["SNP_B"])
    a = pd.to_numeric(df[snp_a], errors="coerce").astype(np.float32).fillna(df[snp_a].astype(float).mean())
    b = pd.to_numeric(df[snp_b], errors="coerce").astype(np.float32).fillna(df[snp_b].astype(float).mean())
    ab = (a * b).values
    X = np.vstack([a.values, b.values, ab]).T
    return X, (snp_a, snp_b)

def auc_from_logit(X, y):
    clf = LogisticRegression(max_iter=2000, solver="liblinear")
    clf.fit(X, y)
    proba = clf.predict_proba(X)[:,1]
    return roc_auc_score(y, proba), average_precision_score(y, proba)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--pheno", default="case")
    ap.add_argument("--fighi_results", required=True)
    ap.add_argument("--plink_epi", required=True)
    ap.add_argument("--mdr_pairs", required=True)
    ap.add_argument("--out_json", required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--top_k_fighi", type=int, default=50)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    y = df[args.pheno].astype(int).values

    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=args.seed)
    (tr_idx, te_idx) = next(sss.split(df, y))
    df_tr, df_te = df.iloc[tr_idx], df.iloc[te_idx]
    y_tr, y_te = y[tr_idx], y[te_idx]

    out = {}

    # ---- FIGHI
    X_tr, feats = design_from_fighi(df_tr, args.fighi_results, top_k=args.top_k_fighi)
    if X_tr is not None:
        # rebuild on test
        X_te_list = []
        for cols in feats:
            X_te_list.append(build_feature_from_tuple(df_te, list(cols)))
        X_te = np.vstack(X_te_list).T
        auc_tr, ap_tr = auc_from_logit(X_tr, y_tr)
        auc_te, ap_te = auc_from_logit(X_te, y_te)
        out["FIGHI"] = {"features": [list(c) for c in feats],
                        "AUC_train": float(auc_tr), "AUC_test": float(auc_te),
                        "AP_train": float(ap_tr), "AP_test": float(ap_te)}
    else:
        out["FIGHI"] = {"error": "no features parsed"}

    # ---- PLINK
    X_tr, pair = design_from_plink_pair(df_tr, args.plink_epi)
    if X_tr is not None:
        # mirror on test
        a = pd.to_numeric(df_te[pair[0]], errors="coerce").astype(np.float32).fillna(df_te[pair[0]].astype(float).mean())
        b = pd.to_numeric(df_te[pair[1]], errors="coerce").astype(np.float32).fillna(df_te[pair[1]].astype(float).mean())
        ab = (a*b).values
        X_te = np.vstack([a.values, b.values, ab]).T
        auc_tr, ap_tr = auc_from_logit(X_tr, y_tr)
        auc_te, ap_te = auc_from_logit(X_te, y_te)
        out["PLINK"] = {"pair": list(pair),
                        "AUC_train": float(auc_tr), "AUC_test": float(auc_te),
                        "AP_train": float(ap_tr), "AP_test": float(ap_te)}
    else:
        out["PLINK"] = {"error": "no epi file"}

    # ---- MDR
    if os.path.exists(args.mdr_pairs):
        m = pd.read_csv(args.mdr_pairs)
        if len(m):
            snp1, snp2 = m.iloc[0]["SNP1"], m.iloc[0]["SNP2"]
            Xi_tr = df_tr[[snp1, snp2]].apply(pd.to_numeric, errors="coerce").fillna(df_tr[[snp1, s*np2]].mean()).values
            Xi_te = df_te[[snp1, snp2]].apply(pd.to_numeric, errors="coerce").fillna(df_te[[snp1, s*np2]].mean()).values
            # Fit MDR on train, predict on test
            mdl = MDR()
            mdl.fit(Xi_tr, y_tr)
            yhat_tr = mdl.predict(Xi_tr)
            yhat_te = mdl.predict(Xi_te)
            out["MDR"] = {"pair": [snp1, snp2],
                          "AUC_train": float(roc_auc_score(y_tr, yhat_tr)),
                          "AUC_test": float(roc_auc_score(y_te, yhat_te))}
        else:
            out["MDR"] = {"error": "empty mdr_pairs.csv"}
    else:
        out["MDR"] = {"error": "missing mdr_pairs.csv"}

    with open(args.out_json, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
