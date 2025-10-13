#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate fighi_feature_scores.csv with Gene and Pathway.

Priority:
  1) SNP→Gene via cS2G (combined_cS2G.tsv or *.SGscore.gz directory)
  2) Else user-supplied map CSV/TSV with cols: SNP,Gene (or SNP,GENE,rsid,gene_symbol)
  3) Pathway via g:Profiler (online) OR local GMTs (offline)

Outputs (written next to --feature_csv unless --outdir given):
  - fighi_feature_scores_annotated.csv
  - fighi_pathway_enrichment.csv
  - annotate_log.txt
"""

import os, sys, csv, gzip, json, time, math, argparse, shutil
from typing import Dict, List, Tuple, Optional, Iterable
from collections import defaultdict, Counter

import pandas as pd
import numpy as np

# ======== small, memory-safe utils ========
def _open_gz_or_txt(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode, newline="")

def _write_log(outdir: str, lines: Iterable[str]):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "annotate_log.txt"), "w") as f:
        for ln in lines:
            f.write(ln.rstrip() + "\n")

def _find_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    low = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None

def _load_feature_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    need = _find_col(df, ["SNP"])
    if need is None:
        raise ValueError("feature CSV must contain column 'SNP'")
    return df

def _read_gmt(path: str) -> Dict[str, set]:
    d = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3: continue
            d[parts[0]] = set(parts[2:])
    return d

# ======== cS2G loaders (OOM-safe; streaming) ========
def _load_cs2g_combined(path: str, min_score: float = 0.05) -> Dict[str, str]:
    """
    Expect a file like combined_cS2G.tsv with columns including 'rsid' and 'symbol' (or 'gene_symbol')
    and 'score' (or 'cS2G').
    Keep the *best* gene per rsid above min_score.
    """
    best = {}         # rsid -> (score, gene)
    with _open_gz_or_txt(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        # tolerant names
        rs_col  = next((c for c in ["rsid", "SNP", "snp"] if c in idx), None)
        g_col   = next((c for c in ["symbol", "gene_symbol", "GENE", "gene"] if c in idx), None)
        sc_col  = next((c for c in ["score","cS2G","SGscore"] if c in idx), None)
        if rs_col is None or g_col is None or sc_col is None:
            raise ValueError("combined_cS2G.tsv must contain rsid + symbol + score-like columns.")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rs, gene = parts[idx[rs_col]], parts[idx[g_col]]
            try:
                s = float(parts[idx[sc_col]])
            except Exception:
                continue
            if s < min_score: 
                continue
            if (rs not in best) or (s > best[rs][0]):
                best[rs] = (s, gene)
    return {rs: g for rs, (s, g) in best.items()}

def _load_cs2g_dir(cs2g_dir: str, min_score: float = 0.05) -> Dict[str, str]:
    """
    Scan *.SGscore.gz files OR any *.tsv/*.gz containing columns rsid+symbol+score.
    Keep best gene per SNP (highest score ≥ min_score).
    """
    best = {}
    files = [os.path.join(cs2g_dir, f) for f in os.listdir(cs2g_dir)
             if f.endswith(".gz") or f.endswith(".tsv") or f.endswith(".txt")]
    for fpath in files:
        try:
            with _open_gz_or_txt(fpath, "rt") as f:
                header = f.readline().rstrip("\n").split("\t")
                idx = {h: i for i, h in enumerate(header)}
                rs_col  = next((c for c in ["rsid","SNP","snp"] if c in idx), None)
                g_col   = next((c for c in ["symbol","gene_symbol","GENE","gene"] if c in idx), None)
                sc_col  = next((c for c in ["score","SGscore","cS2G"] if c in idx), None)
                if rs_col is None or g_col is None or sc_col is None:
                    continue
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    rs, gene = parts[idx[rs_col]], parts[idx[g_col]]
                    try:
                        s = float(parts[idx[sc_col]])
                    except Exception:
                        continue
                    if s < min_score: continue
                    if (rs not in best) or (s > best[rs][0]):
                        best[rs] = (s, gene)
        except Exception:
            continue
    return {rs: g for rs,(s,g) in best.items()}

# ======== user mapping loader ========
def _load_user_map(path: str) -> Dict[str, str]:
    """
    Accept CSV/TSV with columns like: SNP,Gene OR rsid,gene_symbol, etc.
    """
    sep = "\t" if path.endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep)
    rs_col = _find_col(df, ["SNP","rsid","RSID"])
    g_col  = _find_col(df, ["Gene","GENE","gene","symbol","SYMBOL","gene_symbol"])
    if rs_col is None or g_col is None:
        raise ValueError("Mapping file must contain SNP/rsid and Gene/gene_symbol columns.")
    mp = {}
    for rs, g in zip(df[rs_col].astype(str), df[g_col].astype(str)):
        if rs not in mp:
            mp[rs] = g
    return mp

# ======== pathway assignment ========
def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    if pvals.size == 0: return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    m = float(len(pvals))
    q = np.empty_like(ranked)
    min_q = 1.0
    for i in range(len(ranked)-1, -1, -1):
        min_q = min(min_q, ranked[i] * m / (i+1))
        q[i] = min_q
    out = np.empty_like(q)
    out[order] = q
    return out

def _enrich_gprofiler(genes: List[str], organism="hsapiens") -> Tuple[pd.DataFrame, Dict[str, List[str]]]:
    """
    Online enrichment with g:Profiler (robust import + chunking).
    Returns (table, term->gene-list dict)
    """
    try:
        import importlib, subprocess, sys as _sys
        try:
            import gprofiler
        except Exception:
            subprocess.check_call([_sys.executable, "-m", "pip", "install", "-q", "gprofiler-official"])
        from gprofiler import GProfiler
    except Exception as e:
        raise RuntimeError(f"g:Profiler is unavailable ({e}). Use --gmt_files instead.")

    gp = GProfiler(return_dataframe=True, user_agent="FIGHI")
    # do in one go; gprofiler handles up to ~5–10k list sizes fine
    res = gp.profile(organism=organism, query=genes)
    if res is None or res.empty:
        return pd.DataFrame(), {}
    # build term->gene membership
    term_to_genes = {}
    for term, glist in zip(res["name"], res["intersection"]):
        gset = [g.strip() for g in str(glist).split(",")]
        term_to_genes[term] = gset
    # Add our own FDR if not present
    if "p_value" in res.columns and "p_value" not in res.columns:
        res["FDR"] = _bh_fdr(res["p_value"].values)
    elif "p_value" in res.columns:
        res["FDR"] = _bh_fdr(res["p_value"].values)
    elif "p_value" not in res.columns and "p_value" in res.columns:
        res["FDR"] = res["p_value"]
    else:
        # gprofiler uses 'p_value' column; keep defensive
        res["FDR"] = _bh_fdr(res["p_value"].values) if "p_value" in res.columns else np.nan
    return res, term_to_genes

def _enrich_local_gmts(genes: List[str], gmt_files: List[str]) -> Tuple[pd.DataFrame, Dict[str, List[str]]]:
    """
    Offline hypergeometric using local GMT files.
    """
    from math import comb
    from scipy.stats import hypergeom

    # union universe
    db = {}
    for gmt in gmt_files:
        if not os.path.exists(gmt): continue
        for k,v in _read_gmt(gmt).items():
            db[k] = set(v)
    if not db:
        return pd.DataFrame(), {}

    universe = set().union(*db.values())
    genes = [g for g in genes if g in universe]
    if not genes:
        return pd.DataFrame(), {}
    N = len(universe)
    n = len(genes)
    gset = set(genes)

    rows = []
    term_to_genes = {}
    for term, members in db.items():
        K = len(members)
        if K == 0: continue
        k = len(gset & members)
        if k == 0: continue
        p = hypergeom.sf(k-1, N, K, n)
        rows.append((term, p, k, K, n, N))
        term_to_genes[term] = sorted(list(gset & members))

    if not rows:
        return pd.DataFrame(), {}
    res = pd.DataFrame(rows, columns=["term","p_value","k","K","n","N"]).sort_values("p_value")
    res["FDR"] = _bh_fdr(res["p_value"].values)
    return res, term_to_genes

def _assign_gene_to_pathway(genes: List[str], enrichment_df: pd.DataFrame, term_to_genes: Dict[str, List[str]]) -> Dict[str, str]:
    """
    For each gene, pick the best (lowest FDR) term that contains it.
    """
    assign = {}
    if enrichment_df is None or enrichment_df.empty:
        return assign
    # ensure we have an FDR column
    fdr_col = "FDR" if "FDR" in enrichment_df.columns else None
    if fdr_col is None and "p_value" in enrichment_df.columns:
        enrichment_df = enrichment_df.copy()
        enrichment_df["FDR"] = _bh_fdr(enrichment_df["p_value"].values)
        fdr_col = "FDR"
    if fdr_col is None:
        return assign

    enrichment_df = enrichment_df.sort_values(fdr_col)
    # prefer most significant terms
    for _, row in enrichment_df.iterrows():
        term = row.get("term", row.get("name", "term"))
        members = set(term_to_genes.get(term, []))
        for g in genes:
            if g in members and g not in assign:
                assign[g] = term
    return assign

# ======== main ========
def main():
    ap = argparse.ArgumentParser(description="Annotate FIGHI features with Gene & Pathway.")
    ap.add_argument("--feature_csv", required=True, help="fighi_feature_scores.csv")
    ap.add_argument("--outdir", help="Output directory (default: feature_csv dir)")
    # SNP->Gene sources (highest available wins)
    ap.add_argument("--cs2g_combined", help="Path to combined_cS2G.tsv")
    ap.add_argument("--cs2g_dir", help="Directory containing cS2G *.SGscore.gz files")
    ap.add_argument("--snp_gene_map", help="User mapping CSV/TSV with columns SNP, Gene (or synonyms)")
    ap.add_argument("--min_cs2g", type=float, default=0.05, help="Min cS2G score to accept")
    # Pathways
    ap.add_argument("--use_gprofiler", action="store_true", help="Use g:Profiler online enrichment")
    ap.add_argument("--gmt_files", nargs="*", help="Local GMT files for offline enrichment")
    ap.add_argument("--organism", default="hsapiens")
    args = ap.parse_args()

    outdir = args.outdir or os.path.dirname(os.path.abspath(args.feature_csv)) or "."
    os.makedirs(outdir, exist_ok=True)
    logs = []

    # 1) Load features
    df = _load_feature_table(args.feature_csv)
    snps = df["SNP"].astype(str).tolist()
    logs.append(f"[INFO] Loaded {len(snps)} SNPs from {args.feature_csv}")

    # 2) SNP→Gene mapping
    rs2gene: Dict[str,str] = {}
    if args.cs2g_combined and os.path.exists(args.cs2g_combined):
        logs.append(f"[INFO] Using cS2G combined: {args.cs2g_combined}")
        rs2gene = _load_cs2g_combined(args.cs2g_combined, min_score=args.min_cs2g)
    elif args.cs2g_dir and os.path.isdir(args.cs2g_dir):
        logs.append(f"[INFO] Using cS2G directory: {args.cs2g_dir}")
        rs2gene = _load_cs2g_dir(args.cs2g_dir, min_score=args.min_cs2g)
    elif args.snp_gene_map and os.path.exists(args.snp_gene_map):
        logs.append(f"[INFO] Using user map: {args.snp_gene_map}")
        rs2gene = _load_user_map(args.snp_gene_map)
    else:
        logs.append("[WARN] No cS2G or user map provided; Gene column will remain empty.")

    # attach genes
    genes = []
    for rs in snps:
        g = rs2gene.get(rs, "")
        genes.append(g)
    df["Gene"] = genes

    # 3) Pathway enrichment + assignment
    unique_genes = sorted(set([g for g in genes if g]))
    assign_g2p = {}
    enr_table = pd.DataFrame()

    if unique_genes:
        if args.use_gprofiler:
            try:
                enr_table, term_to_genes = _enrich_gprofiler(unique_genes, organism=args.organism)
            except Exception as e:
                logs.append(f"[ERROR] g:Profiler failed ({e}); trying GMTs if provided.")
                enr_table, term_to_genes = (pd.DataFrame(), {})
        if (enr_table is None or enr_table.empty) and args.gmt_files:
            enr_table, term_to_genes = _enrich_local_gmts(unique_genes, args.gmt_files)

        if enr_table is not None and not enr_table.empty:
            # normalize expected columns
            if "term" not in enr_table.columns and "name" in enr_table.columns:
                enr_table = enr_table.rename(columns={"name":"term"})
            # write enrichment table
            enr_out = os.path.join(outdir, "fighi_pathway_enrichment.csv")
            enr_table.to_csv(enr_out, index=False)
            logs.append(f"[INFO] Pathway enrichment written: {enr_out}")
            # assign per gene
            assign_g2p = _assign_gene_to_pathway(unique_genes, enr_table, term_to_genes)
        else:
            logs.append("[WARN] No pathway enrichment available; Pathway column will be blank.")
    else:
        logs.append("[WARN] No genes mapped; skipping enrichment.")

    # Fill per SNP pathway (use gene’s assigned pathway if present)
    df["Pathway"] = [assign_g2p.get(g, "") if g else "" for g in df["Gene"].astype(str)]

    # 4) Save annotated features
    anno_path = os.path.join(outdir, "fighi_feature_scores_annotated.csv")
    df.to_csv(anno_path, index=False)
    logs.append(f"[INFO] Annotated features written: {anno_path}")

    # 5) Write log
    _write_log(outdir, logs)
    print("\n".join(logs))

if __name__ == "__main__":
    main()

