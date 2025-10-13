# fighi_ext/report.py
import os, json, pickle, numpy as np, pandas as pd, matplotlib.pyplot as plt
from typing import Dict, List

def _maf(x: np.ndarray) -> float:
    x = x[~np.isnan(x)]
    if x.size == 0: return 0.0
    p = np.clip(x, 0, 2).mean()/2.0
    return float(min(p, 1-p))

def write_feature_scores(outdir: str, results: List[Dict], X: pd.DataFrame):
    os.makedirs(outdir, exist_ok=True)
    # aggregate FI by SNP
    fi_main = {c:0.0 for c in X.columns}
    fi_inter = {c:0.0 for c in X.columns}
    for r in results:
        S, gain = r["S"], float(r["fi_gain"])
        if len(S) == 1:
            fi_main[S[0]] += gain
        else:
            for s in S:
                fi_inter[s] += gain/len(S)
    rows=[]
    for c in X.columns:
        rows.append([c, fi_main[c]+fi_inter[c], fi_main[c], fi_inter[c], _maf(X[c].values)])
    df = pd.DataFrame(rows, columns=["SNP","FI_total","FI_main","FI_interact","MAF"])
    df = df.sort_values("FI_total", ascending=False).reset_index(drop=True)
    df["Rank"] = np.arange(1, len(df)+1)
    df["Gene"] = ""
    df["Pathway"] = ""
    path = os.path.join(outdir, "fighi_feature_scores.csv")
    df.to_csv(path, index=False)
    return path

def write_summary(outdir: str, search_result: Dict, X: pd.DataFrame, y: np.ndarray, max_order: int, version="0.3.2"):
    res = search_result["results"]
    mean_fi = float(np.mean([r["fi_gain"] for r in res])) if res else 0.0
    info = {
      "n_samples": int(len(y)),
      "n_snps": int(X.shape[1]),
      "phenotype": "case",
      "mean_FI_gain": mean_fi,
      "max_order": int(max_order),
      "n_significant_interactions": int(len(res)),
      "mean_stability": None,
      "runtime_minutes": None,
      "version": version
    }
    with open(os.path.join(outdir, "fighi_summary.json"), "w") as f:
        json.dump(info, f, indent=2)
    return info

def save_model(outdir: str, search_result: Dict):
    path = os.path.join(outdir, "fighi_model.pkl")
    with open(path, "wb") as f:
        pickle.dump(search_result, f, protocol=4)
    return path

def write_log(outdir: str, X: pd.DataFrame, y: np.ndarray, search_result: Dict):
    path = os.path.join(outdir, "fighi_log.txt")
    with open(path, "w") as f:
        f.write("[INFO] Input summary\n")
        f.write(f"  samples={len(y)}, snps={X.shape[1]}\n")
        kept = search_result.get("retained", {})
        f.write("[INFO] Interactions kept per order:\n")
        for k in sorted(kept):
            f.write(f"  K={k}: {len(kept[k])}\n")
        f.write("\nTop interactions:\n")
        for r in sorted(search_result["results"], key=lambda d: -d["fi_gain"])[:25]:
            f.write(f"  K={r['order']}  S={r['S']}  FI={r['fi_gain']:.6g}\n")
    return path

def _hist(ax, vals, title, xlabel):
    ax.hist(vals, bins=30)
    ax.set_title(title); ax.set_xlabel(xlabel); ax.set_ylabel("Count")

def _matrix_from_pairs(snps, results):
    # simple heatmap: pairwise FI sum
    idx = {s:i for i,s in enumerate(snps)}
    M = np.zeros((len(snps), len(snps)), float)
    for r in results:
        S, g = r["S"], float(r["fi_gain"])
        if len(S)==2:
            i, j = idx[S[0]], idx[S[1]]
            M[i,j]+=g; M[j,i]+=g
    return M

def make_plots(outdir: str, search_result: Dict, X: pd.DataFrame):
    os.makedirs(os.path.join(outdir,"plots"), exist_ok=True)
    gains = [float(r["fi_gain"]) for r in search_result["results"]]
    # 1. FI gain hist
    fig, ax = plt.subplots(figsize=(6,4))
    _hist(ax, gains, "FI Gain Distribution", "FI Gain")
    fig.savefig(os.path.join(outdir,"plots","FI_gain_distribution.png"), dpi=200); plt.close(fig)
    # 2. Interaction heatmap (pairwise)
    snps = list(X.columns[:150])  # keep readable
    M = _matrix_from_pairs(snps, search_result["results"])
    fig, ax = plt.subplots(figsize=(7,6))
    im = ax.imshow(M[:len(snps),:len(snps)], aspect="auto")
    fig.colorbar(im, ax=ax)
    ax.set_title("Interaction heatmap (pairwise FI)"); ax.set_xlabel("SNP"); ax.set_ylabel("SNP")
    fig.savefig(os.path.join(outdir,"plots","interaction_heatmap.png"), dpi=200); plt.close(fig)
    # 3. FI vs order
    byK = {}
    for r in search_result["results"]:
        byK.setdefault(r["order"], []).append(float(r["fi_gain"]))
    K = sorted(byK)
    means = [np.mean(byK[k]) for k in K] if K else []
    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(K, means, marker="o"); ax.set_xlabel("Order K"); ax.set_ylabel("Mean FI"); ax.set_title("Fisher Information vs Order")
    fig.savefig(os.path.join(outdir,"plots","fisher_information_vs_order.png"), dpi=200); plt.close(fig)
    # 4. Convergence proxy (cumulative info)
    cumI = search_result.get("cumI", {})
    Kc = sorted(cumI)
    vals = [cumI[k] for k in Kc]
    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(Kc, vals, marker="o")
    ax.set_xlabel("Order K"); ax.set_ylabel("Cumulative info"); ax.set_title("Convergence diagnostics")
    fig.savefig(os.path.join(outdir,"plots","convergence_diagnostics.png"), dpi=200); plt.close(fig)

