import numpy as np
import pandas as pd
from itertools import combinations
from typing import Dict, List, Tuple, Optional


import numpy as np
import pandas as pd
from itertools import combinations
from typing import Dict, List, Tuple, Optional

# --- robust imports: module or script mode ---
try:
    # package mode (python -m fighi_ext.run_cli)
    from .utils import standardize, maf, product_feature, bh_fdr
    from .apriori import apriori_generate, prune_by_subsets, prune_by_variance
    from .fisher import score_test_fi_gain_logistic, score_test_fi_gain_linear, batch_fi_gain_logistic
    from .adaptive import information_ratio, planner_max_K
    from .hypergraph import export_csv, export_gml, export_hyper_json, export_cytoscape_cyjs
except Exception:
    # script mode (python pipeline.py or run_cli.py in same folder)
    import os, sys
    sys.path.append(os.path.dirname(__file__))
    from utils import standardize, maf, product_feature, bh_fdr
    from apriori import apriori_generate, prune_by_subsets, prune_by_variance
    from fisher import score_test_fi_gain_logistic, score_test_fi_gain_linear, batch_fi_gain_logistic
    from adaptive import information_ratio, planner_max_K
    from hypergraph import export_csv, export_gml, export_hyper_json, export_cytoscape_cyjs
# ------------------------------------------------



#from .utils import standardize, maf, product_feature, bh_fdr
#from .apriori import apriori_generate, prune_by_subsets, prune_by_variance
#from .fisher import score_test_fi_gain_logistic, score_test_fi_gain_linear, batch_fi_gain_logistic
#from .adaptive import information_ratio, planner_max_K
#from .hypergraph import export_csv, export_gml, export_hyper_json, export_cytoscape_cyjs
#from .hypergraph import export_csv, export_gml, export_hyper_json, export_cytoscape_cyjs, export_graphml

class FIGHIPipeline:
    def __init__(self, trait_type="binary", min_maf=0.05, min_var=1e-6, n_perm=200, rng=0):
        assert trait_type in ("binary", "linear")
        self.trait_type = trait_type
        self.min_maf = min_maf
        self.min_var = min_var
        self.n_perm = n_perm
        self.rng = rng
        self.N = None

    def screen(self, X: pd.DataFrame, y: np.ndarray, top_m: int=200) -> List[str]:
        Z = standardize(X)
        # simple correlation or point-biserial for binary
        scores = {c: abs(np.corrcoef(Z[c].values, y)[0,1]) for c in Z.columns}
        ranked = sorted(scores.items(), key=lambda kv: -np.nan_to_num(kv[1]))
        kept = []
        for c,_ in ranked:
            if maf(X[c]) >= self.min_maf:
                kept.append(c)
            if len(kept) >= top_m:
                break
        return kept

    def generate_candidates(self, atoms: List[str], prev_level: List[Tuple[str,...]], k: int):
        if k == 2:
            cands = [tuple(sorted(t)) for t in combinations(atoms, 2)]
        else:
            cands = apriori_generate(prev_level)
            cands = prune_by_subsets(cands, set(prev_level))
        return cands

    def evaluate_candidates(self, Xz: pd.DataFrame, y: np.ndarray, cands: List[Tuple[str,...]]):
        results = []
        X_base = np.zeros((len(y),0))
        for S in cands:
            x_new = product_feature(Xz, S)
            if np.var(x_new) < self.min_var:
                continue
            if self.trait_type == "binary":
                fi_gain, bh, info = score_test_fi_gain_logistic(y, X_base, x_new)
            else:
                fi_gain, bh, info = score_test_fi_gain_linear(y, X_base, x_new)
            results.append(dict(S=S, fi_gain=fi_gain, beta_hat=bh, info=info))
        return sorted(results, key=lambda d: -d["fi_gain"])

    def adaptive_search(self, X: pd.DataFrame, y: np.ndarray, max_order: int=5,
                        lambda_by_k: Optional[Dict[int,float]]=None,
                        eps_ratio: float=0.95, target_OR: float=1.3):
        self.N = len(y)
        Xz = standardize(X.copy())
        atoms = self.screen(X, y)
        K_feas, I_req = planner_max_K(N=len(y), target_OR=target_OR,
                                      maf_low=0.1, maf_high=0.4, n_tuples=len(atoms))
        results_all = []
        retained = {}
        cumI = {}
        for k in range(2, min(max_order, max(2,K_feas))+1):
            cands = self.generate_candidates(atoms, retained.get(k-1, []), k)
            cands = prune_by_variance(Xz, cands, self.min_var)
            res = self.evaluate_candidates(Xz, y, cands)
            lam = (lambda_by_k or {}).get(k, 0.0)
            kept = [r for r in res if r["fi_gain"] > lam]
            retained[k] = [r["S"] for r in kept]
            cumI[k] = sum(r["info"] for r in kept)
            for r in kept:
                r["order"] = k
            results_all.extend(kept)
            if k>=3 and information_ratio(cumI, k-1) > eps_ratio:
                break
        return dict(results=results_all, cumI=cumI, retained=retained, K_feasible=K_feas, I_req=I_req)

    def export(self, outdir: str, search_result: Dict):
        import os
        os.makedirs(outdir, exist_ok=True)
        rows = []
        nodes = []
        node_index = {}
        for r in search_result["results"]:
            S = tuple(r["S"])
            for s in S:
                if s not in node_index:
                    node_index[s] = len(nodes)
                    nodes.append(s)
            rows.append(["|".join(S), r["order"], r["fi_gain"], r.get("pval", np.nan), r["beta_hat"], r["info"]])
        export_csv(os.path.join(outdir,"fighi_results.csv"),
                   rows, header=["hyperedge","order","fi_gain","pval","beta_hat","info"])
        hyperedges = []
        for idx, r in enumerate(search_result["results"]):
            members = [node_index[s] for s in r["S"]]
            hyperedges.append((idx, members, float(r["fi_gain"])))
        export_graphml(os.path.join(outdir,"fighi_hypergraph.graphml"), nodes, hyperedges)
        return {
            "csv": os.path.join(outdir,"fighi_results.csv"),
            "gml": os.path.join(outdir,"fighi_hypergraph.gml"),
            "graphml": os.path.join(outdir,"fighi_hypergraph.graphml"),
            "hyper": os.path.join(outdir,"fighi_hypergraph.hyper"),
            "cyjs": os.path.join(outdir,"fighi_cytoscape.cyjs"),
        }
        #export_gml(os.path.join(outdir,"fighi_hypergraph.gml"), nodes, hyperedges)
        #export_hyper_json(os.path.join(outdir,"fighi_hypergraph.hyper"), nodes, hyperedges)
        #export_cytoscape_cyjs(os.path.join(outdir,"fighi_cytoscape.cyjs"), nodes, hyperedges)
        #return {"csv": os.path.join(outdir,"fighi_results.csv"),
        #        "gml": os.path.join(outdir,"fighi_hypergraph.gml"),
        #        "hyper": os.path.join(outdir,"fighi_hypergraph.hyper"),
        #        "cyjs": os.path.join(outdir,"fighi_cytoscape.cyjs")}
