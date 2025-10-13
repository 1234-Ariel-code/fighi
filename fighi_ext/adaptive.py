from typing import Dict

def information_ratio(cum_I_by_order: Dict[int, float], K: int) -> float:
    num = sum(v for k,v in cum_I_by_order.items() if k<=K)
    den = sum(v for k,v in cum_I_by_order.items() if k<=K+1)
    if den <= 0:
        return 1.0
    return num/den

def planner_max_K(N: int, target_OR: float, maf_low: float, maf_high: float, n_tuples: int,
                  alpha=0.05, power=0.8):
    from math import log
    from scipy.stats import norm
    z = norm.ppf(1-alpha/2.0) + norm.ppf(power)
    beta = log(max(target_OR, 1.0))
    I_req = (z/beta)**2 if beta>0 else float("inf")
    p = (maf_low+maf_high)/2.0
    per_tuple_I = (p*(1-p))
    K = 1
    while K < 10 and n_tuples * (per_tuple_I**K) * N >= I_req:
        K += 1
    return K-1, I_req
