import numpy as np

# --- robust import for GLMs ---
try:
    from .glm import LogisticGLM, LinearGLM
except Exception:
    import os, sys
    sys.path.append(os.path.dirname(__file__))
    from glm import LogisticGLM, LinearGLM
# --------------------------------


#from .glm import LogisticGLM, LinearGLM

def score_test_fi_gain_logistic(y, X_base, x_new):
    n = y.shape[0]
    if X_base.size == 0:
        eta = np.zeros(n)
    else:
        base = LogisticGLM().fit(X_base, y)
        eta = X_base @ base.coef_
    p = 1/(1+np.exp(-eta))
    W = np.clip(p*(1-p), 1e-9, None)
    u = float(np.dot(x_new, (y - p)))
    i_ss = float(np.dot(x_new * W, x_new))
    beta_1step = u / max(i_ss, 1e-12)
    fi_gain = 0.5 * (beta_1step**2) * i_ss
    return float(fi_gain), float(beta_1step), float(i_ss)

def score_test_fi_gain_linear(y, X_base, x_new):
    n = y.shape[0]
    if X_base.size == 0:
        resid = y - y.mean()
        sigma2 = float(np.var(resid, ddof=1))
    else:
        base = LinearGLM().fit(X_base, y)
        resid = y - X_base @ base.coef_
        sigma2 = base.sigma2_
    u = float(np.dot(x_new, resid) / sigma2)
    i_ss = float(np.dot(x_new, x_new) / sigma2)
    beta_1step = u / max(i_ss, 1e-12)
    fi_gain = 0.5 * (beta_1step**2) * i_ss
    return float(fi_gain), float(beta_1step), float(i_ss)

def batch_fi_gain_logistic(y, X_base, X_new_mat, use_gpu=False):
    """Compute FI-Gain for many candidates at once; optional cupy/numba speedup."""
    n, C = X_new_mat.shape
    if X_base.size == 0:
        eta = np.zeros(n)
    else:
        base = LogisticGLM().fit(X_base, y)
        eta = X_base @ base.coef_
    p = 1/(1+np.exp(-eta))
    W = np.clip(p*(1-p), 1e-9, None)
    # CPU vectorized
    y_minus_p = (y - p).reshape(-1,1)
    U = (X_new_mat * y_minus_p).sum(axis=0)
    I = (X_new_mat * (W.reshape(-1,1)) * X_new_mat).sum(axis=0)
    beta = U / np.maximum(I, 1e-12)
    fi_gain = 0.5 * (beta**2) * I
    return fi_gain, beta, I
