import numpy as np

class LogisticGLM:
    """Minimal IRLS logistic regression with weights and info matrix."""
    def __init__(self, max_iter=100, tol=1e-8):
        self.max_iter = max_iter
        self.tol = tol
        self.coef_ = None
        self.fitted_ = False

    @staticmethod
    def sigmoid(eta):
        return 1.0 / (1.0 + np.exp(-eta))

    def fit(self, X, y, beta0=None):
        n, p = X.shape
        beta = np.zeros(p) if beta0 is None else beta0.copy()
        for _ in range(self.max_iter):
            eta = X @ beta
            p = self.sigmoid(eta)
            W = p * (1.0 - p)
            W = np.clip(W, 1e-9, None)
            g = X.T @ (y - p)
            XTW = X.T * W
            H = XTW @ X
            try:
                step = np.linalg.solve(H, g)
            except np.linalg.LinAlgError:
                H = H + 1e-6 * np.eye(H.shape[0])
                step = np.linalg.solve(H, g)
            beta_new = beta + step
            if np.linalg.norm(step) < self.tol:
                beta = beta_new
                break
            beta = beta_new
        self.coef_ = beta
        self.W_ = W
        self.p_ = p
        self.H_ = H
        self.fitted_ = True
        return self

    def observed_information(self, X):
        if not self.fitted_:
            raise RuntimeError("Fit the model first.")
        p = self.sigmoid(X @ self.coef_)
        W = np.clip(p * (1 - p), 1e-9, None)
        XTW = X.T * W
        return XTW @ X, W, p

class LinearGLM:
    """Gaussian GLM with identity link (ordinary least squares)."""
    def __init__(self):
        self.fitted_ = False

    def fit(self, X, y):
        XtX = X.T @ X
        Xty = X.T @ y
        try:
            beta = np.linalg.solve(XtX, Xty)
        except np.linalg.LinAlgError:
            beta = np.linalg.lstsq(X, y, rcond=None)[0]
        self.coef_ = beta
        resid = y - X @ beta
        self.sigma2_ = float((resid @ resid) / max(1, X.shape[0] - X.shape[1]))
        self.H_ = XtX / self.sigma2_
        self.fitted_ = True
        return self

    def observed_information(self, X):
        if not self.fitted_:
            raise RuntimeError("Fit the model first.")
        return (X.T @ X) / self.sigma2_, np.ones(X.shape[0]) / self.sigma2_, None
