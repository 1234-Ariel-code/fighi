from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import combinations

import numpy as np
from scipy.special import expit
from scipy.stats import chi2


@dataclass(slots=True)
class ScoreResult:
    statistic: float
    fi_gain: float
    beta_hat: float
    information: float
    p_value: float
    null_converged: bool
    null_iterations: int

    def to_dict(self) -> dict:
        return asdict(self)


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, with NaNs preserved."""
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    values = np.clip(p[valid], 0.0, 1.0)
    if not values.size:
        return out
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out[valid] = restored
    return out


def adjust_pvalues(p_values: np.ndarray, method: str) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    if method == "fdr_bh":
        return bh_adjust(p)
    if method == "bonferroni":
        valid_n = max(1, int(np.isfinite(p).sum()))
        return np.where(np.isfinite(p), np.minimum(1.0, p * valid_n), np.nan)
    if method == "none":
        return p.copy()
    raise ValueError(f"Unknown correction: {method}")


def _logistic_loglik(y: np.ndarray, eta: np.ndarray) -> float:
    return float(np.sum(y * -np.logaddexp(0.0, -eta) + (1.0 - y) * -np.logaddexp(0.0, eta)))


def fit_logistic_null(
    design: np.ndarray,
    y: np.ndarray,
    max_iter: int = 100,
    tolerance: float = 1e-8,
    ridge: float = 1e-8,
) -> tuple[np.ndarray, np.ndarray, bool, int]:
    """Fit a small nuisance logistic model with damped Newton updates."""
    n_columns = design.shape[1]
    beta = np.zeros(n_columns, dtype=float)
    prevalence = float(np.clip(y.mean(), 1e-6, 1 - 1e-6))
    if n_columns and np.allclose(design[:, 0], 1.0):
        beta[0] = np.log(prevalence / (1.0 - prevalence))
    penalty = np.full(n_columns, ridge, dtype=float)
    if n_columns and np.allclose(design[:, 0], 1.0):
        penalty[0] = 0.0

    converged = False
    iterations = 0
    for iterations in range(1, max_iter + 1):
        eta = np.clip(design @ beta, -35.0, 35.0)
        probability = expit(eta)
        weights = np.clip(probability * (1.0 - probability), 1e-9, None)
        score = design.T @ (y - probability) - penalty * beta
        information = (design.T * weights) @ design + np.diag(penalty)
        try:
            step = np.linalg.solve(information, score)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(information, score, rcond=None)[0]

        current_ll = _logistic_loglik(y, eta) - 0.5 * float(np.dot(penalty * beta, beta))
        scale = 1.0
        while scale > 1e-6:
            proposal = beta + scale * step
            proposal_eta = np.clip(design @ proposal, -35.0, 35.0)
            proposal_ll = _logistic_loglik(y, proposal_eta) - 0.5 * float(
                np.dot(penalty * proposal, proposal)
            )
            if proposal_ll >= current_ll - 1e-10:
                break
            scale *= 0.5
        beta = beta + scale * step
        if np.linalg.norm(scale * step, ord=np.inf) <= tolerance:
            converged = True
            break

    probability = expit(np.clip(design @ beta, -35.0, 35.0))
    return beta, probability, converged, iterations


def _lower_order_design(candidate_matrix: np.ndarray, covariates: np.ndarray | None) -> np.ndarray:
    n_samples, order = candidate_matrix.shape
    columns = [np.ones(n_samples, dtype=float)]
    if covariates is not None and covariates.size:
        columns.extend(covariates[:, index] for index in range(covariates.shape[1]))
    for subset_order in range(1, order):
        for subset in combinations(range(order), subset_order):
            columns.append(np.prod(candidate_matrix[:, subset], axis=1))
    return np.column_stack(columns)


def _efficient_information(z: np.ndarray, design: np.ndarray, weights: np.ndarray) -> float:
    weighted_design = design * weights[:, None]
    nuisance_information = design.T @ weighted_design
    cross = design.T @ (weights * z)
    try:
        projected = np.linalg.solve(nuisance_information, cross)
    except np.linalg.LinAlgError:
        projected = np.linalg.lstsq(nuisance_information, cross, rcond=None)[0]
    value = float(np.dot(weights * z, z) - np.dot(cross, projected))
    return max(value, 0.0)


def score_interaction(
    candidate_matrix: np.ndarray,
    y: np.ndarray,
    trait: str,
    covariates: np.ndarray | None = None,
) -> ScoreResult:
    """Efficient score test for the highest-order product, conditional on lower terms."""
    matrix = np.asarray(candidate_matrix, dtype=float)
    response = np.asarray(y, dtype=float)
    design = _lower_order_design(matrix, covariates)
    z = np.prod(matrix, axis=1)

    if trait == "binary":
        _, fitted, converged, iterations = fit_logistic_null(design, response)
        residual = response - fitted
        weights = np.clip(fitted * (1.0 - fitted), 1e-9, None)
        score = float(np.dot(z, residual))
        information = _efficient_information(z, design, weights)
    elif trait == "linear":
        beta = np.linalg.lstsq(design, response, rcond=None)[0]
        residual = response - design @ beta
        dof = max(1, len(response) - np.linalg.matrix_rank(design))
        sigma2 = max(float(np.dot(residual, residual) / dof), 1e-12)
        weights = np.full(len(response), 1.0 / sigma2)
        score = float(np.dot(z, residual) / sigma2)
        information = _efficient_information(z, design, weights)
        converged, iterations = True, 1
    else:
        raise ValueError("trait must be binary or linear")

    if not np.isfinite(information) or information <= 1e-12:
        return ScoreResult(0.0, 0.0, 0.0, 0.0, 1.0, converged, iterations)
    statistic = max(0.0, float(score * score / information))
    return ScoreResult(
        statistic=statistic,
        fi_gain=0.5 * statistic,
        beta_hat=float(score / information),
        information=information,
        p_value=float(chi2.sf(statistic, df=1)),
        null_converged=converged,
        null_iterations=iterations,
    )
