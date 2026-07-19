from __future__ import annotations

import math
import time
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from itertools import combinations

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .data import PreparedData, prepare_dataframe
from .errors import CandidateLimitError, InputValidationError
from .statistics import adjust_pvalues, fit_logistic_null, score_interaction


@dataclass(slots=True)
class InteractionResult:
    features: tuple[str, ...]
    order: int
    statistic: float
    fi_gain: float
    beta_hat: float
    information: float
    p_value: float
    q_order: float = math.nan
    q_value: float = math.nan
    significant: bool = False
    selected_for_expansion: bool = False
    stability: float | None = None
    null_converged: bool = True
    null_iterations: int = 0

    @property
    def odds_ratio(self) -> float:
        return float(np.exp(np.clip(self.beta_hat, -50, 50)))

    def to_dict(self, trait: str) -> dict:
        values = asdict(self)
        values["hyperedge"] = "|".join(self.features)
        values["odds_ratio"] = self.odds_ratio if trait == "binary" else None
        return values


@dataclass(slots=True)
class SearchResult:
    trait: str
    config: AnalysisConfig
    atoms: list[str]
    interactions: list[InteractionResult]
    qc: dict
    warnings: list[str]
    runtime_seconds: float
    evaluated_by_order: dict[int, int]
    phenotype_mapping: dict[str, float] | None = None
    stopped_reason: str | None = None

    def to_frame(self) -> pd.DataFrame:
        rows = [item.to_dict(self.trait) for item in self.interactions]
        columns = [
            "hyperedge",
            "features",
            "order",
            "statistic",
            "fi_gain",
            "beta_hat",
            "odds_ratio",
            "information",
            "p_value",
            "q_order",
            "q_value",
            "significant",
            "selected_for_expansion",
            "stability",
            "null_converged",
            "null_iterations",
        ]
        if not rows:
            return pd.DataFrame(columns=columns)
        return (
            pd.DataFrame(rows)
            .loc[:, columns]
            .sort_values(["q_value", "p_value", "fi_gain"], ascending=[True, True, False])
        )

    @property
    def significant(self) -> list[InteractionResult]:
        return [item for item in self.interactions if item.significant]


def _apriori_candidates(previous: Iterable[tuple[str, ...]], order: int) -> list[tuple[str, ...]]:
    previous_set = {tuple(sorted(item)) for item in previous}
    candidates: set[tuple[str, ...]] = set()
    previous_list = sorted(previous_set)
    for left_index, left in enumerate(previous_list):
        for right in previous_list[left_index + 1 :]:
            union = tuple(sorted(set(left) | set(right)))
            if len(union) != order:
                continue
            if all(tuple(part) in previous_set for part in combinations(union, order - 1)):
                candidates.add(union)
    return sorted(candidates)


class FIGHI:
    """High-level interaction search API."""

    def __init__(self, config: AnalysisConfig | None = None):
        self.config = (config or AnalysisConfig()).validate()

    def _residualized_response(self, data: PreparedData) -> np.ndarray:
        y = data.phenotype.astype(float)
        if data.covariates is None or not data.covariates.size:
            return y - y.mean()
        design = np.column_stack([np.ones(len(y)), data.covariates])
        if data.trait == "binary":
            _, fitted, _, _ = fit_logistic_null(design, y)
            return y - fitted
        beta = np.linalg.lstsq(design, y, rcond=None)[0]
        return y - design @ beta

    @staticmethod
    def _subset(data: PreparedData, indices: np.ndarray) -> PreparedData:
        covariates = data.covariates[indices] if data.covariates is not None else None
        sample_ids = [data.sample_ids[index] for index in indices] if data.sample_ids else None
        return PreparedData(
            features=data.features.iloc[indices].reset_index(drop=True),
            standardized=data.standardized.iloc[indices].reset_index(drop=True),
            phenotype=data.phenotype[indices],
            trait=data.trait,
            covariates=covariates,
            covariate_names=data.covariate_names,
            sample_ids=sample_ids,
            qc=data.qc,
            warnings=data.warnings,
            phenotype_mapping=data.phenotype_mapping,
        )

    def _split(self, data: PreparedData) -> tuple[PreparedData, PreparedData, str]:
        fraction = self.config.discovery_fraction
        if fraction == 0:
            return data, data, "exploratory_all_samples"
        rng = np.random.default_rng(self.config.seed)
        if data.trait == "binary":
            discovery_parts = []
            inference_parts = []
            for label in [0, 1]:
                indices = np.flatnonzero(data.phenotype == label)
                rng.shuffle(indices)
                split = int(round(len(indices) * fraction))
                split = min(max(split, 5), len(indices) - 5)
                if split < 5 or len(indices) - split < 5:
                    raise InputValidationError(
                        "Holdout inference needs at least 10 samples per phenotype class. "
                        "Use more data or --discovery-fraction 0 for explicitly exploratory analysis."
                    )
                discovery_parts.append(indices[:split])
                inference_parts.append(indices[split:])
            discovery_indices = np.concatenate(discovery_parts)
            inference_indices = np.concatenate(inference_parts)
        else:
            indices = rng.permutation(len(data.phenotype))
            split = int(round(len(indices) * fraction))
            if split < 20 or len(indices) - split < 20:
                raise InputValidationError(
                    "Holdout inference needs at least 20 discovery and 20 inference samples. "
                    "Use more data or --discovery-fraction 0 for explicitly exploratory analysis."
                )
            discovery_indices, inference_indices = indices[:split], indices[split:]
        return (
            self._subset(data, np.sort(discovery_indices)),
            self._subset(data, np.sort(inference_indices)),
            "independent_holdout",
        )

    def _screen(self, data: PreparedData) -> list[str]:
        columns = data.standardized.columns.to_numpy()
        if self.config.screen_method == "all":
            return columns.tolist()
        residual = self._residualized_response(data)
        matrix = data.standardized.to_numpy(dtype=float)
        denominator = max(float(np.linalg.norm(residual)), 1e-12)
        marginal = np.abs(matrix.T @ residual) / denominator
        variance = data.features.var(axis=0, ddof=0).to_numpy(dtype=float)
        count = min(self.config.top_m, len(columns))
        marginal_order = np.argsort(-np.nan_to_num(marginal))
        variance_order = np.argsort(-np.nan_to_num(variance))
        if self.config.screen_method == "marginal":
            selected = marginal_order[:count]
        elif self.config.screen_method == "variance":
            selected = variance_order[:count]
        else:
            marginal_count = max(1, int(math.ceil(count * 0.75)))
            chosen = list(marginal_order[:marginal_count])
            for index in variance_order:
                if int(index) not in chosen:
                    chosen.append(int(index))
                if len(chosen) >= count:
                    break
            selected = np.asarray(chosen)
        return columns[selected].tolist()

    def _check_candidate_count(self, count: int, order: int) -> None:
        if count > self.config.max_candidates_per_order:
            raise CandidateLimitError(
                f"Order {order} would evaluate {count:,} candidates, above the safety limit "
                f"of {self.config.max_candidates_per_order:,}. Reduce --top-m/--beam-width, "
                "lower --max-order, or deliberately raise --max-candidates-per-order."
            )

    def _evaluate_order(
        self, data: PreparedData, candidates: list[tuple[str, ...]]
    ) -> list[InteractionResult]:
        self._check_candidate_count(len(candidates), len(candidates[0]) if candidates else 2)
        results: list[InteractionResult] = []
        standardized = data.standardized
        for candidate in candidates:
            score = score_interaction(
                standardized.loc[:, list(candidate)].to_numpy(dtype=float),
                data.phenotype,
                data.trait,
                data.covariates,
            )
            results.append(
                InteractionResult(
                    features=candidate,
                    order=len(candidate),
                    statistic=score.statistic,
                    fi_gain=score.fi_gain,
                    beta_hat=score.beta_hat,
                    information=score.information,
                    p_value=score.p_value,
                    null_converged=score.null_converged,
                    null_iterations=score.null_iterations,
                )
            )
        if results:
            adjusted = adjust_pvalues(
                np.asarray([item.p_value for item in results]), self.config.correction
            )
            for item, q_value in zip(results, adjusted, strict=True):
                item.q_order = float(q_value)
        return results

    def _choose_expansion(self, results: list[InteractionResult]) -> list[tuple[str, ...]]:
        ranked = sorted(results, key=lambda item: (item.q_order, item.p_value, -item.fi_gain))
        selected = [item for item in ranked if item.q_order <= self.config.alpha]
        if len(selected) < self.config.beam_width:
            selected_ids = {item.features for item in selected}
            selected.extend(item for item in ranked if item.features not in selected_ids)
        selected = selected[: self.config.beam_width]
        for item in selected:
            item.selected_for_expansion = True
        return [item.features for item in selected]

    def _stability(self, data: PreparedData, interactions: list[InteractionResult]) -> None:
        repeats = self.config.stability_repeats
        if repeats <= 0 or not interactions:
            return
        targets = sorted(interactions, key=lambda item: (item.q_value, item.p_value))[
            : self.config.graph_top
        ]
        counts = np.zeros(len(targets), dtype=int)
        rng = np.random.default_rng(self.config.seed)
        n = len(data.phenotype)
        size = max(20, int(round(n * self.config.stability_fraction)))
        for _ in range(repeats):
            if data.trait == "binary":
                zero = np.flatnonzero(data.phenotype == 0)
                one = np.flatnonzero(data.phenotype == 1)
                n_zero = max(5, int(round(len(zero) * self.config.stability_fraction)))
                n_one = max(5, int(round(len(one) * self.config.stability_fraction)))
                sample = np.concatenate(
                    [rng.choice(zero, n_zero, replace=False), rng.choice(one, n_one, replace=False)]
                )
            else:
                sample = rng.choice(n, min(size, n), replace=False)
            sample_p = []
            for target in targets:
                matrix = data.standardized.loc[:, list(target.features)].to_numpy(dtype=float)[
                    sample
                ]
                covariates = data.covariates[sample] if data.covariates is not None else None
                score = score_interaction(matrix, data.phenotype[sample], data.trait, covariates)
                sample_p.append(score.p_value)
            adjusted = adjust_pvalues(np.asarray(sample_p), self.config.correction)
            counts += adjusted <= self.config.alpha
        for target, count in zip(targets, counts, strict=True):
            target.stability = float(count / repeats)

    def run(self, data: PreparedData) -> SearchResult:
        started = time.perf_counter()
        discovery, inference, inference_mode = self._split(data)
        atoms = self._screen(discovery)
        candidate_orders: list[list[tuple[str, ...]]] = []
        expansion_selected: set[tuple[str, ...]] = set()
        evaluated: dict[int, int] = {}
        stopped_reason: str | None = None
        previous: list[tuple[str, ...]] = []
        for order in range(2, self.config.max_order + 1):
            if order == 2:
                candidate_count = math.comb(len(atoms), 2)
                self._check_candidate_count(candidate_count, order)
                candidates = list(combinations(atoms, 2))
            else:
                candidates = _apriori_candidates(previous, order)
                self._check_candidate_count(len(candidates), order)
            if not candidates:
                stopped_reason = f"No candidates passed hierarchical expansion at order {order}"
                break
            discovery_results = self._evaluate_order(discovery, candidates)
            candidate_orders.append(candidates)
            evaluated[order] = len(discovery_results)
            previous = self._choose_expansion(discovery_results)
            if order < self.config.max_order:
                expansion_selected.update(previous)

        interactions: list[InteractionResult] = []
        for candidates in candidate_orders:
            inference_results = self._evaluate_order(inference, candidates)
            for item in inference_results:
                item.selected_for_expansion = item.features in expansion_selected
            interactions.extend(inference_results)

        if interactions:
            global_adjusted = adjust_pvalues(
                np.asarray([item.p_value for item in interactions]), self.config.correction
            )
            for item, q_value in zip(interactions, global_adjusted, strict=True):
                item.q_value = float(q_value)
                item.significant = bool(q_value <= self.config.alpha)
        self._stability(inference, interactions)
        warnings = list(data.warnings)
        non_converged = sum(not item.null_converged for item in interactions)
        if non_converged:
            warnings.append(
                f"{non_converged} candidate null models did not meet the convergence tolerance"
            )
        if self.config.screen_method != "all":
            warnings.append(
                "Feature screening is a computational filter and may miss pure interactions with weak marginal signal"
            )
        if inference_mode == "exploratory_all_samples":
            warnings.append(
                "Exploratory mode reused samples for discovery and inference; adjusted p-values do not account for data-dependent screening"
            )
        qc = dict(data.qc)
        qc["inference_mode"] = inference_mode
        qc["discovery_samples"] = len(discovery.phenotype)
        qc["inference_samples"] = len(inference.phenotype)
        return SearchResult(
            trait=data.trait,
            config=self.config,
            atoms=atoms,
            interactions=interactions,
            qc=qc,
            warnings=warnings,
            runtime_seconds=time.perf_counter() - started,
            evaluated_by_order=evaluated,
            phenotype_mapping=data.phenotype_mapping,
            stopped_reason=stopped_reason,
        )

    def fit(
        self,
        features: pd.DataFrame,
        phenotype: np.ndarray | pd.Series,
        covariates: pd.DataFrame | None = None,
        allow_non_genotype: bool = False,
    ) -> SearchResult:
        """Run FIGHI from in-memory objects using the same validation as the CLI."""
        features = features.copy()
        features.columns = features.columns.astype(str)
        frame = features.copy()
        frame["__fighi_phenotype__"] = np.asarray(phenotype)
        covariate_names: list[str] = []
        if covariates is not None:
            for name in covariates.columns:
                safe_name = f"__covariate_{name}"
                frame[safe_name] = covariates[name].to_numpy()
                covariate_names.append(safe_name)
        prepared = prepare_dataframe(
            frame,
            "__fighi_phenotype__",
            self.config,
            covariate_columns=covariate_names,
            feature_names=features.columns.astype(str).tolist(),
            allow_non_genotype=allow_non_genotype,
        )
        return self.run(prepared)
