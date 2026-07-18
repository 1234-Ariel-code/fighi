from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Literal

from .errors import InputValidationError


Trait = Literal["auto", "binary", "linear"]
Correction = Literal["fdr_bh", "bonferroni", "none"]
ScreenMethod = Literal["hybrid", "marginal", "variance", "all"]


@dataclass(slots=True)
class AnalysisConfig:
    """Configuration for one reproducible FIGHI analysis."""

    trait: Trait = "auto"
    max_order: int = 2
    top_m: int = 50
    screen_method: ScreenMethod = "hybrid"
    min_maf: float = 0.01
    max_missing: float = 0.10
    min_variance: float = 1e-8
    alpha: float = 0.05
    correction: Correction = "fdr_bh"
    discovery_fraction: float = 0.40
    beam_width: int = 100
    max_candidates_per_order: int = 100_000
    seed: int = 0
    stability_repeats: int = 0
    stability_fraction: float = 0.80
    stability_threshold: float = 0.70
    graph_top: int = 100

    def validate(self) -> "AnalysisConfig":
        if self.trait not in {"auto", "binary", "linear"}:
            raise InputValidationError("trait must be auto, binary, or linear")
        if not 2 <= self.max_order <= 6:
            raise InputValidationError("max_order must be between 2 and 6")
        if self.top_m < 2:
            raise InputValidationError("top_m must be at least 2")
        if self.screen_method not in {"hybrid", "marginal", "variance", "all"}:
            raise InputValidationError("unsupported screen_method")
        if not 0 <= self.min_maf < 0.5:
            raise InputValidationError("min_maf must be in [0, 0.5)")
        if not 0 <= self.max_missing < 1:
            raise InputValidationError("max_missing must be in [0, 1)")
        if not 0 < self.alpha < 1:
            raise InputValidationError("alpha must be in (0, 1)")
        if self.correction not in {"fdr_bh", "bonferroni", "none"}:
            raise InputValidationError("unsupported multiple-testing correction")
        if self.discovery_fraction != 0 and not 0.20 <= self.discovery_fraction <= 0.50:
            raise InputValidationError("discovery_fraction must be 0 or between 0.20 and 0.50")
        if self.beam_width < 1 or self.max_candidates_per_order < 1:
            raise InputValidationError("beam_width and candidate limit must be positive")
        if self.stability_repeats < 0:
            raise InputValidationError("stability_repeats cannot be negative")
        if not 0.5 <= self.stability_fraction < 1:
            raise InputValidationError("stability_fraction must be in [0.5, 1)")
        if not 0 <= self.stability_threshold <= 1:
            raise InputValidationError("stability_threshold must be in [0, 1]")
        return self

    def to_dict(self) -> dict:
        return asdict(self)
