from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import norm

from .errors import InputValidationError
from .utilities import sha256_file, write_json

SCENARIOS = ("null", "main", "pairwise", "pairwise_main", "threeway", "structure")


@dataclass(slots=True)
class SimulationConfig:
    samples: int = 1000
    features: int = 100
    trait: str = "binary"
    scenario: str = "pairwise"
    effect: float = 1.0
    prevalence: float = 0.5
    noise_sd: float = 1.0
    min_maf: float = 0.10
    max_maf: float = 0.45
    ld_rho: float = 0.0
    seed: int = 17

    def validate(self) -> SimulationConfig:
        if self.samples < 20:
            raise InputValidationError("Simulation requires at least 20 samples")
        if self.features < 3:
            raise InputValidationError("Simulation requires at least 3 features")
        if self.trait not in {"binary", "linear"}:
            raise InputValidationError("Simulation trait must be binary or linear")
        if self.scenario not in SCENARIOS:
            raise InputValidationError(f"Unknown scenario: {self.scenario}")
        if not 0.0 < self.min_maf <= self.max_maf < 0.5:
            raise InputValidationError("MAF bounds must satisfy 0 < min-maf <= max-maf < 0.5")
        if not 0.0 < self.prevalence < 1.0:
            raise InputValidationError("Prevalence must be strictly between 0 and 1")
        if self.noise_sd <= 0:
            raise InputValidationError("noise-sd must be positive")
        if not 0.0 <= self.ld_rho < 1.0:
            raise InputValidationError("ld-rho must satisfy 0 <= ld-rho < 1")
        return self


def _genotypes(
    rng: np.random.Generator,
    samples: int,
    mafs: np.ndarray,
    ld_rho: float,
) -> np.ndarray:
    if ld_rho <= 0:
        return rng.binomial(2, mafs, size=(samples, len(mafs))).astype(float)

    latent = np.empty((samples, len(mafs)), dtype=float)
    latent[:, 0] = rng.normal(size=samples)
    innovation = np.sqrt(1.0 - ld_rho**2)
    for index in range(1, len(mafs)):
        latent[:, index] = ld_rho * latent[:, index - 1] + innovation * rng.normal(size=samples)
    uniforms = norm.cdf(latent)
    p_zero = (1.0 - mafs) ** 2
    p_at_most_one = 1.0 - mafs**2
    return (uniforms > p_zero).astype(float) + (uniforms > p_at_most_one).astype(float)


def _standardize(values: np.ndarray) -> np.ndarray:
    scales = values.std(axis=0)
    scales[scales <= 0] = 1.0
    return (values - values.mean(axis=0)) / scales


def _binary_intercept(linear: np.ndarray, prevalence: float) -> float:
    lower, upper = -30.0, 30.0
    for _ in range(100):
        midpoint = (lower + upper) / 2.0
        if float(expit(midpoint + linear).mean()) < prevalence:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def simulate_dataset(config: SimulationConfig) -> tuple[pd.DataFrame, dict]:
    """Generate a reproducible genotype interaction benchmark with known truth."""
    config.validate()
    rng = np.random.default_rng(config.seed)
    mafs = rng.uniform(config.min_maf, config.max_maf, size=config.features)

    ancestry = rng.binomial(1, 0.5, size=config.samples).astype(float)
    if config.scenario == "structure":
        shifts = rng.uniform(-0.08, 0.08, size=config.features)
        base = np.tile(mafs, (config.samples, 1))
        sample_mafs = np.clip(base + (ancestry[:, None] - 0.5) * shifts, 0.02, 0.48)
        genotype = rng.binomial(2, sample_mafs).astype(float)
    else:
        genotype = _genotypes(rng, config.samples, mafs, config.ld_rho)

    standardized = _standardize(genotype)
    age = np.clip(rng.normal(50.0, 12.0, size=config.samples), 18.0, 90.0)
    age_scaled = (age - age.mean()) / age.std()
    ancestry_scaled = (ancestry - ancestry.mean()) / max(ancestry.std(), 1e-12)
    linear = 0.20 * age_scaled
    truth: list[list[str]] = []

    feature_names = [f"rs_sim_{index + 1:05d}" for index in range(config.features)]
    if config.scenario == "main":
        linear = linear + config.effect * standardized[:, 0]
    elif config.scenario == "pairwise":
        linear = linear + config.effect * standardized[:, 0] * standardized[:, 1]
        truth.append(feature_names[:2])
    elif config.scenario == "pairwise_main":
        linear = linear + 0.35 * standardized[:, 0] + 0.35 * standardized[:, 1]
        linear = linear + config.effect * standardized[:, 0] * standardized[:, 1]
        truth.append(feature_names[:2])
    elif config.scenario == "threeway":
        linear = linear + config.effect * np.prod(standardized[:, :3], axis=1)
        truth.append(feature_names[:3])
    elif config.scenario == "structure":
        linear = linear + 0.9 * ancestry_scaled

    if config.trait == "binary":
        intercept = _binary_intercept(linear, config.prevalence)
        probabilities = expit(intercept + linear)
        phenotype = rng.binomial(1, probabilities)
        phenotype_name = "case"
    else:
        phenotype = linear + rng.normal(0.0, config.noise_sd, size=config.samples)
        phenotype_name = "trait"

    frame = pd.DataFrame(
        {
            "IID": [f"SIM{index + 1:07d}" for index in range(config.samples)],
            phenotype_name: phenotype,
            "age": age.round(4),
            "PC1": ancestry_scaled.round(6),
        }
    )
    frame = pd.concat([frame, pd.DataFrame(genotype.astype(int), columns=feature_names)], axis=1)
    metadata = {
        "schema_version": "1.0",
        "config": asdict(config),
        "phenotype_column": phenotype_name,
        "covariate_columns": ["age", "PC1"],
        "feature_count": config.features,
        "sample_count": config.samples,
        "truth_interactions": truth,
        "truth_hyperedges": ["|".join(sorted(item)) for item in truth],
        "causal_main_effects": [feature_names[0]] if config.scenario == "main" else [],
        "maf": {name: float(maf) for name, maf in zip(feature_names, mafs, strict=True)},
    }
    return frame, metadata


def write_simulation(outdir: str | Path, config: SimulationConfig) -> dict[str, str]:
    """Write a compressed data table and complete simulation provenance."""
    target = Path(outdir)
    target.mkdir(parents=True, exist_ok=True)
    frame, metadata = simulate_dataset(config)
    data_path = target / "simulation.tsv.gz"
    frame.to_csv(data_path, sep="\t", index=False, compression="gzip")

    features = [name for name in frame.columns if name.startswith("rs_sim_")]
    candidate_path = target / "candidates.txt"
    candidate_path.write_text("\n".join(features) + "\n", encoding="utf-8")
    sample_path = target / "samples.txt"
    sample_path.write_text(
        "\n".join(f"0\t{sample}" for sample in frame["IID"].astype(str)) + "\n",
        encoding="utf-8",
    )
    metadata["files"] = {
        "data": {"path": data_path.name, "sha256": sha256_file(data_path)},
        "candidates": {"path": candidate_path.name, "sha256": sha256_file(candidate_path)},
        "samples": {"path": sample_path.name, "sha256": sha256_file(sample_path)},
    }
    truth_path = write_json(target / "truth.json", metadata)
    return {
        "data": str(data_path),
        "candidates": str(candidate_path),
        "samples": str(sample_path),
        "truth": str(truth_path),
    }
