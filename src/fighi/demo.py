from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit


def create_demo_dataset(path: str | Path, samples: int = 500, seed: int = 17) -> Path:
    rng = np.random.default_rng(seed)
    feature_names = [f"rs_demo_{index:02d}" for index in range(1, 21)]
    frequencies = np.linspace(0.18, 0.42, len(feature_names))
    matrix = np.column_stack(
        [rng.binomial(2, frequency, size=samples) for frequency in frequencies]
    )
    standardized = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0)
    age = rng.normal(50, 12, size=samples)
    interaction = standardized[:, 2] * standardized[:, 6]
    probability = expit(-0.25 + 1.35 * interaction + 0.018 * (age - 50))
    phenotype = rng.binomial(1, probability)
    frame = pd.DataFrame(matrix, columns=feature_names)
    frame.insert(0, "age", np.round(age, 2))
    frame.insert(0, "case", phenotype)
    frame.insert(0, "IID", [f"demo_{index:04d}" for index in range(1, samples + 1)])
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(target, index=False)
    return target
