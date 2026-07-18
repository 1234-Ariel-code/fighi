import unittest

import numpy as np
import pandas as pd
from scipy.special import expit

from fighi.config import AnalysisConfig
from fighi.errors import CandidateLimitError
from fighi.search import FIGHI


class SearchTests(unittest.TestCase):
    def test_recovers_planted_pair(self):
        rng = np.random.default_rng(10)
        matrix = rng.binomial(2, 0.35, size=(600, 10))
        z = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0)
        y = rng.binomial(1, expit(1.3 * z[:, 1] * z[:, 5]))
        features = pd.DataFrame(matrix, columns=[f"rs{index}" for index in range(10)])
        result = FIGHI(
            AnalysisConfig(trait="binary", screen_method="all", top_m=10, max_order=2)
        ).fit(features, y)
        self.assertIn(("rs1", "rs5"), {item.features for item in result.significant})
        best = result.to_frame().iloc[0]
        self.assertEqual(best.hyperedge, "rs1|rs5")

    def test_candidate_limit_fails_before_search(self):
        rng = np.random.default_rng(3)
        features = pd.DataFrame(
            rng.binomial(2, 0.3, size=(100, 20)), columns=[f"rs{index}" for index in range(20)]
        )
        y = np.tile([0, 1], 50)
        model = FIGHI(
            AnalysisConfig(
                trait="binary",
                screen_method="all",
                max_candidates_per_order=100,
            )
        )
        with self.assertRaises(CandidateLimitError):
            model.fit(features, y)


if __name__ == "__main__":
    unittest.main()

