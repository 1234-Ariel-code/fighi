import unittest

import numpy as np
from scipy.special import expit

from fighi.statistics import bh_adjust, score_interaction


class StatisticsTests(unittest.TestCase):
    def test_bh_adjustment(self):
        adjusted = bh_adjust(np.array([0.01, 0.02, 0.20]))
        np.testing.assert_allclose(adjusted, [0.03, 0.03, 0.20])

    def test_logistic_interaction_score(self):
        rng = np.random.default_rng(4)
        x = rng.normal(size=(800, 2))
        covariate = rng.normal(size=(800, 1))
        probability = expit(1.15 * x[:, 0] * x[:, 1] + 0.25 * covariate[:, 0])
        y = rng.binomial(1, probability)
        result = score_interaction(x, y, "binary", covariate)
        self.assertTrue(result.null_converged)
        self.assertLess(result.p_value, 1e-12)
        self.assertGreater(result.beta_hat, 0)

    def test_linear_interaction_score(self):
        rng = np.random.default_rng(8)
        x = rng.normal(size=(500, 2))
        y = 1.8 * x[:, 0] * x[:, 1] + rng.normal(scale=0.7, size=500)
        result = score_interaction(x, y, "linear")
        self.assertLess(result.p_value, 1e-30)
        self.assertGreater(result.beta_hat, 1.0)


if __name__ == "__main__":
    unittest.main()
