import unittest

import numpy as np
import pandas as pd

from fighi.config import AnalysisConfig
from fighi.data import prepare_dataframe
from fighi.errors import InputValidationError


class DataTests(unittest.TestCase):
    def test_binary_mapping_and_qc(self):
        rng = np.random.default_rng(2)
        frame = pd.DataFrame(
            {
                "IID": [f"i{index}" for index in range(40)],
                "case": np.tile([1, 2], 20),
                "age": rng.normal(50, 5, 40),
                "rs1": rng.binomial(2, 0.3, 40),
                "rs2": rng.binomial(2, 0.4, 40).astype(float),
                "constant": np.ones(40),
            }
        )
        frame.loc[0, "rs2"] = np.nan
        data = prepare_dataframe(
            frame,
            "case",
            AnalysisConfig(screen_method="all", max_missing=0.2),
            id_column="IID",
            covariate_columns=["age"],
        )
        self.assertEqual(data.trait, "binary")
        self.assertEqual(set(data.phenotype), {0.0, 1.0})
        self.assertNotIn("constant", data.features)
        self.assertFalse(data.features.isna().any().any())

    def test_non_genotype_requires_opt_in(self):
        frame = pd.DataFrame(
            {
                "y": np.arange(30, dtype=float),
                "x1": np.linspace(-3, 3, 30),
                "x2": np.linspace(1, 8, 30),
            }
        )
        with self.assertRaises(InputValidationError):
            prepare_dataframe(frame, "y", AnalysisConfig(trait="linear"))


if __name__ == "__main__":
    unittest.main()
