import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from fighi.data import load_table
from fighi.simulation import SimulationConfig, simulate_dataset, write_simulation


class SimulationTests(unittest.TestCase):
    def test_pairwise_simulation_is_reproducible(self):
        config = SimulationConfig(samples=120, features=8, scenario="pairwise", seed=9)
        left, left_truth = simulate_dataset(config)
        right, right_truth = simulate_dataset(config)
        pd.testing.assert_frame_equal(left, right)
        self.assertEqual(left_truth, right_truth)
        self.assertEqual(left_truth["truth_hyperedges"], ["rs_sim_00001|rs_sim_00002"])

    def test_compressed_output_is_fighi_readable(self):
        with tempfile.TemporaryDirectory() as temporary:
            paths = write_simulation(
                temporary,
                SimulationConfig(samples=100, features=7, trait="linear", scenario="threeway"),
            )
            frame = load_table(paths["data"])
            self.assertEqual(frame.shape, (100, 11))
            truth = json.loads(Path(paths["truth"]).read_text(encoding="utf-8"))
            self.assertEqual(truth["truth_hyperedges"], ["rs_sim_00001|rs_sim_00002|rs_sim_00003"])
            self.assertIn("sha256", truth["files"]["data"])


if __name__ == "__main__":
    unittest.main()
