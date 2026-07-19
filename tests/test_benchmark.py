import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from fighi.benchmark import run_benchmark, validate_benchmark_manifest


class BenchmarkTests(unittest.TestCase):
    def test_generic_methods_are_harmonized_and_compared(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "candidates.txt").write_text("rs1\nrs2\nrs3\n", encoding="utf-8")
            (root / "samples.txt").write_text("0 I1\n0 I2\n", encoding="utf-8")
            (root / "truth.json").write_text(
                json.dumps({"truth_hyperedges": ["rs1|rs2"]}), encoding="utf-8"
            )
            writer = (
                "from pathlib import Path; import sys; "
                "Path(sys.argv[1]).write_text(sys.argv[2], encoding='utf-8')"
            )
            methods = []
            for name, rows in [
                ("alpha", "a,b,p\nrs1,rs2,0.001\nrs1,rs3,0.4\n"),
                ("beta", "a,b,p\nrs1,rs2,0.01\nrs2,rs3,0.02\n"),
            ]:
                methods.append(
                    {
                        "name": name,
                        "command": [
                            sys.executable,
                            "-c",
                            writer,
                            "{run_dir}/results.csv",
                            rows,
                        ],
                        "result": {
                            "path": "{run_dir}/results.csv",
                            "format": "generic",
                            "feature_columns": ["a", "b"],
                            "p_column": "p",
                        },
                    }
                )
            methods.append(
                {
                    "name": "gamma_rank_only",
                    "command": [
                        sys.executable,
                        "-c",
                        writer,
                        "{run_dir}/results.csv",
                        "a,b,score\nrs1,rs2,0.9\nrs1,rs3,0.2\n",
                    ],
                    "result": {
                        "path": "{run_dir}/results.csv",
                        "format": "generic",
                        "feature_columns": ["a", "b"],
                        "p_column": None,
                        "score_column": "score",
                        "score_ascending": False,
                    },
                }
            )
            manifest = root / "benchmark.json"
            manifest.write_text(
                json.dumps(
                    {
                        "analysis_id": "unit",
                        "alpha": 0.05,
                        "correction": "fdr_bh",
                        "top_n": 1,
                        "shared_inputs": {
                            "candidate_file": "candidates.txt",
                            "sample_file": "samples.txt",
                            "truth_file": "truth.json",
                        },
                        "methods": methods,
                    }
                ),
                encoding="utf-8",
            )
            validation = validate_benchmark_manifest(manifest)
            self.assertTrue(validation["valid"])
            paths = run_benchmark(manifest, root / "out", plots=False)
            summary = pd.read_csv(paths["method_summary"])
            self.assertEqual(set(summary["method"]), {"alpha", "beta", "gamma_rank_only"})
            self.assertTrue((summary["top_n_recall"] == 1.0).all())
            gamma = pd.read_csv(paths["interactions"]).query("method == 'gamma_rank_only'")
            self.assertTrue(gamma["p_value"].isna().all())
            self.assertTrue(gamma["q_harmonized"].isna().all())
            overlap = pd.read_csv(paths["overlap"])
            self.assertEqual(overlap.loc[0, "jaccard"], 1.0)
            self.assertTrue(Path(paths["report"]).is_file())

    def test_timeout_terminates_method_and_is_reported(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest = root / "timeout.json"
            manifest.write_text(
                json.dumps(
                    {
                        "analysis_id": "timeout",
                        "methods": [
                            {
                                "name": "slow",
                                "command": [
                                    sys.executable,
                                    "-c",
                                    "import time; time.sleep(30)",
                                ],
                                "timeout_seconds": 0.1,
                                "result": {
                                    "path": "{run_dir}/never.csv",
                                    "format": "generic",
                                    "feature_columns": ["a", "b"],
                                    "p_column": "p",
                                },
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            paths = run_benchmark(manifest, root / "out", plots=False)
            runs = pd.read_csv(paths["runs"])
            self.assertEqual(runs.loc[0, "status"], "timeout")
            self.assertLess(runs.loc[0, "wall_seconds"], 5.0)
            self.assertTrue(pd.read_csv(paths["method_summary"]).empty)
            self.assertTrue(pd.read_csv(paths["overlap"]).empty)


if __name__ == "__main__":
    unittest.main()
