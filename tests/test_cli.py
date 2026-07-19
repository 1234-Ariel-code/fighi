import json
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree as ET

from fighi.cli import main


class CliTests(unittest.TestCase):
    def test_demo_end_to_end(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "demo"
            code = main(["demo", "--outdir", str(output), "--samples", "350"])
            self.assertEqual(code, 0)
            analysis = output / "analysis"
            expected = [
                "fighi_results.csv",
                "fighi_significant_interactions.csv",
                "fighi_feature_scores.csv",
                "fighi_summary.json",
                "fighi_model.json",
                "fighi_report.html",
                "fighi_cytoscape.cyjs",
                "fighi_hypergraph.graphml",
            ]
            for name in expected:
                self.assertTrue((analysis / name).is_file(), name)
            summary = json.loads((analysis / "fighi_summary.json").read_text())
            self.assertEqual(summary["software"]["version"], "1.1.0")
            self.assertGreater(summary["evaluated_interactions"], 0)
            ET.parse(analysis / "fighi_hypergraph.graphml")

    def test_nonempty_output_is_protected(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            (output / "keep.txt").write_text("user file", encoding="utf-8")
            code = main(["demo", "--outdir", str(output), "--samples", "100"])
            self.assertEqual(code, 2)
            self.assertEqual((output / "keep.txt").read_text(), "user file")


if __name__ == "__main__":
    unittest.main()
