import tempfile
import unittest
from pathlib import Path

import pandas as pd

from fighi.annotation import annotate_results
from fighi.conversion import convert_tped


class ConversionAnnotationTests(unittest.TestCase):
    def test_tped_conversion(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "toy.tfam").write_text("F1 I1 0 0 1 -9\nF2 I2 0 0 2 -9\n", encoding="utf-8")
            (root / "toy.tped").write_text(
                "1 rs1 0 100 A A A G\n1 rs2 0 200 C T T T\n", encoding="utf-8"
            )
            (root / "toy.phen").write_text("F1 I1 1\nF2 I2 2\n", encoding="utf-8")
            output = convert_tped(
                root / "toy.tped", root / "toy.tfam", root / "toy.csv", root / "toy.phen", "case"
            )
            frame = pd.read_csv(output)
            self.assertEqual(frame.columns.tolist(), ["IID", "case", "rs1", "rs2"])
            self.assertEqual(frame.shape, (2, 4))

    def test_local_annotation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pd.DataFrame(
                {
                    "feature": ["rs1", "rs2", "rs3"],
                    "significant_interaction_count": [1, 1, 0],
                    "fi_contribution": [4.0, 3.0, 0.2],
                }
            ).to_csv(root / "scores.csv", index=False)
            pd.DataFrame({"SNP": ["rs1", "rs2", "rs3"], "Gene": ["A", "B", "C"]}).to_csv(
                root / "map.csv", index=False
            )
            (root / "pathways.gmt").write_text(
                "AB_PATH\tdescription\tA\tB\nOTHER\tdescription\tC\n", encoding="utf-8"
            )
            paths = annotate_results(
                root / "scores.csv", root / "map.csv", root / "out", [root / "pathways.gmt"]
            )
            annotated = pd.read_csv(paths["annotated_features"])
            self.assertEqual(annotated.loc[0, "gene"], "A")
            self.assertTrue(Path(paths["enrichment"]).is_file())


if __name__ == "__main__":
    unittest.main()
