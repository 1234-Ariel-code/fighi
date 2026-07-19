import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from fighi.plink import build_plink2_export_command, prepare_plink_dataset


class PlinkPreparationTests(unittest.TestCase):
    def _bed_fixture(self, root: Path) -> Path:
        prefix = root / "cohort"
        for suffix in (".bed", ".bim", ".fam"):
            prefix.with_suffix(suffix).write_text("fixture\n", encoding="utf-8")
        return prefix

    def test_build_command_is_explicit_and_safe(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prefix = self._bed_fixture(root)
            candidates = root / "candidates.txt"
            candidates.write_text("rs1\nrs2\n", encoding="utf-8")
            command, detected = build_plink2_export_command(
                prefix, candidates, root / "export", threads=3, memory_mb=512
            )
            self.assertEqual(detected, "bed")
            self.assertEqual(command[0], "plink2")
            self.assertIn("--extract", command)
            self.assertIn("include-alt", command)
            self.assertEqual(command[command.index("--threads") + 1], "3")

    def test_fake_plink_export_is_merged_and_hashed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prefix = self._bed_fixture(root)
            candidates = root / "candidates.txt"
            candidates.write_text("rs1\nrs2\n", encoding="utf-8")
            phenotype = root / "phenotype.tsv"
            phenotype.write_text("IID\tcase\tage\nI1\t0\t35\nI2\t1\t52\n", encoding="utf-8")
            fake = root / "plink2"
            fake.write_text(
                "#!/usr/bin/env bash\n"
                "while [ $# -gt 0 ]; do\n"
                '  if [ "$1" = "--out" ]; then out="$2"; shift 2; else shift; fi\n'
                "done\n"
                "printf 'FID IID PAT MAT SEX PHENOTYPE rs1_A rs2_G\\n0 I1 0 0 1 -9 0 2\\n0 I2 0 0 2 -9 1 1\\n' > \"${out}.raw\"\n",
                encoding="utf-8",
            )
            fake.chmod(0o755)
            paths = prepare_plink_dataset(
                prefix,
                candidates,
                root / "prepared",
                phenotype_file=phenotype,
                phenotype_column="case",
                covariate_columns=["age"],
                plink2=str(fake),
            )
            frame = pd.read_csv(paths["data"], sep="\t")
            self.assertEqual(frame.columns.tolist(), ["IID", "case", "age", "rs1", "rs2"])
            self.assertEqual(frame["rs2"].tolist(), [2, 1])
            manifest = json.loads(Path(paths["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "completed")
            self.assertEqual(manifest["exported_candidates"], 2)
            self.assertEqual(len(manifest["output_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
