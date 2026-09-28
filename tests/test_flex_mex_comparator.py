#!/usr/bin/env python3
import subprocess
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("compare_flex_hash_screen_mex.py")


class FlexMexComparatorTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="flex_mex_compare_")
        self.addCleanup(self.tmp.cleanup)
        self.a, self.b = (Path(self.tmp.name) / name for name in ("a", "b"))
        for root in (self.a, self.b):
            self.write_mex(root / "Solo.out/Gene/raw")

    @staticmethod
    def write_mex(path, value=1):
        path.mkdir(parents=True, exist_ok=True)
        (path / "features.tsv").write_text("ENSG1\tGENE1\n")
        (path / "barcodes.tsv").write_text("AAAA_BC006\n")
        (path / "matrix.mtx").write_text(f"%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 {value}\n")

    def compare(self, *args):
        return subprocess.run(["python3", str(SCRIPT), str(self.a), str(self.b), *args],
                              text=True, capture_output=True, timeout=10)

    def test_only_called_samples_required_in_auto_mode(self):
        for root in (self.a, self.b):
            self.write_mex(root / "per_sample/BC006/Gene/filtered")
        result = self.compare()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("emitted_samples=BC006", result.stdout)

    def test_no_called_samples_still_checks_raw(self):
        self.assertEqual(self.compare().returncode, 0)
        self.write_mex(self.b / "Solo.out/Gene/raw", value=2)
        self.assertNotEqual(self.compare().returncode, 0)

    def test_missing_sample_on_one_side_fails(self):
        self.write_mex(self.a / "per_sample/BC006/Gene/filtered")
        result = self.compare()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("per-sample sets differ", result.stderr)

    def test_explicit_missing_sample_fails(self):
        self.assertNotEqual(self.compare("--samples", "BC004").returncode, 0)

    def test_matching_directory_names_do_not_hide_missing_matrix(self):
        for root in (self.a, self.b):
            (root / "per_sample/BC004").mkdir(parents=True)
        self.assertNotEqual(self.compare().returncode, 0)


if __name__ == "__main__":
    unittest.main()
