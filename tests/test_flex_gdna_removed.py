#!/usr/bin/env python3
"""Prevent reintroduction of the withdrawn Flex estimator."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
STAR = Path(os.environ.get("STAR_BIN", ROOT / "core/legacy/source/STAR")).resolve()


class RemovedDiagnosticTests(unittest.TestCase):
    def test_source_removed(self):
        for path in ("flex/source/FlexGdna.cpp", "flex/source/FlexGdna.h",
                     "core/legacy/test/test_flex_gdna.cpp"):
            self.assertFalse((ROOT / path).exists(), path)

    def test_binary_has_no_estimator(self):
        symbols = subprocess.check_output(["nm", "-C", str(STAR)], text=True)
        self.assertNotRegex(symbols, r"FlexGdnaProbeMetadata|flexGdnaEstimate|writeGdnaJson")

    def test_retired_options(self):
        with tempfile.TemporaryDirectory(prefix="star-gdna-removed-") as temp:
            base = [
                str(STAR), "--genomeDir", str(Path(temp) / "missing-index"),
                "--readFilesIn", str(ROOT / "tests/fixtures/trim_qc_fastq_tiny.fastq"),
                str(ROOT / "tests/fixtures/trim_qc_fastq_tiny.fastq"),
                "--soloType", "CB_UMI_Simple", "--soloCBwhitelist", "None",
                "--soloCBlen", "16", "--soloUMIstart", "17", "--soloUMIlen", "12",
                "--soloFeatures", "Gene", "--outSAMtype", "None",
            ]
            for flex in (False, True):
                extra = (["--flex", "yes", "--flexLegacy", "yes", "--soloProbeList",
                          str(ROOT / "tests/fixtures/flex_probe_gene_list_tiny.txt")]
                         if flex else [])
                for option, value, rejected in (
                    ("--soloFlexGdna", "auto", False),
                    ("--soloFlexGdna", "no", False),
                    ("--soloFlexGdna", "yes", True),
                    ("--soloFlexGdna", "invalid", True),
                    ("--soloFlexGdnaProbeSet", "/missing/probes.csv", True),
                ):
                    with self.subTest(flex=flex, option=option, value=value):
                        prefix = str(Path(temp) / f"{flex}-{option}-{value.split('/')[-1]}.")
                        result = subprocess.run(
                            base + extra + [option, value, "--outFileNamePrefix", prefix],
                            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                        self.assertNotEqual(result.returncode, 0, result.stdout)
                        if rejected:
                            self.assertIn("diagnostic has been removed", result.stdout)
                        else:
                            self.assertNotIn("diagnostic has been removed", result.stdout)
                            self.assertRegex(result.stdout, r"could not open genome file")
            self.assertFalse(list(Path(temp).rglob("*gdna*.json")))
            self.assertFalse(list(Path(temp).rglob("*gdna*.tsv")))


if __name__ == "__main__":
    unittest.main()
