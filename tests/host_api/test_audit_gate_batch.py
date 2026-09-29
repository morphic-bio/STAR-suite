import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from audit_gate_batch import TIER_A, audit


class GateAuditTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="star-host-gate-audit-")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.manifest = self.root / "manifest.tsv"
        self.manifest.write_text("# module\tcase\nscrna\tscrna\n")
        self.run = self.root / "run"
        (self.run / "manifest/scrna").mkdir(parents=True)
        (self.run / "tierA").mkdir()
        (self.run / "started_utc").write_text("2026-09-28T01:00:00Z\n")
        (self.run / "finished_utc").write_text("2026-09-28T02:00:00Z\n")
        (self.run / "manifest_status.tsv").write_text("scrna\t0\n")
        (self.run / "manifest/scrna/summary.tsv").write_text("PASS\tscrna\tscrna\n")
        (self.run / "tierA/status.tsv").write_text(
            "".join(f"{case}\t0\n" for case in TIER_A))

    def result(self):
        return audit(self.run, self.manifest)

    def test_complete_execution_does_not_assert_parity(self):
        result = self.result()
        self.assertTrue(result["execution_coverage_passed"])
        self.assertEqual(result["output_parity"], "not_evaluated")
        self.assertEqual(result["performance"], "not_evaluated")

    def test_failed_case_is_not_accepted(self):
        (self.run / "manifest_status.tsv").write_text("scrna\t1\n")
        (self.run / "manifest/scrna/summary.tsv").write_text("FAIL\tscrna\tscrna\n")
        result = self.result()
        self.assertFalse(result["execution_coverage_passed"])
        self.assertEqual(result["failures"]["production"], {"scrna": 1})

    def test_skipped_case_with_zero_exit_is_not_passed(self):
        (self.run / "manifest/scrna/summary.tsv").write_text("SKIP\tscrna\tscrna\n")
        self.assertEqual(self.result()["skipped"], ["scrna"])
        self.assertFalse(self.result()["execution_coverage_passed"])

    def test_successful_wrapper_with_downstream_failure(self):
        (self.run / "CELLBENDER_FAILED.txt").write_text("No denoised output\n")
        self.assertFalse(self.result()["execution_coverage_passed"])

    def test_nested_skip_requires_review(self):
        folder = self.run / "manifest/scrna/scrna_scrna"
        folder.mkdir()
        (folder / "stdout.log").write_text("PASS\treader\nSKIP\tnetwork\tRUN_NETWORK=0\n")
        result = self.result()
        self.assertFalse(result["execution_coverage_passed"])
        self.assertEqual(len(result["nested_skips_to_review"]), 1)

    def test_missing_completion_or_status(self):
        for name in ("finished_utc", "manifest_status.tsv", "tierA/status.tsv"):
            with self.subTest(name=name):
                path = self.run / name
                data = path.read_text()
                path.unlink()
                self.assertFalse(self.result()["execution_coverage_passed"])
                path.write_text(data)

    def test_invalid_completion_timestamp(self):
        for value in ("", "yesterday", "2026-09-28T00:00:00Z", "2026-09-28T02:00:00"):
            with self.subTest(value=value):
                (self.run / "finished_utc").write_text(value)
                self.assertFalse(self.result()["execution_coverage_passed"])

    def test_duplicate_unknown_or_invalid_status(self):
        for row in ("scrna\t0\n", "unknown\t0\n", "bad\n", "bad\tPASS\n"):
            with self.subTest(row=row):
                (self.run / "manifest_status.tsv").write_text("scrna\t0\n" + row)
                self.assertFalse(self.result()["execution_coverage_passed"])

    def test_tier_a_failure(self):
        path = self.run / "tierA/status.tsv"
        path.write_text(path.read_text().replace("\t0", "\t1", 1))
        self.assertFalse(self.result()["execution_coverage_passed"])

    def test_missing_or_duplicate_summary(self):
        path = self.run / "manifest/scrna/summary.tsv"
        path.unlink()
        self.assertFalse(self.result()["execution_coverage_passed"])
        path.write_text("PASS\tscrna\tscrna\nPASS\tscrna\tscrna\n")
        self.assertFalse(self.result()["execution_coverage_passed"])

    def test_summary_exit_disagreement(self):
        (self.run / "manifest/scrna/summary.tsv").write_text("FAIL\tscrna\tscrna\n")
        self.assertFalse(self.result()["execution_coverage_passed"])

    def test_cli_nonzero_and_report_preservation(self):
        report = self.root / "audit.json"
        (self.run / "CELLBENDER_FAILED.txt").touch()
        cmd = [sys.executable, str(Path(__file__).with_name("audit_gate_batch.py")),
               str(self.run), "--manifest", str(self.manifest), "--report", str(report)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 1)
        payload = report.read_text()
        self.assertFalse(json.loads(payload)["execution_coverage_passed"])
        proc = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 2)
        self.assertEqual(report.read_text(), payload)


if __name__ == "__main__":
    unittest.main()
