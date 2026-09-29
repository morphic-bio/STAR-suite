#!/usr/bin/env python3
"""Exercise the smoke's telemetry checks without running STAR or reading FASTQs."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


SCRIPT = Path(__file__).with_name("run_pf_dynamic_permit_100k_smoke.sh")


class PermitValidationTests(unittest.TestCase):
    def validate(self, changes=None, mode="1"):
        values = {
            "enableStarDynamicPermitHooks": "1",
            "dynamicPermitDelta.acquires": "5692",
            "dynamicPermitDelta.workUnits": "369546",
            "dynamicPermitDelta.feature.acquires": "5689",
            "dynamicPermitDelta.feature.workUnits": "200000",
            "dynamicPermitDelta.feature.waitNs": "1295269959",
            "dynamicPermitDelta.map.acquires": "3",
            "dynamicPermitDelta.map.workUnits": "169546",
            "dynamicPermitDelta.map.waitNs": "1844",
        }
        values.update(changes or {})
        with tempfile.TemporaryDirectory(prefix="pf_permit_validation_") as tmp:
            path = Path(tmp) / "api_run.txt"
            path.write_text("".join(f"{key}={value}\n" for key, value in values.items()
                                    if value is not None))
            env = dict(os.environ, STAR_BIN="/nonexistent/STAR",
                       PF_DYNAMIC_100K_MAX_FEATURE_WAIT_NS="20000000000")
            return subprocess.run(["bash", str(SCRIPT), "--validate-api-run",
                                   str(path), mode], env=env, text=True,
                                  capture_output=True, timeout=10)

    def test_batched_work_and_concurrent_map_are_valid(self):
        result = self.validate()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("STAR not executed", result.stdout)

    def test_uncontended_zero_wait_is_valid(self):
        result = self.validate({"dynamicPermitDelta.feature.waitNs": "0"})
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_counters_must_be_present_nonnegative_integers(self):
        for value in (None, "-1", "garbage", "1.5"):
            with self.subTest(value=value):
                self.assertNotEqual(self.validate({"dynamicPermitDelta.map.acquires": value}).returncode, 0)

    def test_feature_work_and_acquisitions_must_be_positive(self):
        for key in ("dynamicPermitDelta.feature.acquires", "dynamicPermitDelta.feature.workUnits",
                    "dynamicPermitDelta.acquires", "dynamicPermitDelta.workUnits"):
            with self.subTest(key=key):
                self.assertNotEqual(self.validate({key: "0"}).returncode, 0)

    def test_wait_limit_still_enforced(self):
        result = self.validate({"dynamicPermitDelta.feature.waitNs": "20000000001"})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("exceeds threshold", result.stderr)

    def test_dynamic_hooks_required(self):
        self.assertNotEqual(self.validate({"enableStarDynamicPermitHooks": "0"}).returncode, 0)

    def test_disabled_mode_rejects_feature_acquisitions(self):
        self.assertNotEqual(self.validate({"enableStarDynamicPermitHooks": "0"}, "0").returncode, 0)

    def test_disabled_mode_without_permit_work(self):
        result = self.validate({"enableStarDynamicPermitHooks": "0",
                                "dynamicPermitDelta.feature.acquires": "0",
                                "dynamicPermitDelta.map.acquires": "0"}, "0")
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_invalid_mode_rejected(self):
        self.assertNotEqual(self.validate(mode="2").returncode, 0)


if __name__ == "__main__":
    unittest.main()
