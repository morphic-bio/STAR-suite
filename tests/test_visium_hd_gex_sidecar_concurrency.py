#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
WRAPPER = ROOT / "scripts" / "run_visium_hd_gex_sidecar_100k.py"
SPEC = importlib.util.spec_from_file_location("visium_hd_gex_sidecar_wrapper", WRAPPER)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import {WRAPPER}")
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def python_command(source: str, *arguments: Path) -> list[str]:
    return [sys.executable, "-c", source, *(str(path) for path in arguments)]


class ProducerThreadBudgetTest(unittest.TestCase):
    def test_serial_is_the_validated_default(self) -> None:
        with mock.patch.object(sys, "argv", [str(WRAPPER)]):
            arguments = MODULE.arguments()
            self.assertEqual(arguments.producer_mode, "serial")
            self.assertEqual(arguments.evidence_mode, "contracts")
            self.assertEqual(arguments.assignment_policy, "all")

    def test_fused_mode_requires_concurrent_producers(self) -> None:
        with mock.patch.object(
            sys, "argv", [str(WRAPPER), "--evidence-mode", "fused"]
        ), self.assertRaises(SystemExit):
            MODULE.arguments()
        with mock.patch.object(
            sys,
            "argv",
            [str(WRAPPER), "--evidence-mode", "fused", "--producer-mode", "concurrent"],
        ):
            arguments = MODULE.arguments()
            self.assertEqual(arguments.evidence_mode, "fused")
            self.assertEqual((arguments.r1_threads, arguments.star_threads), (8, 8))

    def test_concurrent_defaults_split_total_budget(self) -> None:
        self.assertEqual(MODULE.producer_thread_budgets(16, "concurrent", None, None), (8, 8))
        self.assertEqual(MODULE.producer_thread_budgets(16, "concurrent", 6, None), (6, 10))
        self.assertEqual(MODULE.producer_thread_budgets(16, "concurrent", None, 12), (4, 12))

    def test_serial_defaults_reuse_total_budget(self) -> None:
        self.assertEqual(MODULE.producer_thread_budgets(16, "serial", None, None), (16, 16))

    def test_invalid_budgets_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            MODULE.producer_thread_budgets(1, "concurrent", None, None)
        with self.assertRaises(ValueError):
            MODULE.producer_thread_budgets(16, "concurrent", 9, 8)
        with self.assertRaises(ValueError):
            MODULE.producer_thread_budgets(16, "serial", 17, 1)


class ProducerSchedulerTest(unittest.TestCase):
    def driver(self, root: Path) -> object:
        (root / "logs").mkdir()
        return MODULE.Driver(root)

    def test_concurrent_fork_join_preserves_dependency_barrier(self) -> None:
        with tempfile.TemporaryDirectory(prefix="star_gex_producer_concurrency_") as value:
            root = Path(value)
            driver = self.driver(root)
            star_started = root / "star.started"
            r1_done = root / "r1.done"
            h0_done = root / "h0.done"
            star_done = root / "star.done"
            r1 = python_command(
                """
import sys, time
from pathlib import Path
star_started, r1_done = map(Path, sys.argv[1:])
deadline = time.monotonic() + 5
while not star_started.exists():
    if time.monotonic() > deadline:
        raise SystemExit('STAR did not overlap the R1 producer')
    time.sleep(0.01)
r1_done.write_text('complete')
""",
                star_started,
                r1_done,
            )
            star = python_command(
                """
import sys, time
from pathlib import Path
star_started, r1_done, star_done = map(Path, sys.argv[1:])
star_started.write_text('running')
deadline = time.monotonic() + 5
while not r1_done.exists():
    if time.monotonic() > deadline:
        raise SystemExit('R1 did not overlap the STAR producer')
    time.sleep(0.01)
star_done.write_text('complete')
""",
                star_started,
                r1_done,
                star_done,
            )
            h0 = python_command(
                """
import sys
from pathlib import Path
r1_done, h0_done = map(Path, sys.argv[1:])
if not r1_done.exists():
    raise SystemExit('H0 prior started before R1 completion')
h0_done.write_text('complete')
""",
                r1_done,
                h0_done,
            )
            elapsed = MODULE.run_producer_stages(driver, "concurrent", r1, h0, star)
            self.assertGreaterEqual(elapsed, 0)
            self.assertTrue(h0_done.is_file())
            self.assertTrue(star_done.is_file())
            manifest = json.loads((root / "commands.json").read_text())
            self.assertEqual([row["label"] for row in manifest], ["r1_decode", "star_sidecar", "h0_prior"])
            self.assertTrue(all(row["state"] == "complete" for row in manifest))
            records = {row["label"]: row for row in manifest}
            self.assertLess(
                records["r1_decode"]["started_unix_seconds"],
                records["star_sidecar"]["finished_unix_seconds"],
            )
            self.assertLess(
                records["star_sidecar"]["started_unix_seconds"],
                records["r1_decode"]["finished_unix_seconds"],
            )
            self.assertGreaterEqual(
                records["h0_prior"]["started_unix_seconds"],
                records["r1_decode"]["finished_unix_seconds"],
            )

    def test_serial_control_keeps_historical_order(self) -> None:
        with tempfile.TemporaryDirectory(prefix="star_gex_producer_serial_") as value:
            root = Path(value)
            driver = self.driver(root)
            trace = root / "trace.txt"

            def append(label: str) -> list[str]:
                return python_command(
                    """
import sys
from pathlib import Path
with Path(sys.argv[1]).open('a') as handle:
    handle.write(sys.argv[2] + '\\n')
""",
                    trace,
                ) + [label]

            MODULE.run_producer_stages(
                driver, "serial", append("r1"), append("h0"), append("star")
            )
            self.assertEqual(trace.read_text().splitlines(), ["r1", "h0", "star"])

    def test_failure_terminates_sibling_and_suppresses_h0(self) -> None:
        with tempfile.TemporaryDirectory(prefix="star_gex_producer_failure_") as value:
            root = Path(value)
            driver = self.driver(root)
            r1_started = root / "r1.started"
            forbidden_h0 = root / "h0.should_not_exist"
            r1 = python_command(
                """
import sys, time
from pathlib import Path
Path(sys.argv[1]).write_text('running')
time.sleep(30)
""",
                r1_started,
            )
            star = python_command(
                """
import sys, time
from pathlib import Path
started = Path(sys.argv[1])
deadline = time.monotonic() + 5
while not started.exists():
    if time.monotonic() > deadline:
        raise SystemExit(8)
    time.sleep(0.01)
raise SystemExit(7)
""",
                r1_started,
            )
            h0 = python_command("from pathlib import Path; import sys; Path(sys.argv[1]).touch()", forbidden_h0)
            started = time.monotonic()
            with self.assertRaisesRegex(RuntimeError, "star_sidecar failed with exit code 7"):
                MODULE.run_producer_stages(driver, "concurrent", r1, h0, star)
            self.assertLess(time.monotonic() - started, 5)
            self.assertFalse(forbidden_h0.exists())
            manifest = json.loads((root / "commands.json").read_text())
            records = {row["label"]: row for row in manifest}
            self.assertEqual(records["star_sidecar"]["exit_code"], 7)
            self.assertEqual(records["star_sidecar"]["state"], "failed")
            self.assertTrue(records["r1_decode"]["terminated_by_driver"])
            self.assertEqual(records["r1_decode"]["state"], "terminated")


if __name__ == "__main__":
    unittest.main()
