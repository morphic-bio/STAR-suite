#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import mock_consumer_report as mcr  # noqa: E402


class MockConsumerReportTests(unittest.TestCase):
    def test_parse_interface_and_telemetry(self) -> None:
        log_text = """
        Dynamic thread interface enabled: map permits=2 (runThreadN=4, telemetry=on)
        Dynamic thread telemetry: acquires=16, workUnits=128, workBytes=12345, waitMs(total/avg/max)=9.5/0.59375/1.3, workMs(total/avg/max)=100.1/6.25625/11.9
        """

        report = mcr.parse_log_text(log_text)
        self.assertEqual(report["interface"]["map_permits"], 2)
        self.assertEqual(report["interface"]["run_thread_n"], 4)
        self.assertTrue(report["interface"]["telemetry_enabled"])
        self.assertFalse(report["interface"]["variable_threads_enabled"])
        self.assertEqual(report["interface"]["retune_every_acquires"], 0)
        self.assertEqual(report["interface"]["retune_sequence_length"], 0)
        self.assertEqual(report["telemetry"]["acquires"], 16)
        self.assertEqual(report["telemetry"]["retunes"], 0)
        self.assertEqual(report["telemetry"]["retune_every_acquires"], 0)
        self.assertEqual(report["telemetry"]["retune_sequence_length"], 0)
        self.assertIsNone(report["telemetry"]["target_permits"])
        self.assertIsNone(report["telemetry"]["configured_permits"])
        self.assertEqual(report["telemetry"]["retune_trace"], [])
        self.assertEqual(report["telemetry"]["retune_trace_dropped"], 0)
        self.assertEqual(report["telemetry"]["work_units"], 128)
        self.assertAlmostEqual(report["telemetry"]["wait_ms_avg"], 0.59375)

    def test_parse_new_interface_and_telemetry_fields(self) -> None:
        log_text = """
        Dynamic thread interface enabled: map permits=3 (runThreadN=4, telemetry=on, variableThreads=on, retuneEveryAcquires=1, retuneSequenceLength=3)
        Dynamic thread telemetry: acquires=10, retunes=2, retuneEveryAcquires=1, retuneSequenceLength=3, targetPermits=4, configuredPermits=4, retuneTrace=2|4, retuneTraceDropped=0, workUnits=1000, workBytes=2000, waitMs(total/avg/max)=1.0/0.1/0.3, workMs(total/avg/max)=20.0/2.0/5.0
        """

        report = mcr.parse_log_text(log_text)
        self.assertTrue(report["interface"]["variable_threads_enabled"])
        self.assertEqual(report["interface"]["retune_every_acquires"], 1)
        self.assertEqual(report["interface"]["retune_sequence_length"], 3)
        self.assertEqual(report["telemetry"]["retunes"], 2)
        self.assertEqual(report["telemetry"]["retune_every_acquires"], 1)
        self.assertEqual(report["telemetry"]["retune_sequence_length"], 3)
        self.assertEqual(report["telemetry"]["target_permits"], 4)
        self.assertEqual(report["telemetry"]["configured_permits"], 4)
        self.assertEqual(report["telemetry"]["retune_trace"], [2, 4])
        self.assertEqual(report["telemetry"]["retune_trace_dropped"], 0)

    def test_missing_telemetry_validation(self) -> None:
        report = mcr.parse_log_text("some unrelated line")
        errors = mcr.validate_report(
            report=report,
            require_telemetry=True,
            min_acquires=0,
            min_work_units=0,
            min_retunes=0,
        )
        self.assertEqual(errors, ["Telemetry line not found in log."])

    def test_threshold_validation(self) -> None:
        log_text = """
        Dynamic thread telemetry: acquires=1, workUnits=2, workBytes=3, waitMs(total/avg/max)=0.1/0.1/0.1, workMs(total/avg/max)=0.2/0.2/0.2
        """
        report = mcr.parse_log_text(log_text)
        errors = mcr.validate_report(
            report=report,
            require_telemetry=True,
            min_acquires=2,
            min_work_units=3,
            min_retunes=2,
        )
        self.assertEqual(
            errors,
            [
                "acquires=1 is below required minimum 2.",
                "retunes=0 is below required minimum 2.",
                "work_units=2 is below required minimum 3.",
            ],
        )


if __name__ == "__main__":
    unittest.main()
