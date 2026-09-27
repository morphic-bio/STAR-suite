#!/usr/bin/env python3
"""Ensure real-data regression comparisons reject count, cell, and input drift."""
import copy
import unittest

from run_scrna_gex_100k_regression import verify_report


class ReportTests(unittest.TestCase):
    def setUp(self):
        self.expected = {"inputs": {"fixture": "pinned"}, "profiles": {
            "modern": {"qc": {"reads": 100000}, "matrices": {
                "Gene": {"raw": {"umis": 10, "counts_sha256": "counts"},
                         "filtered": {"barcodes_sha256": "cells"}}}}}}

    def test_equal_report(self):
        verify_report(copy.deepcopy(self.expected), self.expected, ["modern"])

    def test_changed_count(self):
        actual = copy.deepcopy(self.expected)
        actual["profiles"]["modern"]["matrices"]["Gene"]["raw"]["umis"] = 0
        with self.assertRaises(AssertionError):
            verify_report(actual, self.expected, ["modern"])

    def test_changed_cell_set(self):
        actual = copy.deepcopy(self.expected)
        actual["profiles"]["modern"]["matrices"]["Gene"]["filtered"]["barcodes_sha256"] = "different"
        with self.assertRaises(AssertionError):
            verify_report(actual, self.expected, ["modern"])

    def test_changed_input(self):
        actual = copy.deepcopy(self.expected)
        actual["inputs"]["fixture"] = "changed"
        with self.assertRaises(AssertionError):
            verify_report(actual, self.expected, ["modern"])

    def test_missing_profile(self):
        actual = copy.deepcopy(self.expected)
        actual["profiles"].clear()
        with self.assertRaises(KeyError):
            verify_report(actual, self.expected, ["modern"])


if __name__ == "__main__":
    unittest.main()
