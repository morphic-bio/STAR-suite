#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "validate_records", ROOT / "scripts/validate_records.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ValidatorTests(unittest.TestCase):
    def test_forbids_operational_field_names(self) -> None:
        self.assertIsNotNone(MODULE.FORBIDDEN_KEYS.search("transfer_endpoint"))
        self.assertIsNotNone(MODULE.FORBIDDEN_KEYS.search("login_host"))

    def test_allows_public_hardware_class(self) -> None:
        self.assertIsNone(MODULE.FORBIDDEN_KEYS.search("hardware_class"))

    def test_forbids_absolute_host_paths(self) -> None:
        pattern = MODULE.FORBIDDEN_VALUES["absolute host path"]
        self.assertIsNotNone(pattern.search("/storage/project/raw.fastq.gz"))
        self.assertIsNone(pattern.search("record://paper/example"))


if __name__ == "__main__":
    unittest.main()

