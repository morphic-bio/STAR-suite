#!/usr/bin/env python3
from __future__ import annotations

import copy
import importlib.util
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "partition_manifest", ROOT / "scripts/validate_partition_manifest.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PartitionManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.fixture = json.loads(
            (ROOT / "tests/fixtures/partition_manifest.json").read_text()
        )

    def test_valid_contiguous_manifest(self) -> None:
        parts = MODULE.validate_manifest(copy.deepcopy(self.fixture))
        self.assertEqual([item["ordinal"] for item in parts], [0, 1])

    def test_gap_is_rejected(self) -> None:
        document = copy.deepcopy(self.fixture)
        document["partitions"][1]["logical_range"]["start"] = 101
        with self.assertRaisesRegex(MODULE.ManifestError, "contiguous start"):
            MODULE.validate_manifest(document)

    def test_reordered_ordinal_is_rejected(self) -> None:
        document = copy.deepcopy(self.fixture)
        document["partitions"][1]["ordinal"] = 3
        with self.assertRaisesRegex(MODULE.ManifestError, "contiguous"):
            MODULE.validate_manifest(document)

    def test_provider_details_are_not_part_of_public_contract(self) -> None:
        document = copy.deepcopy(self.fixture)
        document["provider"] = {"implementation": "private"}
        with self.assertRaisesRegex(MODULE.ManifestError, "unsupported fields"):
            MODULE.validate_manifest(document)


if __name__ == "__main__":
    unittest.main()

