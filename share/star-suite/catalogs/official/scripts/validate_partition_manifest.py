#!/usr/bin/env python3
"""Validate the public, provider-neutral partition-manifest contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any


SCHEMA = "biodepot.partition_manifest/v1"
ROOT_KEYS = {"schema", "dataset_id", "partition_set_id", "source_digest", "partitions"}
PARTITION_KEYS = {
    "id",
    "ordinal",
    "mate1",
    "mate2",
    "logical_range",
    "sha256_mate1",
    "sha256_mate2",
}
RANGE_KEYS = {"unit", "start", "end"}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class ManifestError(ValueError):
    pass


def _require_keys(value: dict[str, Any], allowed: set[str], label: str) -> None:
    unexpected = sorted(set(value) - allowed)
    if unexpected:
        raise ManifestError(f"{label} has unsupported fields: {', '.join(unexpected)}")


def _require_id(value: Any, label: str) -> str:
    if not isinstance(value, str) or not SAFE_ID.fullmatch(value):
        raise ManifestError(f"{label} must be a non-empty portable identifier")
    return value


def _require_path(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value or any(ch in value for ch in "\n\r\t"):
        raise ManifestError(f"{label} must be a non-empty path without control separators")
    return value


def _check_digest(value: Any, label: str) -> None:
    if value is not None and (not isinstance(value, str) or not SHA256.fullmatch(value)):
        raise ManifestError(f"{label} must be a lowercase SHA-256 digest")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_manifest(document: Any, *, check_files: bool = False) -> list[dict[str, Any]]:
    if not isinstance(document, dict):
        raise ManifestError("manifest root must be an object")
    _require_keys(document, ROOT_KEYS, "manifest")
    if document.get("schema") != SCHEMA:
        raise ManifestError(f"schema must be {SCHEMA}")
    _require_id(document.get("dataset_id"), "dataset_id")
    _require_id(document.get("partition_set_id"), "partition_set_id")
    _check_digest(document.get("source_digest"), "source_digest")

    partitions = document.get("partitions")
    if not isinstance(partitions, list) or not partitions:
        raise ManifestError("partitions must be a non-empty array")

    normalized: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    expected_start = 0
    for expected_ordinal, item in enumerate(partitions):
        label = f"partitions[{expected_ordinal}]"
        if not isinstance(item, dict):
            raise ManifestError(f"{label} must be an object")
        _require_keys(item, PARTITION_KEYS, label)
        part_id = _require_id(item.get("id"), f"{label}.id")
        if part_id in seen_ids:
            raise ManifestError(f"duplicate partition id: {part_id}")
        seen_ids.add(part_id)
        if item.get("ordinal") != expected_ordinal:
            raise ManifestError(
                f"{label}.ordinal must be contiguous and equal {expected_ordinal}"
            )

        mate1 = _require_path(item.get("mate1"), f"{label}.mate1")
        mate2 = _require_path(item.get("mate2"), f"{label}.mate2")
        logical_range = item.get("logical_range")
        if not isinstance(logical_range, dict):
            raise ManifestError(f"{label}.logical_range must be an object")
        _require_keys(logical_range, RANGE_KEYS, f"{label}.logical_range")
        if logical_range.get("unit") != "read_pairs":
            raise ManifestError(f"{label}.logical_range.unit must be read_pairs")
        start = logical_range.get("start")
        end = logical_range.get("end")
        if not isinstance(start, int) or isinstance(start, bool) or start < 0:
            raise ManifestError(f"{label}.logical_range.start must be a non-negative integer")
        if not isinstance(end, int) or isinstance(end, bool) or end <= start:
            raise ManifestError(f"{label}.logical_range.end must be greater than start")
        if start != expected_start:
            raise ManifestError(
                f"{label}.logical_range starts at {start}; expected contiguous start {expected_start}"
            )
        expected_start = end

        for key in ("sha256_mate1", "sha256_mate2"):
            _check_digest(item.get(key), f"{label}.{key}")
        if check_files:
            for key, raw_path in (("mate1", mate1), ("mate2", mate2)):
                path = Path(raw_path)
                if not path.is_file():
                    raise ManifestError(f"{label}.{key} is not a local file: {path}")
                digest_key = f"sha256_{key}"
                if item.get(digest_key) and _sha256(path) != item[digest_key]:
                    raise ManifestError(f"{label}.{digest_key} does not match {path}")
        normalized.append(item)
    return normalized


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--check-files", action="store_true")
    parser.add_argument(
        "--emit-tsv",
        action="store_true",
        help="Emit ordinal, id, mate1, mate2, start, and end after validation",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        with args.manifest.open(encoding="utf-8") as stream:
            document = json.load(stream)
        partitions = validate_manifest(document, check_files=args.check_files)
    except (OSError, json.JSONDecodeError, ManifestError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if args.emit_tsv:
        for item in partitions:
            logical_range = item["logical_range"]
            print(
                item["ordinal"],
                item["id"],
                item["mate1"],
                item["mate2"],
                logical_range["start"],
                logical_range["end"],
                sep="\t",
            )
    else:
        print(f"PASS: {args.manifest} ({len(partitions)} partitions)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
