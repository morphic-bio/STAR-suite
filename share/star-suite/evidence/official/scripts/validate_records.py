#!/usr/bin/env python3
"""Validate public STAR Suite provenance records and disclosure boundaries."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any, Iterator

from jsonschema import Draft202012Validator, FormatChecker


ROOT = Path(__file__).resolve().parents[1]
SCHEMA_PATH = ROOT / "schema/record-v1.schema.json"
RECORD_ROOT = ROOT / "records"
FORBIDDEN_KEYS = re.compile(
    r"(?:^|_)(?:account|credential|endpoint|login_host|operator|provider|raw_path|task_id|token)(?:$|_)",
    re.IGNORECASE,
)
FORBIDDEN_VALUES = {
    "absolute host path": re.compile(r"/(?:mnt|storage|home)/[A-Za-z0-9_.-]"),
    "private key": re.compile(r"BEGIN (?:RSA |OPENSSH )?PRIVATE KEY"),
    "GitHub token": re.compile(r"gh[opsu]_[A-Za-z0-9]{20,}"),
    "AWS access key": re.compile(r"AKIA[0-9A-Z]{16}"),
    "transfer-style UUID": re.compile(
        r"\b[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}\b",
        re.IGNORECASE,
    ),
    "unpublished implementation": re.compile(r"in-place\s+shard", re.IGNORECASE),
}


def walk(value: Any, path: str = "$") -> Iterator[tuple[str, str | None, Any]]:
    if isinstance(value, dict):
        for key, child in value.items():
            child_path = f"{path}.{key}"
            yield child_path, key, child
            yield from walk(child, child_path)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            child_path = f"{path}[{index}]"
            yield child_path, None, child
            yield from walk(child, child_path)


def main() -> int:
    schema = json.loads(SCHEMA_PATH.read_text())
    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    failures: list[str] = []
    record_ids: dict[str, Path] = {}
    paths = sorted(RECORD_ROOT.rglob("*.json"))

    for path in paths:
        relative = path.relative_to(ROOT)
        try:
            record = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            failures.append(f"{relative}: invalid JSON: {exc}")
            continue
        for error in sorted(validator.iter_errors(record), key=lambda item: list(item.path)):
            location = ".".join(str(item) for item in error.path) or "$"
            failures.append(f"{relative}:{location}: {error.message}")

        record_id = record.get("record_id")
        if isinstance(record_id, str):
            if record_id in record_ids:
                failures.append(
                    f"{relative}: duplicate record_id also in {record_ids[record_id].relative_to(ROOT)}"
                )
            record_ids[record_id] = path

        for location, key, value in walk(record):
            if key is not None and FORBIDDEN_KEYS.search(key):
                failures.append(f"{relative}:{location}: forbidden public field name")
            if isinstance(value, str):
                for label, pattern in FORBIDDEN_VALUES.items():
                    if pattern.search(value):
                        failures.append(f"{relative}:{location}: contains {label}")

    if failures:
        print("ERROR: public provenance validation failed", file=sys.stderr)
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"PASS: validated {len(paths)} public provenance records")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

