#!/usr/bin/env python3
"""Validate the pinned official recipe and provenance release snapshots."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SHARE = ROOT / "share/star-suite"
METADATA = SHARE / "SNAPSHOTS.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def main() -> int:
    try:
        metadata = json.loads(METADATA.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"invalid snapshot metadata: {exc}")
    if metadata.get("schema") != "starsuite.release_snapshots/v1":
        fail("unexpected snapshot metadata schema")

    recipes = metadata.get("recipes", {})
    catalog_path = SHARE / recipes.get("catalog", "")
    if not catalog_path.is_file():
        fail(f"official catalog missing: {catalog_path}")
    if sha256(catalog_path) != recipes.get("catalog_sha256"):
        fail("official catalog digest does not match SNAPSHOTS.json")
    catalog_text = catalog_path.read_text()
    workflow_count = sum(
        1 for line in catalog_text.splitlines() if line.startswith("  - id: ")
    )
    if workflow_count != recipes.get("workflow_count"):
        fail(
            f"official catalog workflow count is {workflow_count}; "
            f"expected {recipes.get('workflow_count')}"
        )

    provenance = metadata.get("provenance", {})
    schema_path = SHARE / provenance.get("record_schema", "")
    if not schema_path.is_file():
        fail(f"official provenance schema missing: {schema_path}")
    if sha256(schema_path) != provenance.get("record_schema_sha256"):
        fail("official provenance schema digest does not match SNAPSHOTS.json")
    record_root = SHARE / "evidence/official/records"
    record_count = len(list(record_root.rglob("*.json")))
    if record_count != provenance.get("record_count"):
        fail(
            f"official provenance record count is {record_count}; "
            f"expected {provenance.get('record_count')}"
        )

    for snapshot_root in (catalog_path.parent, SHARE / "evidence/official"):
        if any(path.name == ".git" for path in snapshot_root.rglob(".git")):
            fail(f"snapshot contains Git internals: {snapshot_root}")
    print(
        f"PASS: {workflow_count} official recipes and {record_count} public evidence records"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
