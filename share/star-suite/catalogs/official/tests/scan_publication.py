#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKIP = {".git", "__pycache__", ".pytest_cache"}
PATTERNS = {
    "private key": re.compile(r"BEGIN (?:RSA |OPENSSH )?PRIVATE KEY"),
    "GitHub token": re.compile(r"gh[opsu]_[A-Za-z0-9]{20,}"),
    "AWS access key": re.compile(r"AKIA[0-9A-Z]{16}"),
    "workstation path": re.compile(r"/(?:mnt|storage|home)/[A-Za-z0-9_.-]"),
    "unpublished implementation": re.compile(r"in-place\s+shard", re.IGNORECASE),
}

failures: list[str] = []
for path in ROOT.rglob("*"):
    if not path.is_file() or any(part in SKIP for part in path.parts):
        continue
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        continue
    for label, pattern in PATTERNS.items():
        for match in pattern.finditer(text):
            line = text.count("\n", 0, match.start()) + 1
            failures.append(f"{path.relative_to(ROOT)}:{line}: {label}")

if failures:
    raise SystemExit("ERROR: public disclosure scan failed\n" + "\n".join(failures))
print("PASS: public disclosure scan")

