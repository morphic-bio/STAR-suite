#!/usr/bin/env python3
"""
Normalize LSAN output to comparable stack signatures.

Input: an ASAN/LSAN log (stdout/stderr capture).
Output: one normalized stack signature per leak record.

Goal: allow comparing two runs (e.g., stub on/off) for "no new leak stacks".
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


FRAME_RE = re.compile(r"^\s*#\d+\s+0x[0-9a-fA-F]+\s+in\s+(.+?)\s+\(")


def normalize_func(func: str) -> str:
    func = func.strip()
    # Drop template args to reduce incidental diffs.
    func = re.sub(r"<[^<>]*>", "<>", func)
    return func


def parse_signatures(text: str) -> list[str]:
    signatures: list[str] = []

    in_leak = False
    frames: list[str] = []

    def flush():
        nonlocal frames
        if frames:
            # Use first few frames; keep stable "who allocated" identity.
            sig = " | ".join(frames[:6])
            signatures.append(sig)
        frames = []

    for line in text.splitlines():
        if "ERROR: LeakSanitizer: detected memory leaks" in line:
            in_leak = True
            continue
        if not in_leak:
            continue
        if line.startswith("SUMMARY: AddressSanitizer:"):
            flush()
            break
        if line.startswith("Direct leak of ") or line.startswith("Indirect leak of "):
            flush()
            continue

        m = FRAME_RE.match(line)
        if m:
            frames.append(normalize_func(m.group(1)))

    flush()
    return signatures


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <asan.log>", file=sys.stderr)
        return 2
    path = Path(sys.argv[1])
    text = path.read_text(errors="replace")
    sigs = parse_signatures(text)
    for s in sigs:
        print(s)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

