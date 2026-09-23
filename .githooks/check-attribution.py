#!/usr/bin/env python3
"""Reject Anthropic co-author trailers in messages and pushed history."""
import re
import subprocess
import sys
from pathlib import Path

PATTERN = re.compile(rb"^[ \t]*Co-Authored-By:[^\n]*anthropic[^\n]*", re.M | re.I)

def main():
    mode = sys.argv[1]
    if mode == "message":
        if PATTERN.search(Path(sys.argv[2]).read_bytes()):
            sys.exit("Commit rejected: remove the Anthropic co-author trailer.")
        return
    if mode != "push":
        sys.exit("Expected message or push mode")
    tips = set()
    for line in sys.stdin:
        fields = line.split()
        if len(fields) != 4:
            sys.exit("Malformed pre-push input")
        oid = fields[1]
        if set(oid) == {"0"}:
            continue
        if not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", oid):
            sys.exit("Invalid pushed object ID")
        tips.add(oid)
    if not tips:
        return
    result = subprocess.run(
        ["git", "log", "--no-show-signature", "--format=%H%x00%B%x00", "--stdin"],
        input=("\n".join(sorted(tips)) + "\n").encode(),
        stdout=subprocess.PIPE, check=True,
    )
    fields = result.stdout.split(b"\0")
    for pos in range(0, len(fields) - 1, 2):
        if PATTERN.search(fields[pos + 1]):
            oid = fields[pos].strip().decode()
            sys.exit("Push rejected: commit " + oid + " retains an Anthropic co-author trailer. "
                     "Move this branch onto the rewritten history; see docs/HISTORY_REWRITE_20260913.md.")

if __name__ == "__main__":
    main()
