#!/usr/bin/env python3
"""Parse STAR Log.out for VmRSS after Velocyto RAM markers (sorted-replay + integrated-hash spill).

Extracts VmRSS (kB) from the first ~1200 chars after each marker occurrence, then prints:
  MAX_VM_RSS_KB=<int>  — maximum across all markers in all given logs (backward compatible)
  N_SAMPLES=<int>      — total marker hits counted
  PER_LOG_MAX_VM_RSS_KB[<run_subdir>]=<int> — max within that log (parent directory name of Log.out)

Use PER_LOG lines to compare Stage 1 vs Stage 2 when global MAX is dominated by one run.

If a log has no Velocyto-specific markers (e.g. small GEX-only runs), fall back to scanning the
entire log for VmRSS lines and using those values so harnesses can still emit SUMMARY.txt.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

MARKERS = (
    "RAM after Velocyto sorted-replay materialization:",
    "RAM after Velocyto integrated-hash spill (all records staged to disk):",
)
RSS_RE = re.compile(r"VmRSS:\s*(\d+)\s*kB")


def rss_values_after_markers(log_text: str) -> list[int]:
    found: list[int] = []
    for marker in MARKERS:
        idx = 0
        while True:
            j = log_text.find(marker, idx)
            if j < 0:
                break
            chunk = log_text[j : j + 1200]
            m = RSS_RE.search(chunk)
            if m:
                found.append(int(m.group(1)))
            idx = j + len(marker)
    return found


def all_vm_rss_kb(log_text: str) -> list[int]:
    """Every VmRSS reading in the log (STAR prints /proc/self/status style lines)."""
    return [int(m.group(1)) for m in RSS_RE.finditer(log_text)]


def rss_values_for_log(log_text: str) -> tuple[list[int], bool]:
    """Return (values, used_marker_path)."""
    marker_vals = rss_values_after_markers(log_text)
    if marker_vals:
        return marker_vals, True
    fallback = all_vm_rss_kb(log_text)
    return fallback, False


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "logs",
        nargs="+",
        type=Path,
        help="STAR Log.out files (e.g. det_sort_t1/Log.out det_hash_t1/Log.out)",
    )
    args = p.parse_args()

    global_vals: list[int] = []
    per_log: list[tuple[str, list[int]]] = []
    any_fallback = False

    for path in args.logs:
        if not path.is_file():
            print(f"report_velocyto_sorted_replay_rss: missing {path}", file=sys.stderr)
            return 1
        text = path.read_text(encoding="utf-8", errors="replace")
        vals, used_markers = rss_values_for_log(text)
        if vals and not used_markers:
            any_fallback = True
        tag = path.parent.name if path.parent.name else path.name
        per_log.append((tag, vals))
        global_vals.extend(vals)

    if not global_vals:
        print(
            "report_velocyto_sorted_replay_rss: no Velocyto RAM markers and no VmRSS lines in logs",
            file=sys.stderr,
        )
        return 2

    if any_fallback:
        print(
            "report_velocyto_sorted_replay_rss: note: no Velocyto RAM markers in at least one log; "
            "using max VmRSS sampled from full Log.out (fallback)",
            file=sys.stderr,
        )

    method = "log_vmRSS_fallback" if any_fallback else "velocyto_markers"
    print(f"RSS_METHOD={method}")
    print(f"MAX_VM_RSS_KB={max(global_vals)}")
    print(f"N_SAMPLES={len(global_vals)}")
    for tag, vals in per_log:
        if not vals:
            print(f"PER_LOG_MAX_VM_RSS_KB[{tag}]=NA")
            continue
        print(f"PER_LOG_MAX_VM_RSS_KB[{tag}]={max(vals)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
