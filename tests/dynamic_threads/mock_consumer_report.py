#!/usr/bin/env python3
"""Mock consumer report generator for STAR dynamic-thread telemetry logs."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List

INTERFACE_RE = re.compile(
    r"Dynamic thread interface enabled: map permits=(\d+) "
    r"\(runThreadN=(\d+), telemetry=(on|off)"
    r"(?:, variableThreads=(on|off))?"
    r"(?:, retuneEveryAcquires=(\d+), retuneSequenceLength=(\d+))?"
    r"\)"
)

TELEMETRY_RE = re.compile(
    r"Dynamic thread telemetry: acquires=(\d+)"
    r"(?:, retunes=(\d+), retuneEveryAcquires=(\d+), retuneSequenceLength=(\d+), targetPermits=(\d+), configuredPermits=(\d+))?"
    r"(?:, retuneTrace=([^,]+), retuneTraceDropped=(\d+))?"
    r", workUnits=(\d+), workBytes=(\d+), "
    r"waitMs\(total/avg/max\)=([0-9eE+\-\.]+)/([0-9eE+\-\.]+)/([0-9eE+\-\.]+), "
    r"workMs\(total/avg/max\)=([0-9eE+\-\.]+)/([0-9eE+\-\.]+)/([0-9eE+\-\.]+)"
)


def _parse_retune_trace(trace_raw: str | None) -> List[int]:
    if trace_raw is None:
        return []
    trace_raw = trace_raw.strip()
    if trace_raw == "" or trace_raw == "-":
        return []
    return [int(token) for token in trace_raw.split("|")]


def parse_log_text(log_text: str) -> Dict[str, Any]:
    report: Dict[str, Any] = {"interface": None, "telemetry": None}

    for line in log_text.splitlines():
        interface_match = INTERFACE_RE.search(line)
        if interface_match:
            report["interface"] = {
                "map_permits": int(interface_match.group(1)),
                "run_thread_n": int(interface_match.group(2)),
                "telemetry_enabled": interface_match.group(3) == "on",
                "variable_threads_enabled": interface_match.group(4) == "on",
                "retune_every_acquires": (
                    int(interface_match.group(5))
                    if interface_match.group(5) is not None
                    else 0
                ),
                "retune_sequence_length": (
                    int(interface_match.group(6))
                    if interface_match.group(6) is not None
                    else 0
                ),
            }

        telemetry_match = TELEMETRY_RE.search(line)
        if telemetry_match:
            retunes = (
                int(telemetry_match.group(2))
                if telemetry_match.group(2) is not None
                else 0
            )
            retune_every_acquires = (
                int(telemetry_match.group(3))
                if telemetry_match.group(3) is not None
                else 0
            )
            retune_sequence_length = (
                int(telemetry_match.group(4))
                if telemetry_match.group(4) is not None
                else 0
            )
            target_permits = (
                int(telemetry_match.group(5))
                if telemetry_match.group(5) is not None
                else None
            )
            configured_permits = (
                int(telemetry_match.group(6))
                if telemetry_match.group(6) is not None
                else None
            )
            retune_trace = _parse_retune_trace(telemetry_match.group(7))
            retune_trace_dropped = (
                int(telemetry_match.group(8))
                if telemetry_match.group(8) is not None
                else 0
            )
            report["telemetry"] = {
                "acquires": int(telemetry_match.group(1)),
                "retunes": retunes,
                "retune_every_acquires": retune_every_acquires,
                "retune_sequence_length": retune_sequence_length,
                "target_permits": target_permits,
                "configured_permits": configured_permits,
                "retune_trace": retune_trace,
                "retune_trace_dropped": retune_trace_dropped,
                "work_units": int(telemetry_match.group(9)),
                "work_bytes": int(telemetry_match.group(10)),
                "wait_ms_total": float(telemetry_match.group(11)),
                "wait_ms_avg": float(telemetry_match.group(12)),
                "wait_ms_max": float(telemetry_match.group(13)),
                "work_ms_total": float(telemetry_match.group(14)),
                "work_ms_avg": float(telemetry_match.group(15)),
                "work_ms_max": float(telemetry_match.group(16)),
            }

    return report


def validate_report(
    report: Dict[str, Any],
    require_telemetry: bool,
    min_acquires: int,
    min_work_units: int,
    min_retunes: int,
) -> List[str]:
    errors: List[str] = []
    telemetry = report.get("telemetry")

    if require_telemetry and telemetry is None:
        errors.append("Telemetry line not found in log.")
        return errors

    if telemetry is None:
        return errors

    if telemetry["acquires"] < min_acquires:
        errors.append(
            f"acquires={telemetry['acquires']} is below required minimum {min_acquires}."
        )

    if telemetry["retunes"] < min_retunes:
        errors.append(
            f"retunes={telemetry['retunes']} is below required minimum {min_retunes}."
        )

    if telemetry["work_units"] < min_work_units:
        errors.append(
            f"work_units={telemetry['work_units']} is below required minimum {min_work_units}."
        )

    return errors


def render_summary(report: Dict[str, Any]) -> str:
    lines: List[str] = []
    interface = report.get("interface")
    telemetry = report.get("telemetry")

    if interface is None:
        lines.append("interface: not found")
    else:
        lines.append(
            "interface: map_permits={map_permits}, run_thread_n={run_thread_n}, telemetry_enabled={telemetry_enabled}, variable_threads_enabled={variable_threads_enabled}, retune_every_acquires={retune_every_acquires}, retune_sequence_length={retune_sequence_length}".format(
                **interface
            )
        )

    if telemetry is None:
        lines.append("telemetry: not found")
    else:
        lines.append(
            "telemetry: acquires={acquires}, retunes={retunes}, retune_every_acquires={retune_every_acquires}, retune_sequence_length={retune_sequence_length}, target_permits={target_permits}, configured_permits={configured_permits}, retune_trace={retune_trace}, retune_trace_dropped={retune_trace_dropped}, work_units={work_units}, work_bytes={work_bytes}".format(
                **telemetry
            )
        )
        lines.append(
            "timing_ms: wait(total/avg/max)={wait_ms_total:.6f}/{wait_ms_avg:.6f}/{wait_ms_max:.6f}".format(
                **telemetry
            )
        )
        lines.append(
            "timing_ms: work(total/avg/max)={work_ms_total:.6f}/{work_ms_avg:.6f}/{work_ms_max:.6f}".format(
                **telemetry
            )
        )

    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a mock consumer report from STAR dynamic-thread telemetry lines."
    )
    parser.add_argument("--log", required=True, help="Path to STAR Log.out file.")
    parser.add_argument(
        "--json-out",
        default="",
        help="Optional JSON output path for machine-readable report.",
    )
    parser.add_argument(
        "--summary-out",
        default="",
        help="Optional plain-text summary output path.",
    )
    parser.add_argument(
        "--require-telemetry",
        action="store_true",
        help="Fail if the telemetry line is missing.",
    )
    parser.add_argument(
        "--min-acquires",
        type=int,
        default=0,
        help="Minimum expected telemetry acquires count.",
    )
    parser.add_argument(
        "--min-work-units",
        type=int,
        default=0,
        help="Minimum expected telemetry work_units count.",
    )
    parser.add_argument(
        "--min-retunes",
        type=int,
        default=0,
        help="Minimum expected telemetry retunes count.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    log_path = Path(args.log)
    if not log_path.is_file():
        print(f"ERROR: log file does not exist: {log_path}", file=sys.stderr)
        return 2

    report = parse_log_text(log_path.read_text(encoding="utf-8", errors="replace"))
    errors = validate_report(
        report=report,
        require_telemetry=args.require_telemetry,
        min_acquires=args.min_acquires,
        min_work_units=args.min_work_units,
        min_retunes=args.min_retunes,
    )

    summary_text = render_summary(report)
    print(summary_text, end="")

    if args.summary_out:
        Path(args.summary_out).write_text(summary_text, encoding="utf-8")

    if args.json_out:
        Path(args.json_out).write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
