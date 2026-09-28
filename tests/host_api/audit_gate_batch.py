#!/usr/bin/env python3
"""Audit a run_gs1.sh batch without executing or modifying any workload.

This checks execution coverage, not biological output parity or performance.
Identical failures on two binaries remain failures, not passing gate evidence.
"""

import argparse
import csv
from datetime import datetime
import json
from pathlib import Path
import re


TIER_A = (
    "run_solo_smoke", "plain_scrna_exact_counts", "flex_gdna_removed",
    "run_scrna_sidecar_off_golden", "run_spatial_r1_tap_guard",
    "test_visium_hd_gex_sidecar_concurrency", "test_snp_mask_build_smoke",
    "run_flex_tiny_public_smoke", "run_molecule_first_native_smoke",
    "run_adapter_clip_synthetic_test", "run_transcriptvb_scatter_gather_smoke",
    "run_trim_qc_merge_smoke", "run_star_trim_qc_smoke",
)


def table(path):
    with path.open() as stream:
        return [row for row in csv.reader(stream, delimiter="\t")
                if row and not row[0].startswith("#")]


def read_status(path, expected, issues):
    statuses = {}
    try:
        rows = table(path)
    except OSError as exc:
        issues.append(f"Cannot read {path.name}: {exc}")
        rows = []
    for row in rows:
        if len(row) < 2:
            issues.append(f"Malformed status in {path.name}: {row!r}")
            continue
        name = row[0]
        if name in statuses:
            issues.append(f"Duplicate status in {path.name}: {name}")
            continue
        try:
            statuses[name] = int(row[1])
        except ValueError:
            issues.append(f"Invalid exit status for {name}: {row[1]!r}")
    for name in sorted(set(expected) - statuses.keys()):
        issues.append(f"Missing status in {path.name}: {name}")
    for name in sorted(statuses.keys() - set(expected)):
        issues.append(f"Unexpected status in {path.name}: {name}")
    return statuses


def audit(run, manifest):
    issues = []
    rows = table(manifest)
    cases = [row[1] for row in rows if row[0] != "multiome"]
    if not cases or len(cases) != len(set(cases)):
        raise ValueError("Manifest must contain distinct non-multiome cases")
    times = {}
    for marker in ("started_utc", "finished_utc"):
        try:
            stamp = datetime.fromisoformat(
                (run / marker).read_text().strip().replace("Z", "+00:00"))
            if stamp.tzinfo is None:
                raise ValueError("timestamp must include a timezone")
            times[marker] = stamp
        except (OSError, ValueError) as exc:
            issues.append(f"Missing or invalid {marker}: {exc}")
    if len(times) == 2 and times["finished_utc"] < times["started_utc"]:
        issues.append("finished_utc precedes started_utc")

    production = read_status(run / "manifest_status.tsv", cases, issues)
    tier_a = read_status(run / "tierA/status.tsv", TIER_A, issues)
    skipped = []
    for case in cases:
        path = run / "manifest" / case / "summary.tsv"
        try:
            records = table(path)
        except OSError as exc:
            issues.append(f"Missing summary for {case}: {exc}")
            continue
        terminal = [row for row in records if row[0] in ("PASS", "FAIL", "SKIP")]
        if len(terminal) != 1 or len(terminal[0]) < 3 or terminal[0][2] != case:
            issues.append(f"Missing, duplicated, or mismatched terminal summary: {case}")
            continue
        state = terminal[0][0]
        if state == "SKIP":
            skipped.append(case)
        elif (production.get(case) == 0) != (state == "PASS"):
            issues.append(f"Status/summary disagreement: {case}")

    failures = {
        "production": {k: v for k, v in production.items() if v != 0},
        "tier_a": {k: v for k, v in tier_a.items() if v != 0},
    }
    # A successful wrapper can deliberately continue after CellBender fails.
    # Its failure marker is authoritative for downstream coverage.
    downstream_failures = sorted(str(p.relative_to(run))
                                 for p in run.rglob("CELLBENDER_FAILED.txt"))
    nested_skips = []
    logs = sorted(set(run.glob("manifest/*/*/*out.log"))
                  | set(run.glob("manifest/*/*/*err.log"))
                  | set(run.glob("tierA/*/*out.log"))
                  | set(run.glob("tierA/*/*err.log")))
    for log in logs:
        with log.open(errors="replace") as stream:
            for number, line in enumerate(stream, 1):
                if re.match(r"^\s*(?:SKIP(?:\s|:|$)|\[SKIP\])", line):
                    nested_skips.append({"path": str(log.relative_to(run)),
                                         "line": number, "message": line.strip()})
    passed = not (issues or skipped or nested_skips or downstream_failures
                  or failures["production"] or failures["tier_a"])
    return {
        "schema": "star-host-api-gs1-execution-audit/v1",
        "run": str(run), "manifest": str(manifest),
        "execution_coverage_passed": passed,
        "output_parity": "not_evaluated",
        "performance": "not_evaluated",
        "expected": {"production": len(cases), "tier_a": len(TIER_A)},
        "recorded": {"production": len(production), "tier_a": len(tier_a)},
        "failures": failures, "skipped": skipped,
        "nested_skips_to_review": nested_skips,
        "downstream_failure_markers": downstream_failures,
        "issues": issues,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path,
                        help="New JSON report; existing reports are never overwritten")
    args = parser.parse_args()
    try:
        report = audit(args.run, args.manifest)
        payload = json.dumps(report, indent=2) + "\n"
        if args.report:
            with args.report.open("x") as stream:
                stream.write(payload)
    except (OSError, ValueError, IndexError) as exc:
        parser.exit(2, f"ERROR: {exc}\n")
    print(payload, end="")
    return 0 if report["execution_coverage_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
