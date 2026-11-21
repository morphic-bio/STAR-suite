#!/usr/bin/env python3
"""Check that SLAM --trimScope per-file keeps thread 0's temporary output.

Without automatic trimming, per-file and first-file scopes map the same reads the
same way. A baseline first-file run is therefore the reference for the alignments,
gene counts, junctions, unmapped reads and chimeric junctions of each per-file run.
One thread makes thread 0 map every file, which is the case that loses data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from batch_bam_handles_regression import NORMAL_LIMIT, SLAM, bam_digest, run_star, sha256_file, star_args  # noqa: E402

CASES = {
    "one_thread_two_files": {"files": 2, "threads": 1},
    "four_threads_two_files": {"files": 2, "threads": 4},
    "one_thread_two_files_auto_trim": {"files": 2, "threads": 1, "extra": ("--autoTrim", "variance"), "reference": False},
}
BAM = "Aligned.sortedByCoord.out.bam"
EXACT_FILES = ("ReadsPerGene.out.tab", "SJ.out.tab")
LINE_SET_FILES = ("Unmapped.out.mate1", "Unmapped.out.mate2", "Chimeric.out.junction")


def line_set_digest(path: Path) -> str:
    total = 0
    with path.open("rb") as handle:
        for line in handle:
            if not line.startswith(b"#"):
                value = int.from_bytes(hashlib.blake2b(line, digest_size=16).digest(), "big")
                total = (total + value) & ((1 << 128) - 1)
    return f"{total:032x}"


def final_log_value(path: Path, label: str) -> int:
    for line in path.read_text().splitlines():
        if label in line:
            return int(line.split("|")[1])
    raise ValueError(f"{label} not found in {path}")


def summarize(work: Path) -> dict:
    bam = bam_digest(work / BAM)
    summary = {"bam": bam, "files": {}}
    for name in EXACT_FILES:
        summary["files"][name] = sha256_file(work / name)
    for name in LINE_SET_FILES:
        if (work / name).is_file():
            summary["files"][name] = line_set_digest(work / name)
    primary = subprocess.run(
        ["samtools", "view", "-c", "-F", "256", str(work / BAM)], capture_output=True, text=True
    )
    log = work / "Log.final.out"
    summary["primary_records"] = int(primary.stdout.strip() or -1)
    summary["expected_primary_records"] = 2 * (
        final_log_value(log, "Uniquely mapped reads number") + final_log_value(log, "Number of reads mapped to multiple loci")
    )
    return summary


def run_case(label: str, binary: Path, name: str, scope: str, options, env: dict[str, str]) -> dict:
    case = CASES[name]
    run_dir = options.out / "run"
    star = run_dir / "STAR"
    work = run_dir / "work"
    if star.is_symlink():
        star.unlink()
    star.symlink_to(binary.resolve())
    scenario = {
        "mode": "merged",
        "samples": case["files"],
        "threads": case["threads"],
        "extra": SLAM + case.get("extra", ()) + ("--trimScope", scope),
    }
    final = options.out / label / name
    final.parent.mkdir(parents=True, exist_ok=True)
    log = options.out / label / f"{name}.log"
    args = star_args(scenario, options.reads_dir, options.genome_dir, work, 0)
    result = run_star(star, args, work, log, NORMAL_LIMIT, env)
    result["sanitizer_error"] = b"ERROR: AddressSanitizer" in log.read_bytes()
    if result["exit"] == 0:
        result.update(summarize(work))
    shutil.move(str(work), str(final))
    return result


def integrity_problems(label: str, result: dict) -> list[str]:
    if result["exit"] != 0:
        return [f"{label} exit {result['exit']}"]
    problems = []
    if result["sanitizer_error"]:
        problems.append(f"{label} sanitizer error")
    if not result["bam"]["readable"]:
        problems.append(f"{label} BAM unreadable")
    elif result["bam"]["coordinate_order_ok"] is False:
        problems.append(f"{label} BAM out of coordinate order")
    if result["primary_records"] != result["expected_primary_records"]:
        problems.append(f"{label} has {result['primary_records']} primary records, expected {result['expected_primary_records']}")
    return problems


def reference_problems(label: str, result: dict, reference: dict) -> list[str]:
    problems = []
    if (result["bam"]["records"], result["bam"]["record_set"]) != (reference["bam"]["records"], reference["bam"]["record_set"]):
        problems.append(f"{label} BAM records differ from reference")
    for name in sorted(set(result["files"]) | set(reference["files"])):
        if result["files"].get(name) != reference["files"].get(name):
            problems.append(f"{label} {name} differs from reference")
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--asan", type=Path)
    parser.add_argument("--genome-dir", type=Path, required=True)
    parser.add_argument("--reads-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--skip-baseline-per-file", action="append", default=[], choices=sorted(CASES))
    options = parser.parse_args()
    if options.out.exists():
        raise SystemExit(f"use a fresh output directory: {options.out}")
    (options.out / "run").mkdir(parents=True)

    plain = dict(os.environ)
    asan = {**plain, "ASAN_OPTIONS": "detect_leaks=0:halt_on_error=1:exitcode=99"}
    report = {}
    failed = False
    for name, case in CASES.items():
        row = {}
        problems = []
        reference = None
        if case.get("reference", True):
            reference = run_case("reference_first_scope", options.baseline, name, "first", options, plain)
            row["reference_first_scope"] = reference
            problems += integrity_problems("reference", reference)
        if name not in options.skip_baseline_per_file:
            baseline = run_case("baseline_per_file", options.baseline, name, "per-file", options, plain)
            row["baseline_per_file"] = baseline
            row["baseline_per_file_problems"] = integrity_problems("baseline", baseline)
            if reference and baseline["exit"] == 0 and baseline["bam"]["readable"]:
                row["baseline_per_file_problems"] += reference_problems("baseline", baseline, reference)
        checks = [("candidate", options.candidate, plain), ("asan", options.asan, asan)]
        for label, binary, env in checks:
            if not binary:
                continue
            result = run_case(label, binary, name, "per-file", options, env)
            row[label] = result
            found = integrity_problems(label, result)
            if not found and reference and reference["exit"] == 0:
                found += reference_problems(label, result, reference)
            problems += found
        row["problems"] = problems
        failed = failed or bool(problems)
        report[name] = row
        before = row.get("baseline_per_file_problems", "not rerun")
        print(f"{name}: {'PASS' if not problems else 'FAIL ' + '; '.join(problems)} | unpatched per-file: {before}", flush=True)
        (options.out / "results.json").write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    shutil.rmtree(options.out / "run")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
