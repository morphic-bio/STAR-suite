#!/usr/bin/env python3
"""Compare STAR builds on single-file and multi-file runs.

Every scenario runs once per build, serially, from the same work path, so the
command lines recorded in outputs match. The baseline runs with a high open-file
limit so that batch runs can finish. The candidate and sanitizer builds run at
the normal limit of 1024.

The test data are the six paired FASTQs of the public nf-core GSE110004 test set
and a STAR index built with GeneCounts annotation for its small genome.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

SAMPLES = ("SRR6357070", "SRR6357071", "SRR6357072", "SRR6357076", "SRR6357077", "SRR6357078")
SORTED = ("BAM", "SortedByCoordinate")
UNSORTED = ("BAM", "Unsorted")
BYSJOUT = ("--outFilterType", "BySJout", "--outSAMunmapped", "Within")
SLAM = ("--slamQuantMode", "1", "--slamStrandness", "Sense")
# The default sort limit follows genome size, which is too low for six merged files here.
LARGE_SORT_RAM = ("--limitBAMsortRAM", "4000000000")

SCENARIOS = {
    "single_sorted": {"mode": "single", "samples": 1},
    "single_unsorted": {"mode": "single", "samples": 1, "sam": UNSORTED},
    "single_twopass": {"mode": "single", "samples": 1, "extra": ("--twopassMode", "Basic")},
    "single_bysjout": {"mode": "single", "samples": 1, "extra": BYSJOUT},
    "single_batch": {"mode": "batch", "samples": 1},
    "single_slam": {"mode": "single", "samples": 1, "extra": SLAM},
    "merged_6_sorted": {"mode": "merged", "samples": 6, "extra": LARGE_SORT_RAM},
    "merged_6_bysjout": {"mode": "merged", "samples": 6, "extra": BYSJOUT + LARGE_SORT_RAM},
    "batch_2_sorted": {"mode": "batch", "samples": 2},
    "batch_6_sorted": {"mode": "batch", "samples": 6},
    "batch_6_sorted_16t": {"mode": "batch", "samples": 6, "threads": 16},
    "batch_6_unmapped_within": {"mode": "batch", "samples": 6, "extra": ("--outSAMunmapped", "Within")},
    "batch_6_unsorted": {"mode": "batch", "samples": 6, "sam": UNSORTED},
    "batch_3_slam": {"mode": "batch", "samples": 3, "extra": SLAM},
}
# The baseline also runs these at the normal limit, to show the failure being fixed.
LIMIT_CONTROLS = ("batch_6_sorted", "batch_6_sorted_16t")
IGNORED_SUFFIXES = ("Log.out", "Log.progress.out", "Log.std.out")
TIMING_LINES = (
    b"Started job on", b"Started mapping on", b"Finished on", b"Mapping speed",
    b"total ms", b"us/call", b"thread-seconds", b"MB/s", b"(ms)",
)
HIGH_LIMIT = 65536
NORMAL_LIMIT = 1024


def star_args(scenario: dict, reads: Path, genome: Path, work: Path, read_limit: int) -> list[str]:
    names = SAMPLES[: scenario["samples"]]
    mate1 = ",".join(str(reads / f"{name}_1.fastq.gz") for name in names)
    mate2 = ",".join(str(reads / f"{name}_2.fastq.gz") for name in names)
    args = [
        "--runMode", "alignReads",
        "--genomeDir", str(genome),
        "--readFilesIn", mate1, mate2,
        "--readFilesCommand", "zcat",
        "--outFileNamePrefix", f"{work}/",
        "--runThreadN", str(scenario.get("threads", 4)),
        "--quantMode", "GeneCounts",
        "--outSAMtype", *scenario.get("sam", SORTED),
    ]
    if scenario["mode"] == "batch":
        args += ["--batchMode", "1", "--outFileNamePrefixAuto", "1"]
    if read_limit:
        args += ["--readMapNumber", str(read_limit)]
    return args + list(scenario.get("extra", ()))


def watch_open_files(pid: int, work: Path, peaks: dict[int, int], stop: threading.Event) -> None:
    alignments = work / "alignments"
    while not stop.is_set():
        try:
            count = len(os.listdir(f"/proc/{pid}/fd"))
        except OSError:
            return
        try:
            started = len(os.listdir(alignments))
        except OSError:
            started = 0
        peaks[started] = max(peaks.get(started, 0), count)
        stop.wait(0.01)


def run_star(star: Path, args: list[str], work: Path, log: Path, limit: int, env: dict[str, str]) -> dict:
    def set_limit() -> None:
        resource.setrlimit(resource.RLIMIT_NOFILE, (limit, limit))

    work.mkdir(parents=True)
    peaks: dict[int, int] = {}
    stop = threading.Event()
    start = time.monotonic()
    with log.open("wb") as handle:
        process = subprocess.Popen(
            [str(star), *args], cwd=work, stdout=handle, stderr=subprocess.STDOUT, env=env, preexec_fn=set_limit
        )
        watcher = threading.Thread(target=watch_open_files, args=(process.pid, work, peaks, stop), daemon=True)
        watcher.start()
        code = process.wait()
    stop.set()
    watcher.join()
    return {
        "exit": code,
        "seconds": round(time.monotonic() - start, 1),
        "nofile": limit,
        "peak_open_files": max(peaks.values(), default=0),
        "peak_by_samples_started": {str(key): value for key, value in sorted(peaks.items())},
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def final_log_digest(path: Path) -> str:
    kept = [line for line in path.read_bytes().splitlines() if not any(mark in line for mark in TIMING_LINES)]
    return hashlib.sha256(b"\n".join(kept)).hexdigest()


def bam_digest(path: Path) -> dict:
    """Hash header and records independently of record order, and check coordinate order."""
    # --no-PG keeps samtools from adding its own command line, which holds this file's path.
    header = subprocess.run(["samtools", "view", "--no-PG", "-H", str(path)], check=True, capture_output=True).stdout
    coordinate = b"SO:coordinate" in header.splitlines()[0] if header else False
    rank = {}
    for line in header.splitlines():
        if line.startswith(b"@SQ"):
            name = next(field[3:] for field in line.split(b"\t") if field.startswith(b"SN:"))
            rank[name] = len(rank)
    count = total = 0
    in_order = True
    previous = (-1, -1)
    process = subprocess.Popen(["samtools", "view", str(path)], stdout=subprocess.PIPE)
    for line in process.stdout:
        count += 1
        value = int.from_bytes(hashlib.blake2b(line, digest_size=16).digest(), "big")
        total = (total + value) & ((1 << 128) - 1)
        if coordinate:
            fields = line.split(b"\t", 4)
            key = (rank.get(fields[2], len(rank)), int(fields[3]))
            in_order = in_order and key >= previous
            previous = key
    return {
        "readable": process.wait() == 0,
        "header": hashlib.sha256(header).hexdigest(),
        "records": count,
        "record_set": f"{total:032x}",
        "coordinate_order_ok": in_order if coordinate else None,
    }


def fingerprint(work: Path) -> dict:
    files: dict[str, object] = {}
    bam_bytes: dict[str, str] = {}
    for path in sorted(item for item in work.rglob("*") if item.is_file()):
        relative = path.relative_to(work)
        if any(part.startswith("_STARtmp") for part in relative.parts):
            continue
        if path.name.endswith(IGNORED_SUFFIXES):
            continue
        key = str(relative)
        if path.name.endswith("Log.final.out"):
            files[key] = final_log_digest(path)
        elif path.suffix == ".bam" and path.stat().st_size:
            files[key] = bam_digest(path)
            bam_bytes[key] = sha256_file(path)
        else:
            files[key] = sha256_file(path)
    return {"files": files, "bam_bytes": bam_bytes}


def differences(left: dict, right: dict) -> list[str]:
    a, b = left["files"], right["files"]
    return sorted(key for key in set(a) | set(b) if a.get(key) != b.get(key))


def run_build(label: str, binary: Path, name: str, limit: int, options, env: dict[str, str]) -> dict:
    run_dir = options.out / "run"
    work = run_dir / "work"
    star = run_dir / "STAR"
    if star.is_symlink():
        star.unlink()
    star.symlink_to(binary.resolve())
    final = options.out / label / name
    if final.exists():
        raise SystemExit(f"output already exists: {final}")
    final.parent.mkdir(parents=True, exist_ok=True)
    log = options.out / label / f"{name}.log"
    args = star_args(SCENARIOS[name], options.reads_dir, options.genome_dir, work, options.read_limit)
    result = run_star(star, args, work, log, limit, env)
    text = log.read_bytes()
    result["sanitizer_error"] = b"ERROR: AddressSanitizer" in text
    if result["exit"] == 0:
        result.update(fingerprint(work))
    else:
        tail = [line for line in text.decode(errors="replace").splitlines() if line.strip()]
        result["error_tail"] = tail[-4:]
    shutil.move(str(work), str(final))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--asan", type=Path)
    parser.add_argument("--genome-dir", type=Path, required=True)
    parser.add_argument("--reads-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--scenario", action="append", choices=sorted(SCENARIOS))
    parser.add_argument("--read-limit", type=int, default=0, help="Map only this many reads, for a harness check.")
    options = parser.parse_args()
    if options.out.exists():
        raise SystemExit(f"use a fresh output directory: {options.out}")
    (options.out / "run").mkdir(parents=True)

    plain_env = dict(os.environ)
    asan_env = {**plain_env, "ASAN_OPTIONS": "detect_leaks=0:halt_on_error=1:exitcode=99"}
    names = options.scenario or list(SCENARIOS)
    report: dict[str, dict] = {}
    failed = False
    for name in names:
        row = {"baseline": run_build("baseline", options.baseline, name, HIGH_LIMIT, options, plain_env)}
        if name in LIMIT_CONTROLS and options.candidate:
            row["baseline_at_limit"] = run_build(
                "baseline_at_limit", options.baseline, name, NORMAL_LIMIT, options, plain_env
            )
        checks = [("candidate", options.candidate, plain_env), ("asan", options.asan, asan_env)]
        for label, binary, env in checks:
            if binary:
                row[label] = run_build(label, binary, name, NORMAL_LIMIT, options, env)

        problems = []
        base = row["baseline"]
        if base["exit"] != 0:
            problems.append(f"baseline exit {base['exit']}")
        for label in ("candidate", "asan"):
            result = row.get(label)
            if result is None:
                continue
            if result["exit"] != 0:
                problems.append(f"{label} exit {result['exit']}")
            elif base["exit"] == 0:
                changed = differences(base, result)
                result["changed_files"] = changed
                if changed:
                    problems.append(f"{label} changed {len(changed)} files")
            if result["sanitizer_error"]:
                problems.append(f"{label} sanitizer error")
            for digest in (result.get("files") or {}).values():
                if isinstance(digest, dict) and not digest["readable"]:
                    problems.append(f"{label} BAM unreadable")
                elif isinstance(digest, dict) and digest["coordinate_order_ok"] is False:
                    problems.append(f"{label} BAM out of coordinate order")
        row["problems"] = problems
        failed = failed or bool(problems)
        report[name] = row
        peak = {label: row[label]["peak_open_files"] for label in row if label != "problems"}
        print(f"{name}: {'PASS' if not problems else 'FAIL ' + '; '.join(problems)} peak_open={peak}", flush=True)
        (options.out / "results.json").write_text(json.dumps(report, indent=1, sort_keys=True) + "\n")
    shutil.rmtree(options.out / "run")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
