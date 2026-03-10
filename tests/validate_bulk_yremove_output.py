#!/usr/bin/env python3
import argparse
import gzip
import os
import subprocess
import sys
from pathlib import Path


def fail(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(2)


def read_fastq_names(path: Path) -> set[str]:
    names: set[str] = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for idx, line in enumerate(handle):
            if idx % 4 != 0:
                continue
            header = line.strip()
            if not header.startswith("@"):
                fail(f"Unexpected FASTQ header in {path}: {header!r}")
            name = header[1:].split()[0]
            if name.endswith("/1") or name.endswith("/2"):
                name = name[:-2]
            names.add(name)
    return names


def read_bam_names(path: Path) -> set[str]:
    proc = subprocess.run(
        ["samtools", "view", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    names = set()
    for line in proc.stdout.splitlines():
        if not line:
            continue
        names.add(line.split("\t", 1)[0])
    return names


def locate_fastq_pair(root: Path, token: str) -> tuple[Path, Path]:
    candidates = []
    for search_root in (root, root / "y_separated"):
        if not search_root.exists():
            continue
        patterns = [
            (f"*_{token}_R1_001.fastq*", f"*_{token}_R2_001.fastq*"),
            (f"{token}_reads.mate1.fastq*", f"{token}_reads.mate2.fastq*"),
        ]
        for p1, p2 in patterns:
            r1 = sorted(search_root.glob(p1))
            r2 = sorted(search_root.glob(p2))
            if len(r1) == 1 and len(r2) == 1:
                candidates.append((r1[0], r2[0]))
    if len(candidates) != 1:
        fail(f"Could not uniquely locate FASTQ pair for token {token!r} under {root}")
    return candidates[0]


def locate_bam(root: Path, token: str) -> Path:
    candidates = [
        root / f"Aligned.sortedByCoord.out_{token}.bam",
        root / f"Aligned.out_{token}.bam",
    ]
    for path in candidates:
        if path.exists():
            return path
    fail(f"Could not locate BAM for token {token!r} under {root}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate BAM-to-FASTQ parity for Y/noY emission")
    parser.add_argument("--outdir", required=True, help="STAR output directory to validate")
    parser.add_argument("--require-y-reads", action="store_true", help="Fail if no Y BAM reads are present")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    if not outdir.is_dir():
        fail(f"Output directory not found: {outdir}")

    y_bam = locate_bam(outdir, "Y")
    noy_bam = locate_bam(outdir, "noY")
    y_r1, y_r2 = locate_fastq_pair(outdir, "Y")
    noy_r1, noy_r2 = locate_fastq_pair(outdir, "noY")

    print(f"Y BAM: {y_bam}")
    print(f"noY BAM: {noy_bam}")
    print(f"Y FASTQ R1/R2: {y_r1} | {y_r2}")
    print(f"noY FASTQ R1/R2: {noy_r1} | {noy_r2}")

    y_bam_names = read_bam_names(y_bam)
    noy_bam_names = read_bam_names(noy_bam)
    y_r1_names = read_fastq_names(y_r1)
    y_r2_names = read_fastq_names(y_r2)
    noy_r1_names = read_fastq_names(noy_r1)
    noy_r2_names = read_fastq_names(noy_r2)

    ok = True

    print(f"Y BAM unique names: {len(y_bam_names)}")
    print(f"noY BAM unique names: {len(noy_bam_names)}")
    print(f"Y FASTQ R1 unique names: {len(y_r1_names)}")
    print(f"Y FASTQ R2 unique names: {len(y_r2_names)}")
    print(f"noY FASTQ R1 unique names: {len(noy_r1_names)}")
    print(f"noY FASTQ R2 unique names: {len(noy_r2_names)}")

    if args.require_y_reads and not y_bam_names:
        print("FAIL: expected at least one Y BAM read", file=sys.stderr)
        ok = False

    if y_r1_names != y_r2_names:
        print(f"FAIL: Y FASTQ mate name sets differ by {len(y_r1_names ^ y_r2_names)}", file=sys.stderr)
        ok = False
    if noy_r1_names != noy_r2_names:
        print(f"FAIL: noY FASTQ mate name sets differ by {len(noy_r1_names ^ noy_r2_names)}", file=sys.stderr)
        ok = False

    y_missing = y_bam_names - y_r1_names
    noy_missing = noy_bam_names - noy_r1_names
    if y_missing:
        print(f"FAIL: Y BAM names missing from Y FASTQ: {len(y_missing)}", file=sys.stderr)
        ok = False
    if noy_missing:
        print(f"FAIL: noY BAM names missing from noY FASTQ: {len(noy_missing)}", file=sys.stderr)
        ok = False

    if y_bam_names & noy_r1_names:
        print(f"FAIL: Y BAM names leaked into noY FASTQ: {len(y_bam_names & noy_r1_names)}", file=sys.stderr)
        ok = False
    if noy_bam_names & y_r1_names:
        print(f"FAIL: noY BAM names leaked into Y FASTQ: {len(noy_bam_names & y_r1_names)}", file=sys.stderr)
        ok = False

    if y_r1_names & noy_r1_names:
        print(f"FAIL: Y and noY FASTQ outputs overlap: {len(y_r1_names & noy_r1_names)}", file=sys.stderr)
        ok = False

    if y_bam_names & noy_bam_names:
        print(f"FAIL: Y and noY BAM outputs overlap: {len(y_bam_names & noy_bam_names)}", file=sys.stderr)
        ok = False

    if ok:
        print("PASS: BAM-to-FASTQ Y/noY parity checks succeeded")
        return 0
    print("FAIL: BAM-to-FASTQ Y/noY parity checks failed", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
