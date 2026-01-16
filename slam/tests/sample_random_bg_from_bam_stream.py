#!/usr/bin/env python3
"""
Reservoir-sample random 1bp loci from a BAM stream (SAM records on stdin).

Motivation:
- A naive "head -n N" after coordinate-sorted BAM output biases strongly to chr1.
- This script performs deterministic reservoir sampling over unique (chrom,pos0)
  loci seen in the stream, optionally excluding a provided BED set.

Usage:
  samtools view -F 4 -s 1.05 in.bam | \\
    python3 sample_random_bg_from_bam_stream.py --exclude union.bed --n 30000 --seed 1 > random.bed

Input:
  SAM records on stdin (no header required).

Output:
  BED 3-col (chrom, start0, end0) to stdout, sorted by (chrom,pos0).
  Summary stats to stderr.
"""

from __future__ import annotations

import argparse
import random
import sys
from typing import Set, Tuple, List


def load_exclude_bed(path: str) -> Set[Tuple[str, int]]:
    excl: Set[Tuple[str, int]] = set()
    with open(path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
            except ValueError:
                continue
            excl.add((chrom, start))
    return excl


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exclude", default="", help="Optional BED to exclude (chrom,start) loci")
    ap.add_argument("--n", type=int, default=30000, help="Number of loci to sample")
    ap.add_argument("--seed", type=int, default=1, help="Deterministic RNG seed")
    args = ap.parse_args()

    exclude: Set[Tuple[str, int]] = set()
    if args.exclude:
        exclude = load_exclude_bed(args.exclude)

    rng = random.Random(args.seed)

    seen: Set[Tuple[str, int]] = set()
    reservoir: List[Tuple[str, int]] = []
    eligible_unique = 0
    unmapped_or_bad = 0

    for line in sys.stdin:
        if not line or line[0] == "@":
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 5:
            unmapped_or_bad += 1
            continue
        chrom = parts[2]
        if chrom == "*" or chrom == "":
            unmapped_or_bad += 1
            continue
        try:
            pos1 = int(parts[3])
        except ValueError:
            unmapped_or_bad += 1
            continue
        pos0 = pos1 - 1
        if pos0 < 0:
            unmapped_or_bad += 1
            continue

        key = (chrom, pos0)
        if key in exclude:
            continue
        if key in seen:
            continue
        seen.add(key)
        eligible_unique += 1

        if len(reservoir) < args.n:
            reservoir.append(key)
        else:
            j = rng.randrange(eligible_unique)
            if j < args.n:
                reservoir[j] = key

    reservoir.sort(key=lambda x: (x[0], x[1]))
    for chrom, pos0 in reservoir:
        sys.stdout.write(f"{chrom}\t{pos0}\t{pos0+1}\n")

    sys.stderr.write(
        f"sample_random_bg_from_bam_stream: exclude={len(exclude)} eligible_unique={eligible_unique} sampled={len(reservoir)} bad={unmapped_or_bad}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

