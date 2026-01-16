#!/usr/bin/env python3
"""
Convert GEDI *.snpdata (Location is 0-based) into a deduplicated BED goal set.

GEDI snpdata columns (tab-separated, header starts with 'Location'):
  Location   Coverage   Mismatches   Pval

Filtering for "GEDI defaults" goal set:
  - minCov (default 6)
  - minRatio (default 0.3) where ratio = mismatches / coverage
  - maxPval (default 0.001) using GEDI's own pval column

Output BED:
  chrom  start  end
with start = pos0, end = pos0+1 (1bp interval), sorted and deduplicated.

Chrom naming:
  - Keep 'chr' prefix if present
  - Otherwise, prefix 'chr' only for {1-22, X, Y, M, MT}
  - Leave contigs like KI270734.1 untouched (do NOT prepend chr)
"""

from __future__ import annotations

import argparse
import re
import sys
from typing import Iterable, Tuple, Set, List


def normalize_chrom(chrom: str) -> str:
    if chrom.startswith("chr"):
        return chrom
    # Our STAR reference uses chrM (not chrMT). Normalize both M and MT to chrM.
    if chrom in {"MT", "M"}:
        return "chrM"
    if chrom in {"X", "Y"}:
        return "chr" + chrom
    if re.fullmatch(r"[0-9]+", chrom):
        return "chr" + chrom
    # Non-canonical contig, do not add chr
    return chrom


def iter_snpdata_lines(path: str) -> Iterable[Tuple[str, float, float, float]]:
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("Location"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            loc = parts[0]
            try:
                cov = float(parts[1])
                mm = float(parts[2])
                pval = float(parts[3])
            except ValueError:
                continue
            yield loc, cov, mm, pval


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--snpdata", required=True, help="Input GEDI .snpdata")
    ap.add_argument("--out", required=True, help="Output BED (3-column, sorted, uniq)")
    ap.add_argument("--minCov", type=float, default=6.0)
    ap.add_argument("--minRatio", type=float, default=0.3)
    ap.add_argument("--maxPval", type=float, default=0.001)
    args = ap.parse_args()

    loci: Set[Tuple[str, int]] = set()
    total = 0
    kept = 0

    for loc, cov, mm, pval in iter_snpdata_lines(args.snpdata):
        total += 1
        if cov < args.minCov:
            continue
        if cov <= 0:
            continue
        ratio = mm / cov
        if ratio < args.minRatio:
            continue
        if pval >= args.maxPval:
            continue

        if ":" not in loc:
            continue
        chrom_raw, pos_raw = loc.split(":", 1)
        try:
            pos0 = int(pos_raw)
        except ValueError:
            continue

        chrom = normalize_chrom(chrom_raw)
        loci.add((chrom, pos0))
        kept += 1

    # Sort
    def chrom_key(c: str) -> Tuple[int, str]:
        # Put chr1..chr22, chrX, chrY, chrM/chrMT first; others lexicographically after.
        if c.startswith("chr"):
            core = c[3:]
            if core.isdigit():
                return (0, f"{int(core):03d}")
            if core in {"X", "Y"}:
                return (0, {"X": "023", "Y": "024"}[core])
            if core in {"M", "MT"}:
                return (0, "025")
        return (1, c)

    sorted_loci: List[Tuple[str, int]] = sorted(loci, key=lambda x: (chrom_key(x[0]), x[1]))

    with open(args.out, "w") as out:
        for chrom, pos0 in sorted_loci:
            out.write(f"{chrom}\t{pos0}\t{pos0+1}\n")

    sys.stderr.write(
        f"convert_gedi_snpdata_to_goal_bed: total_lines={total} kept_after_filters={kept} unique_loci={len(sorted_loci)} out={args.out}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

