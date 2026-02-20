#!/usr/bin/env python3
"""Inspect per-barcode CRISPR feature totals in STAR and Cell Ranger raw MEX.

Useful for drilling into specific discrepancy examples (e.g. FOXD3 barcode).
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Set, Tuple


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def load_lines(path: Path) -> List[str]:
    with open_text(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def strip_barcode_suffix(bc: str) -> str:
    return bc[:-2] if bc.endswith("-1") else bc


def load_translation(path: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    nxt_to_tru: Dict[str, str] = {}
    tru_to_nxt: Dict[str, str] = {}
    with open_text(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            nxt, tru = parts[0], parts[1]
            nxt_to_tru[nxt] = tru
            tru_to_nxt[tru] = nxt
    return nxt_to_tru, tru_to_nxt


def parse_matrix_subset(
    matrix_path: Path,
    needed_rows: Set[int],
    needed_cols: Set[int],
) -> Dict[Tuple[int, int], int]:
    counts: Dict[Tuple[int, int], int] = defaultdict(int)
    with open_text(matrix_path) as fh:
        dims_seen = False
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            if not dims_seen:
                dims_seen = True
                continue

            row_s, col_s, val_s = line.split()
            row = int(row_s)
            col = int(col_s)
            if row not in needed_rows or col not in needed_cols:
                continue
            value = int(round(float(val_s)))
            if value:
                counts[(row, col)] += value
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inspect selected TRU barcodes across STAR and CR CRISPR features."
    )
    parser.add_argument("--barcodes-tru", default="", help="Comma-separated TRU barcodes")
    parser.add_argument("--barcodes-file", default="", help="Optional file (one TRU barcode per line)")
    parser.add_argument("--star-m1-dir", required=True, help="STAR m=1 output dir")
    parser.add_argument("--cr-mex-dir", required=True, help="CR raw_feature_bc_matrix dir")
    parser.add_argument("--translation-file", required=True, help="2-column NXT/TRU translation whitelist")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--feature-type", default="CRISPR Guide Capture", help="CR feature type")
    parser.add_argument("--top-n", type=int, default=12, help="Top N features in summary per source/barcode")
    args = parser.parse_args()

    barcodes: Set[str] = set()
    if args.barcodes_tru.strip():
        barcodes.update(b.strip() for b in args.barcodes_tru.split(",") if b.strip())
    if args.barcodes_file:
        barcodes.update(load_lines(Path(args.barcodes_file)))
    if not barcodes:
        raise ValueError("No barcodes provided. Use --barcodes-tru and/or --barcodes-file.")

    star_m1_dir = Path(args.star_m1_dir)
    cr_mex_dir = Path(args.cr_mex_dir)
    translation_file = Path(args.translation_file)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    _, tru_to_nxt = load_translation(translation_file)

    # STAR indices
    star_features = load_lines(star_m1_dir / "features.txt")
    star_barcodes = load_lines(star_m1_dir / "barcodes.txt")
    star_row_to_feature = {i: f for i, f in enumerate(star_features, start=1)}
    star_barcode_to_col = {b: i for i, b in enumerate(star_barcodes, start=1)}

    # CR indices (CRISPR only)
    cr_row_to_feature: Dict[int, str] = {}
    with gzip.open(cr_mex_dir / "features.tsv.gz", "rt") as fh:
        for i, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[2] == args.feature_type:
                cr_row_to_feature[i] = parts[1]

    cr_barcode_to_col: Dict[str, int] = {}
    with gzip.open(cr_mex_dir / "barcodes.tsv.gz", "rt") as fh:
        for i, line in enumerate(fh, start=1):
            cr_barcode_to_col[strip_barcode_suffix(line.strip())] = i

    # Resolve requested barcodes to STAR NXT barcodes
    barcode_pairs: List[Tuple[str, str]] = []
    for bc_tru in sorted(barcodes):
        bc_nxt = tru_to_nxt.get(bc_tru, bc_tru)
        barcode_pairs.append((bc_tru, bc_nxt))

    wanted_star_cols = {
        star_barcode_to_col[bc_nxt]
        for _bc_tru, bc_nxt in barcode_pairs
        if bc_nxt in star_barcode_to_col
    }
    wanted_cr_cols = {
        cr_barcode_to_col[bc_tru]
        for bc_tru, _bc_nxt in barcode_pairs
        if bc_tru in cr_barcode_to_col
    }

    star_counts = parse_matrix_subset(
        star_m1_dir / "matrix.mtx",
        set(star_row_to_feature.keys()),
        wanted_star_cols,
    )
    cr_counts = parse_matrix_subset(
        cr_mex_dir / "matrix.mtx.gz",
        set(cr_row_to_feature.keys()),
        wanted_cr_cols,
    )

    # Build feature-level records
    records: List[Tuple[str, str, str, int]] = []
    # (barcode_tru, source, feature, count)

    for bc_tru, bc_nxt in barcode_pairs:
        if bc_nxt in star_barcode_to_col:
            cidx = star_barcode_to_col[bc_nxt]
            for ridx, feature in star_row_to_feature.items():
                val = star_counts.get((ridx, cidx), 0)
                if val:
                    records.append((bc_tru, "STAR", feature, val))

        if bc_tru in cr_barcode_to_col:
            cidx = cr_barcode_to_col[bc_tru]
            for ridx, feature in cr_row_to_feature.items():
                val = cr_counts.get((ridx, cidx), 0)
                if val:
                    records.append((bc_tru, "CR", feature, val))

    records.sort(key=lambda x: (x[0], x[1], -x[3], x[2]))

    out_tsv = out_dir / "BARCODE_INSPECTION.tsv"
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["barcode_tru", "source", "feature", "count"])
        writer.writerows(records)

    # Human-readable summary
    out_summary = out_dir / "BARCODE_INSPECTION_SUMMARY.txt"
    by_key: Dict[Tuple[str, str], List[Tuple[str, int]]] = defaultdict(list)
    for bc_tru, source, feature, count in records:
        by_key[(bc_tru, source)].append((feature, count))

    with open(out_summary, "w") as fh:
        fh.write(f"generated_utc={datetime.now(timezone.utc).isoformat()}\n")
        fh.write(f"star_m1_dir={star_m1_dir}\n")
        fh.write(f"cr_mex_dir={cr_mex_dir}\n")
        fh.write(f"translation_file={translation_file}\n")
        fh.write(f"feature_type={args.feature_type}\n")
        fh.write(f"top_n={args.top_n}\n")
        for bc_tru, bc_nxt in barcode_pairs:
            fh.write("\n")
            fh.write(f"barcode_tru={bc_tru}\n")
            fh.write(f"barcode_nxt={bc_nxt}\n")
            for source in ("STAR", "CR"):
                vals = sorted(by_key.get((bc_tru, source), []), key=lambda x: (-x[1], x[0]))
                total = sum(v for _, v in vals)
                fh.write(f"[{source}] nonzero_features={len(vals)} total_counts={total}\n")
                for feature, count in vals[: args.top_n]:
                    fh.write(f"  {feature}\t{count}\n")

    print(out_tsv)
    print(out_summary)


if __name__ == "__main__":
    main()
