#!/usr/bin/env python3
"""
Evaluate a gray-zone admission policy using a known EM-like reference label set.

This harness is intentionally post hoc:
- the base run defines the raw matrix, filtered cutoff, and baseline labels
- the EM-like run provides "known" admission labels in the gray zone
- CR filtered barcodes are the external target for overlap metrics

It does NOT change counts. It only asks whether swapping admission labels inside
the gray zone improves filtered-cell overlap.
"""

from __future__ import annotations

import argparse
import gzip
from array import array
from pathlib import Path
from typing import Iterable, Set, Tuple


def open_maybe_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def normalize_barcode(raw: str) -> str:
    bc = raw.strip()
    if not bc:
        return ""
    if "-" in bc:
        base, suffix = bc.rsplit("-", 1)
        if suffix.isdigit():
            return base
    return bc


def read_barcodes(path: Path) -> Set[str]:
    out = set()
    with open_maybe_gz(path) as handle:
        for raw in handle:
            bc = normalize_barcode(raw)
            if bc:
                out.add(bc)
    return out


def read_sorted_umi_counts(path: Path) -> list[int]:
    vals = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            s = raw.strip()
            if s:
                vals.append(int(s))
    return vals


def count_lines(path: Path) -> int:
    n = 0
    with open(path, "r", encoding="utf-8") as handle:
        for _ in handle:
            n += 1
    return n


def stream_barcode_totals(raw_dir: Path) -> Tuple[array, int]:
    barcodes_path = raw_dir / "barcodes.tsv"
    matrix_path = raw_dir / "matrix.mtx"
    if not barcodes_path.exists():
        raise SystemExit(f"Missing raw barcodes.tsv: {barcodes_path}")
    if not matrix_path.exists():
        raise SystemExit(f"Missing raw matrix.mtx: {matrix_path}")

    n_barcodes = count_lines(barcodes_path)
    totals = array("Q", [0]) * n_barcodes

    with open(matrix_path, "r", encoding="utf-8") as handle:
        dims_seen = False
        for raw in handle:
            if raw.startswith("%"):
                continue
            if not dims_seen:
                dims = raw.split()
                if len(dims) != 3:
                    raise SystemExit(f"Malformed matrix dimensions line in {matrix_path}: {raw.strip()}")
                n_cols = int(dims[1])
                if n_cols != n_barcodes:
                    raise SystemExit(
                        f"Barcode count mismatch: {n_barcodes} lines in {barcodes_path}, "
                        f"matrix reports {n_cols} columns"
                    )
                dims_seen = True
                continue

            parts = raw.split()
            if len(parts) != 3:
                continue
            col = int(parts[1]) - 1
            val = int(float(parts[2]))
            totals[col] += val

    return totals, n_barcodes


def build_gray_zone_set(raw_dir: Path, totals: array, low_umi: int, high_umi: int) -> Set[str]:
    barcodes_path = raw_dir / "barcodes.tsv"
    out = set()
    with open(barcodes_path, "r", encoding="utf-8") as handle:
        for idx, raw in enumerate(handle):
            total = totals[idx]
            if low_umi <= total <= high_umi:
                bc = normalize_barcode(raw)
                if bc:
                    out.add(bc)
    return out


def jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 1.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 1.0


def summarize_overlap(label: str, sample: Set[str], cr: Set[str]) -> dict[str, str]:
    inter = len(sample & cr)
    union = len(sample | cr)
    return {
        "label": label,
        "cells": str(len(sample)),
        "common_cr": str(inter),
        "union_cr": str(union),
        "jaccard_cr": f"{(inter / union) if union else 1.0:.6f}",
    }


def write_barcode_list(path: Path, barcodes: Iterable[str]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for bc in sorted(barcodes):
            handle.write(f"{bc}\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-raw-dir", type=Path, required=True, help="Base run raw GeneFull MEX dir")
    ap.add_argument("--base-filtered-barcodes", type=Path, required=True, help="Base run filtered barcodes.tsv")
    ap.add_argument("--base-umi-sorted", type=Path, required=True, help="Base run UMIperCellSorted.txt")
    ap.add_argument("--em-filtered-barcodes", type=Path, required=True, help="EM-like reference filtered barcodes.tsv")
    ap.add_argument("--cr-filtered-barcodes", type=Path, required=True, help="CR filtered barcodes.tsv(.gz)")
    ap.add_argument("--window-before", type=int, default=2000, help="Ranks above the base filtered cutoff to include in the gray zone")
    ap.add_argument("--window-after", type=int, default=2000, help="Ranks below the base filtered cutoff to include in the gray zone")
    ap.add_argument("--outdir", type=Path, required=True, help="Output directory")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    base_filtered = read_barcodes(args.base_filtered_barcodes)
    em_filtered = read_barcodes(args.em_filtered_barcodes)
    cr_filtered = read_barcodes(args.cr_filtered_barcodes)
    umi_sorted = read_sorted_umi_counts(args.base_umi_sorted)

    base_n = len(base_filtered)
    if base_n == 0:
        raise SystemExit("Base filtered barcode set is empty")
    if len(umi_sorted) < base_n:
        raise SystemExit(
            f"UMIperCellSorted has only {len(umi_sorted)} entries, fewer than base filtered cells {base_n}"
        )

    high_rank = max(1, base_n - args.window_before)
    low_rank = min(len(umi_sorted), base_n + args.window_after)
    high_umi = umi_sorted[high_rank - 1]
    low_umi = umi_sorted[low_rank - 1]

    totals, n_barcodes = stream_barcode_totals(args.base_raw_dir)
    gray_zone = build_gray_zone_set(args.base_raw_dir, totals, low_umi, high_umi)

    base_outside = base_filtered - gray_zone
    hybrid_filtered = base_outside | (em_filtered & gray_zone)

    base_summary = summarize_overlap("base", base_filtered, cr_filtered)
    em_summary = summarize_overlap("em_reference", em_filtered, cr_filtered)
    hybrid_summary = summarize_overlap("hybrid_gray_em", hybrid_filtered, cr_filtered)

    base_gray = base_filtered & gray_zone
    em_gray = em_filtered & gray_zone
    hybrid_gray = hybrid_filtered & gray_zone
    added = hybrid_filtered - base_filtered
    removed = base_filtered - hybrid_filtered

    summary_path = args.outdir / "summary.tsv"
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write("metric\tvalue\n")
        out.write(f"base_filtered_cells\t{len(base_filtered)}\n")
        out.write(f"em_filtered_cells\t{len(em_filtered)}\n")
        out.write(f"hybrid_filtered_cells\t{len(hybrid_filtered)}\n")
        out.write(f"cr_filtered_cells\t{len(cr_filtered)}\n")
        out.write(f"raw_barcodes\t{n_barcodes}\n")
        out.write(f"gray_zone_window_before\t{args.window_before}\n")
        out.write(f"gray_zone_window_after\t{args.window_after}\n")
        out.write(f"gray_zone_rank_high\t{high_rank}\n")
        out.write(f"gray_zone_rank_low\t{low_rank}\n")
        out.write(f"gray_zone_umi_high\t{high_umi}\n")
        out.write(f"gray_zone_umi_low\t{low_umi}\n")
        out.write(f"gray_zone_barcodes\t{len(gray_zone)}\n")
        out.write(f"base_gray_barcodes\t{len(base_gray)}\n")
        out.write(f"em_gray_barcodes\t{len(em_gray)}\n")
        out.write(f"hybrid_gray_barcodes\t{len(hybrid_gray)}\n")
        out.write(f"hybrid_added_vs_base\t{len(added)}\n")
        out.write(f"hybrid_removed_vs_base\t{len(removed)}\n")
        out.write(f"hybrid_added_in_cr\t{len(added & cr_filtered)}\n")
        out.write(f"hybrid_added_not_in_cr\t{len(added - cr_filtered)}\n")
        out.write(f"hybrid_removed_in_cr\t{len(removed & cr_filtered)}\n")
        out.write(f"hybrid_removed_not_in_cr\t{len(removed - cr_filtered)}\n")
        for row in (base_summary, em_summary, hybrid_summary):
            for key, value in row.items():
                if key == "label":
                    continue
                out.write(f"{row['label']}_{key}\t{value}\n")

    write_barcode_list(args.outdir / "gray_zone_barcodes.tsv", gray_zone)
    write_barcode_list(args.outdir / "hybrid_filtered_barcodes.tsv", hybrid_filtered)
    write_barcode_list(args.outdir / "hybrid_added_vs_base.tsv", added)
    write_barcode_list(args.outdir / "hybrid_removed_vs_base.tsv", removed)
    write_barcode_list(args.outdir / "hybrid_added_in_cr.tsv", added & cr_filtered)
    write_barcode_list(args.outdir / "hybrid_removed_in_cr.tsv", removed & cr_filtered)

    print(summary_path.read_text(), end="")
    print(f"Wrote gray-zone admission evaluation to: {args.outdir}")


if __name__ == "__main__":
    main()
