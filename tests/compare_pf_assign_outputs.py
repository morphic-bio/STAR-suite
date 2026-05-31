#!/usr/bin/env python3
"""Canonical parity check for assignBarcodes output directories."""

import argparse
import csv
import gzip
from pathlib import Path
from typing import Dict, List, Tuple


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def resolve_file(directory: Path, basename: str) -> Path:
    plain = directory / basename
    gz = directory / f"{basename}.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise FileNotFoundError(f"missing {basename}(.gz) in {directory}")


def read_lines(directory: Path, basename: str) -> List[str]:
    path = resolve_file(directory, basename)
    with open_text(path) as handle:
        return [line.rstrip("\n") for line in handle]


def read_matrix(directory: Path) -> Tuple[List[str], List[str], Dict[Tuple[int, str], int]]:
    features = read_lines(directory, "features.txt")
    barcodes = read_lines(directory, "barcodes.txt")
    matrix_path = resolve_file(directory, "matrix.mtx")
    counts: Dict[Tuple[int, str], int] = {}
    with open_text(matrix_path) as handle:
        header = handle.readline().strip()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"{matrix_path}: invalid MatrixMarket header")
        dims_seen = False
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            if not dims_seen:
                dims_seen = True
                continue
            row = int(parts[0])
            col = int(parts[1])
            value = int(float(parts[2]))
            if row <= 0 or row > len(features):
                raise ValueError(f"{matrix_path}: row {row} outside feature list")
            if col <= 0 or col > len(barcodes):
                raise ValueError(f"{matrix_path}: column {col} outside barcode list")
            key = (row, barcodes[col - 1])
            counts[key] = counts.get(key, 0) + value
    return features, barcodes, counts


def read_csv_rows(directory: Path, basename: str) -> Tuple[Tuple[str, ...], List[Tuple[str, ...]]]:
    path = resolve_file(directory, basename)
    with open_text(path) as handle:
        reader = csv.reader(handle)
        try:
            header = tuple(next(reader))
        except StopIteration:
            return tuple(), []
        rows = [tuple(row) for row in reader]
    return header, sorted(rows)


def comparable_stats(directory: Path) -> List[str]:
    try:
        lines = read_lines(directory, "stats.txt")
    except FileNotFoundError:
        return []
    return [line for line in lines if line.strip()]


def compare_equal(label: str, left, right, errors: List[str]) -> None:
    if left != right:
        errors.append(label)


def compare_outputs(left_dir: Path, right_dir: Path) -> List[str]:
    errors: List[str] = []
    left_features, left_barcodes, left_counts = read_matrix(left_dir)
    right_features, right_barcodes, right_counts = read_matrix(right_dir)

    compare_equal("features.txt", left_features, right_features, errors)
    compare_equal("barcode set", set(left_barcodes), set(right_barcodes), errors)
    compare_equal("matrix counts by feature row and barcode", left_counts, right_counts, errors)
    compare_equal(
        "deduped_counts_histograms.txt",
        read_lines(left_dir, "deduped_counts_histograms.txt"),
        read_lines(right_dir, "deduped_counts_histograms.txt"),
        errors,
    )
    compare_equal(
        "feature_per_cell.csv sorted rows",
        read_csv_rows(left_dir, "feature_per_cell.csv"),
        read_csv_rows(right_dir, "feature_per_cell.csv"),
        errors,
    )
    compare_equal("stats.txt", comparable_stats(left_dir), comparable_stats(right_dir), errors)
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare assignBarcodes outputs independent of barcode row order."
    )
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    args = parser.parse_args()

    errors = compare_outputs(args.left, args.right)
    if errors:
        print("FAIL: assignBarcodes outputs differ")
        for error in errors:
            print(f"  - {error}")
        return 1
    print(f"PASS: assignBarcodes canonical parity {args.left} {args.right}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
