#!/usr/bin/env python3
"""Classify STAR m=1 rescue pairs that are missing in Cell Ranger raw matrix.

Input: STAR_M1_DELTA_VS_CR.tsv
Output:
- STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.tsv
- STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.summary.txt
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
        return [line.rstrip("\n") for line in fh if line.strip()]


def strip_barcode_suffix(bc: str) -> str:
    return bc[:-2] if bc.endswith("-1") else bc


def partner_feature(feature: str) -> str:
    if feature.endswith("_A"):
        return f"{feature[:-2]}_B"
    if feature.endswith("_B"):
        return f"{feature[:-2]}_A"
    return ""


def parse_needed_pairs(delta_table: Path):
    rows = []
    needed_cr: Set[Tuple[str, str]] = set()
    needed_star: Set[Tuple[str, str]] = set()
    needed_cr_barcodes: Set[str] = set()

    with open(delta_table, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {
            "feature",
            "barcode_nxt",
            "barcode_tru",
            "star_m0_count",
            "star_m1_count",
            "star_m1_minus_m0",
            "cr_raw_count",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns in {delta_table}: {sorted(missing)}")

        for row in reader:
            if int(row["cr_raw_count"]) != 0:
                continue

            feature = row["feature"]
            bc_nxt = row["barcode_nxt"]
            bc_tru = row["barcode_tru"]
            partner = partner_feature(feature)

            rec = {
                "feature": feature,
                "barcode_nxt": bc_nxt,
                "barcode_tru": bc_tru,
                "star_m0": int(row["star_m0_count"]),
                "star_m1": int(row["star_m1_count"]),
                "delta": int(row["star_m1_minus_m0"]),
                "partner_feature": partner,
            }
            rows.append(rec)

            needed_cr.add((feature, bc_tru))
            needed_star.add((feature, bc_nxt))
            needed_cr_barcodes.add(bc_tru)
            if partner:
                needed_cr.add((partner, bc_tru))
                needed_star.add((partner, bc_nxt))

    return rows, needed_cr, needed_star, needed_cr_barcodes


def load_cr_feature_row_map(cr_mex_dir: Path, feature_type: str) -> Dict[str, int]:
    features_path = cr_mex_dir / "features.tsv.gz"
    feature_to_row: Dict[str, int] = {}
    with gzip.open(features_path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[2] == feature_type:
                feature_to_row[parts[1]] = i
    return feature_to_row


def load_cr_barcode_col_map(cr_mex_dir: Path) -> Dict[str, int]:
    barcodes_path = cr_mex_dir / "barcodes.tsv.gz"
    barcode_to_col: Dict[str, int] = {}
    with gzip.open(barcodes_path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            barcode_to_col[strip_barcode_suffix(line.strip())] = i
    return barcode_to_col


def load_star_feature_row_map(star_dir: Path) -> Dict[str, int]:
    features = load_lines(star_dir / "features.txt")
    return {name: idx for idx, name in enumerate(features, start=1)}


def load_star_barcode_col_map(star_dir: Path) -> Dict[str, int]:
    barcodes = load_lines(star_dir / "barcodes.txt")
    return {name: idx for idx, name in enumerate(barcodes, start=1)}


def parse_subset_matrix(
    matrix_path: Path,
    needed_coords: Dict[Tuple[int, int], Tuple[str, str]],
) -> Dict[Tuple[str, str], int]:
    counts: Dict[Tuple[str, str], int] = defaultdict(int)
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
            key = (int(row_s), int(col_s))
            pair = needed_coords.get(key)
            if pair is None:
                continue
            value = int(round(float(val_s)))
            if value:
                counts[pair] += value
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify CR-missing STAR m=1 rescue pairs.")
    parser.add_argument("--delta-table", required=True, help="STAR_M1_DELTA_VS_CR.tsv")
    parser.add_argument("--star-m1-dir", required=True, help="STAR m=1 output directory")
    parser.add_argument("--cr-mex-dir", required=True, help="Cell Ranger raw_feature_bc_matrix directory")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--feature-type", default="CRISPR Guide Capture", help="Feature type in CR MEX")
    args = parser.parse_args()

    delta_table = Path(args.delta_table)
    star_m1_dir = Path(args.star_m1_dir)
    cr_mex_dir = Path(args.cr_mex_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows, needed_cr, needed_star, needed_cr_barcodes = parse_needed_pairs(delta_table)

    cr_feature_to_row = load_cr_feature_row_map(cr_mex_dir, args.feature_type)
    cr_barcode_to_col = load_cr_barcode_col_map(cr_mex_dir)

    star_feature_to_row = load_star_feature_row_map(star_m1_dir)
    star_barcode_to_col = load_star_barcode_col_map(star_m1_dir)

    needed_cr_coords: Dict[Tuple[int, int], Tuple[str, str]] = {}
    for feature, barcode in needed_cr:
        row_idx = cr_feature_to_row.get(feature)
        col_idx = cr_barcode_to_col.get(barcode)
        if row_idx is None or col_idx is None:
            continue
        needed_cr_coords[(row_idx, col_idx)] = (feature, barcode)

    needed_star_coords: Dict[Tuple[int, int], Tuple[str, str]] = {}
    for feature, barcode in needed_star:
        row_idx = star_feature_to_row.get(feature)
        col_idx = star_barcode_to_col.get(barcode)
        if row_idx is None or col_idx is None:
            continue
        needed_star_coords[(row_idx, col_idx)] = (feature, barcode)

    cr_counts = parse_subset_matrix(cr_mex_dir / "matrix.mtx.gz", needed_cr_coords)
    star_counts = parse_subset_matrix(star_m1_dir / "matrix.mtx", needed_star_coords)

    out_rows = []
    class_pairs: Dict[str, int] = defaultdict(int)
    class_delta: Dict[str, int] = defaultdict(int)

    for rec in sorted(rows, key=lambda r: (-r["delta"], r["feature"], r["barcode_tru"])):
        feature = rec["feature"]
        barcode_tru = rec["barcode_tru"]
        barcode_nxt = rec["barcode_nxt"]
        partner = rec["partner_feature"]

        cr_self = cr_counts.get((feature, barcode_tru), 0)
        cr_partner = cr_counts.get((partner, barcode_tru), 0) if partner else 0
        star_partner = star_counts.get((partner, barcode_nxt), 0) if partner else 0

        barcode_in_cr_raw = 1 if barcode_tru in cr_barcode_to_col else 0

        if not barcode_in_cr_raw:
            classification = "barcode_absent_from_cr_matrix"
        elif cr_partner > 0:
            classification = "shift_to_partner_feature"
        else:
            classification = "no_partner_signal_in_cr"

        class_pairs[classification] += 1
        class_delta[classification] += rec["delta"]

        out_rows.append(
            [
                feature,
                barcode_nxt,
                barcode_tru,
                rec["delta"],
                rec["star_m0"],
                rec["star_m1"],
                cr_self,
                partner,
                star_partner,
                cr_partner,
                barcode_in_cr_raw,
                classification,
            ]
        )

    out_tsv = out_dir / "STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.tsv"
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "feature",
                "barcode_nxt",
                "barcode_tru",
                "delta",
                "star_m0",
                "star_m1",
                "cr_self",
                "partner_feature",
                "star_partner",
                "cr_partner",
                "barcode_in_cr_raw",
                "classification",
            ]
        )
        writer.writerows(out_rows)

    missing_pairs = len(out_rows)
    missing_delta_sum = sum(int(r[3]) for r in out_rows)

    out_summary = out_dir / "STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.summary.txt"
    with open(out_summary, "w") as fh:
        fh.write(f"generated_utc={datetime.now(timezone.utc).isoformat()}\n")
        fh.write(f"delta_table={delta_table}\n")
        fh.write(f"star_m1_dir={star_m1_dir}\n")
        fh.write(f"cr_mex_dir={cr_mex_dir}\n")
        fh.write(f"feature_type={args.feature_type}\n")
        fh.write(f"missing_pairs={missing_pairs}\n")
        fh.write(f"missing_delta_sum={missing_delta_sum}\n")
        for key in sorted(class_pairs.keys()):
            fh.write(f"{key}_pairs={class_pairs[key]}\n")
            fh.write(f"{key}_delta_sum={class_delta[key]}\n")

    print(out_tsv)
    print(out_summary)


if __name__ == "__main__":
    main()
