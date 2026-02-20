#!/usr/bin/env python3
"""Compare STAR exact matches (m=0) against Cell Ranger raw/calls.

Outputs:
- STAR_EXACT_VS_CR.tsv
- STAR_EXACT_VS_CR_SUMMARY.txt
- STAR_EXACT_VS_CR_MISSING_TOP.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def load_lines(path: Path) -> List[str]:
    with open_text(path) as fh:
        return [line.rstrip("\n") for line in fh if line.strip()]


def parse_matrix_market(path: Path) -> Dict[Tuple[int, int], int]:
    counts: Dict[Tuple[int, int], int] = defaultdict(int)
    with open_text(path) as fh:
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
            value = int(round(float(val_s)))
            if value:
                counts[(row, col)] += value
    return counts


def load_translation_nxt_to_tru(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open_text(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 2:
                continue
            nxt, tru = parts[0], parts[1]
            mapping[nxt] = tru
    return mapping


def strip_barcode_suffix(bc: str) -> str:
    return bc[:-2] if bc.endswith("-1") else bc


def parse_feature_call(call_value: str) -> Iterable[str]:
    call_value = (call_value or "").strip()
    if not call_value or call_value.lower() == "none":
        return []
    return [tok for tok in call_value.split("|") if tok and tok.lower() != "none"]


def load_cr_calls_set(calls_csv: Path) -> Set[Tuple[str, str]]:
    called: Set[Tuple[str, str]] = set()
    with open_text(calls_csv) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            bc = strip_barcode_suffix(row.get("cell_barcode", "").strip())
            if not bc:
                continue
            for feature in parse_feature_call(row.get("feature_call", "")):
                called.add((feature, bc))
    return called


def load_cr_raw_feature_counts(
    cr_mex_dir: Path,
    feature_type: str,
) -> Dict[Tuple[str, str], int]:
    features_path = cr_mex_dir / "features.tsv.gz"
    barcodes_path = cr_mex_dir / "barcodes.tsv.gz"
    matrix_path = cr_mex_dir / "matrix.mtx.gz"

    row_to_feature: Dict[int, str] = {}
    with gzip.open(features_path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[2] == feature_type:
                row_to_feature[i] = parts[1]

    col_to_barcode: Dict[int, str] = {}
    with gzip.open(barcodes_path, "rt") as fh:
        for i, line in enumerate(fh, start=1):
            col_to_barcode[i] = strip_barcode_suffix(line.strip())

    counts: Dict[Tuple[str, str], int] = defaultdict(int)
    with gzip.open(matrix_path, "rt") as fh:
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
            if row not in row_to_feature:
                continue
            feature = row_to_feature[row]
            barcode = col_to_barcode[int(col_s)]
            value = int(val_s)
            if value:
                counts[(feature, barcode)] += value
    return counts


def partner_feature(feature: str) -> str:
    if feature.endswith("_A"):
        return f"{feature[:-2]}_B"
    if feature.endswith("_B"):
        return f"{feature[:-2]}_A"
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare STAR exact (m=0) pairs vs Cell Ranger raw/calls."
    )
    parser.add_argument("--star-dir", required=True, help="STAR m=0 output directory")
    parser.add_argument("--cr-mex-dir", required=True, help="Cell Ranger raw_feature_bc_matrix directory")
    parser.add_argument("--cr-calls-csv", required=True, help="Cell Ranger protospacer_calls_per_cell.csv")
    parser.add_argument("--translation-file", required=True, help="2-column NXT/TRU translation whitelist")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--feature-type", default="CRISPR Guide Capture", help="Feature type in CR MEX")
    parser.add_argument("--top-missing", type=int, default=200, help="Rows to keep in STAR_EXACT_VS_CR_MISSING_TOP.tsv")
    args = parser.parse_args()

    star_dir = Path(args.star_dir)
    cr_mex_dir = Path(args.cr_mex_dir)
    cr_calls_csv = Path(args.cr_calls_csv)
    translation_file = Path(args.translation_file)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    star_features = load_lines(star_dir / "features.txt")
    star_barcodes = load_lines(star_dir / "barcodes.txt")
    star_counts = parse_matrix_market(star_dir / "matrix.mtx")

    nxt_to_tru = load_translation_nxt_to_tru(translation_file)
    cr_raw_counts = load_cr_raw_feature_counts(cr_mex_dir, args.feature_type)
    cr_calls = load_cr_calls_set(cr_calls_csv)

    rows: List[Tuple[str, str, str, int, int, int]] = []
    for (row_idx, col_idx), star_count in star_counts.items():
        feature = star_features[row_idx - 1]
        bc_nxt = star_barcodes[col_idx - 1]
        bc_tru = nxt_to_tru.get(bc_nxt, bc_nxt)
        cr_raw = cr_raw_counts.get((feature, bc_tru), 0)
        in_cr_call = 1 if (feature, bc_tru) in cr_calls else 0
        rows.append((feature, bc_nxt, bc_tru, star_count, cr_raw, in_cr_call))

    rows.sort(key=lambda r: (-r[3], r[0], r[2]))

    out_tsv = out_dir / "STAR_EXACT_VS_CR.tsv"
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "feature",
                "barcode_nxt",
                "barcode_tru",
                "star_exact_count",
                "cr_raw_count",
                "in_cr_call_list",
            ]
        )
        writer.writerows(rows)

    missing = [r for r in rows if r[4] == 0]
    top_missing = missing[: args.top_missing]

    out_top = out_dir / "STAR_EXACT_VS_CR_MISSING_TOP.tsv"
    with open(out_top, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "feature",
                "barcode_nxt",
                "barcode_tru",
                "star_exact_count",
                "cr_raw_count",
                "in_cr_call_list",
                "partner_feature",
                "cr_partner_count",
                "missing_class",
            ]
        )
        for feature, bc_nxt, bc_tru, star_count, cr_raw, in_call in top_missing:
            partner = partner_feature(feature)
            cr_partner = cr_raw_counts.get((partner, bc_tru), 0) if partner else 0
            if partner and cr_partner > 0:
                mclass = "shift_to_partner_feature"
            else:
                mclass = "no_partner_signal"
            writer.writerow(
                [
                    feature,
                    bc_nxt,
                    bc_tru,
                    star_count,
                    cr_raw,
                    in_call,
                    partner,
                    cr_partner,
                    mclass,
                ]
            )

    pairs_total = len(rows)
    star_exact_sum = sum(r[3] for r in rows)
    pairs_cr_raw_gt0 = sum(1 for r in rows if r[4] > 0)
    pairs_cr_raw_eq0 = pairs_total - pairs_cr_raw_gt0
    star_sum_cr_raw_gt0 = sum(r[3] for r in rows if r[4] > 0)
    star_sum_cr_raw_eq0 = star_exact_sum - star_sum_cr_raw_gt0
    pairs_cr_call_yes = sum(1 for r in rows if r[5] == 1)
    pairs_cr_call_no = pairs_total - pairs_cr_call_yes

    missing_shift_pairs = 0
    missing_shift_sum = 0
    missing_no_partner_pairs = 0
    missing_no_partner_sum = 0
    for feature, _bc_nxt, bc_tru, star_count, cr_raw, _in_call in missing:
        if cr_raw != 0:
            continue
        partner = partner_feature(feature)
        if partner and cr_raw_counts.get((partner, bc_tru), 0) > 0:
            missing_shift_pairs += 1
            missing_shift_sum += star_count
        else:
            missing_no_partner_pairs += 1
            missing_no_partner_sum += star_count

    out_summary = out_dir / "STAR_EXACT_VS_CR_SUMMARY.txt"
    with open(out_summary, "w") as fh:
        fh.write(f"generated_utc={datetime.now(timezone.utc).isoformat()}\n")
        fh.write(f"star_dir={star_dir}\n")
        fh.write(f"cr_mex_dir={cr_mex_dir}\n")
        fh.write(f"cr_calls_csv={cr_calls_csv}\n")
        fh.write(f"translation_file={translation_file}\n")
        fh.write(f"feature_type={args.feature_type}\n")
        fh.write(f"pairs_total={pairs_total}\n")
        fh.write(f"star_exact_sum={star_exact_sum}\n")
        fh.write(f"pairs_cr_raw_gt0={pairs_cr_raw_gt0}\n")
        fh.write(f"pairs_cr_raw_eq0={pairs_cr_raw_eq0}\n")
        fh.write(f"pairs_cr_raw_gt0_pct={100.0 * pairs_cr_raw_gt0 / pairs_total:.4f}\n")
        fh.write(f"pairs_cr_raw_eq0_pct={100.0 * pairs_cr_raw_eq0 / pairs_total:.4f}\n")
        fh.write(f"star_sum_cr_raw_gt0={star_sum_cr_raw_gt0}\n")
        fh.write(f"star_sum_cr_raw_eq0={star_sum_cr_raw_eq0}\n")
        fh.write(f"star_sum_cr_raw_eq0_pct={100.0 * star_sum_cr_raw_eq0 / star_exact_sum:.4f}\n")
        fh.write(f"pairs_cr_call_yes={pairs_cr_call_yes}\n")
        fh.write(f"pairs_cr_call_no={pairs_cr_call_no}\n")
        fh.write(f"missing_shift_to_partner_pairs={missing_shift_pairs}\n")
        fh.write(f"missing_shift_to_partner_star_sum={missing_shift_sum}\n")
        fh.write(f"missing_no_partner_signal_pairs={missing_no_partner_pairs}\n")
        fh.write(f"missing_no_partner_signal_star_sum={missing_no_partner_sum}\n")

    print(out_tsv)
    print(out_summary)
    print(out_top)


if __name__ == "__main__":
    main()
