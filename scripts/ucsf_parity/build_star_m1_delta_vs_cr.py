#!/usr/bin/env python3
"""Build STAR m=1 rescue delta table and compare against Cell Ranger raw/call outputs.

This reproduces the analysis table used for UCSF parity debugging:
- STAR_M1_DELTA_VS_CR.tsv
- STAR_M1_DELTA_VS_CR_SUMMARY.txt
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


def parse_matrix_market(path: Path) -> Dict[Tuple[int, int], int]:
    counts: Dict[Tuple[int, int], int] = defaultdict(int)
    with open_text(path) as fh:
        dims_seen = False
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            if not dims_seen:
                parts = line.split()
                if len(parts) != 3:
                    raise ValueError(f"Invalid MatrixMarket dimensions line in {path}: {line}")
                dims_seen = True
                continue

            parts = line.split()
            if len(parts) != 3:
                raise ValueError(f"Invalid MatrixMarket entry line in {path}: {line}")
            row = int(parts[0])
            col = int(parts[1])
            value = int(round(float(parts[2])))
            if value:
                counts[(row, col)] += value
    return counts


def to_named_counts(
    counts: Dict[Tuple[int, int], int],
    features: List[str],
    barcodes: List[str],
) -> Dict[Tuple[str, str], int]:
    named: Dict[Tuple[str, str], int] = defaultdict(int)
    for (row_idx, col_idx), value in counts.items():
        feature = features[row_idx - 1]
        barcode = barcodes[col_idx - 1]
        named[(feature, barcode)] += value
    return named


def load_lines(path: Path) -> List[str]:
    with open_text(path) as fh:
        return [line.rstrip("\n") for line in fh if line.strip()]


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


def load_cr_raw_feature_counts(
    cr_mex_dir: Path,
    feature_type: str,
) -> Dict[Tuple[str, str], int]:
    features_path = cr_mex_dir / "features.tsv.gz"
    barcodes_path = cr_mex_dir / "barcodes.tsv.gz"
    matrix_path = cr_mex_dir / "matrix.mtx.gz"

    if not features_path.exists() or not barcodes_path.exists() or not matrix_path.exists():
        raise FileNotFoundError(
            "CR raw MEX is missing one or more files: features.tsv.gz, barcodes.tsv.gz, matrix.mtx.gz"
        )

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


def parse_feature_call(call_value: str) -> Iterable[str]:
    call_value = (call_value or "").strip()
    if not call_value or call_value.lower() == "none":
        return []
    return [tok for tok in call_value.split("|") if tok and tok.lower() != "none"]


def load_cr_calls_set(calls_csv: Path) -> Set[Tuple[str, str]]:
    if not calls_csv.exists():
        raise FileNotFoundError(f"CR calls CSV not found: {calls_csv}")

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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build STAR m=1 rescue delta table vs Cell Ranger raw/calls."
    )
    parser.add_argument("--star-m1-dir", required=True, help="STAR m=1 output directory (features.txt/barcodes.txt/matrix.mtx)")
    parser.add_argument("--star-m0-dir", required=True, help="STAR m=0 output directory (features.txt/barcodes.txt/matrix.mtx)")
    parser.add_argument("--cr-mex-dir", required=True, help="Cell Ranger raw_feature_bc_matrix directory")
    parser.add_argument("--cr-calls-csv", required=True, help="Cell Ranger protospacer_calls_per_cell.csv")
    parser.add_argument("--translation-file", required=True, help="2-column NXT/TRU barcode translation file (.txt or .txt.gz)")
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument("--feature-type", default="CRISPR Guide Capture", help="Feature type to read from CR MEX")
    parser.add_argument("--min-delta", type=int, default=1, help="Minimum STAR (m1-m0) delta to include (default: 1)")
    args = parser.parse_args()

    star_m1_dir = Path(args.star_m1_dir)
    star_m0_dir = Path(args.star_m0_dir)
    cr_mex_dir = Path(args.cr_mex_dir)
    cr_calls_csv = Path(args.cr_calls_csv)
    translation_file = Path(args.translation_file)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    star_features_m1 = load_lines(star_m1_dir / "features.txt")
    star_barcodes_m1 = load_lines(star_m1_dir / "barcodes.txt")
    star_features_m0 = load_lines(star_m0_dir / "features.txt")
    star_barcodes_m0 = load_lines(star_m0_dir / "barcodes.txt")

    star_m1_counts = parse_matrix_market(star_m1_dir / "matrix.mtx")
    star_m0_counts = parse_matrix_market(star_m0_dir / "matrix.mtx")
    star_m1_named = to_named_counts(star_m1_counts, star_features_m1, star_barcodes_m1)
    star_m0_named = to_named_counts(star_m0_counts, star_features_m0, star_barcodes_m0)

    nxt_to_tru = load_translation_nxt_to_tru(translation_file)
    cr_raw_counts = load_cr_raw_feature_counts(cr_mex_dir, args.feature_type)
    cr_called_pairs = load_cr_calls_set(cr_calls_csv)

    rows: List[Tuple[str, str, str, int, int, int, int, int]] = []
    for (feature, bc_nxt), m1_count in star_m1_named.items():
        m0_count = star_m0_named.get((feature, bc_nxt), 0)
        delta = m1_count - m0_count
        if delta < args.min_delta:
            continue

        bc_tru = nxt_to_tru.get(bc_nxt, bc_nxt)

        cr_raw = cr_raw_counts.get((feature, bc_tru), 0)
        in_cr_call = 1 if (feature, bc_tru) in cr_called_pairs else 0

        rows.append((feature, bc_nxt, bc_tru, m0_count, m1_count, delta, cr_raw, in_cr_call))

    rows.sort(key=lambda r: (-r[5], r[0], r[2]))

    out_tsv = out_dir / "STAR_M1_DELTA_VS_CR.tsv"
    with open(out_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "feature",
                "barcode_nxt",
                "barcode_tru",
                "star_m0_count",
                "star_m1_count",
                "star_m1_minus_m0",
                "cr_raw_count",
                "in_cr_call_list",
            ]
        )
        writer.writerows(rows)

    pairs_total = len(rows)
    delta_sum = sum(r[5] for r in rows)
    pairs_cr_raw_gt0 = sum(1 for r in rows if r[6] > 0)
    pairs_cr_raw_eq0 = pairs_total - pairs_cr_raw_gt0
    delta_sum_cr_raw_gt0 = sum(r[5] for r in rows if r[6] > 0)
    delta_sum_cr_raw_eq0 = delta_sum - delta_sum_cr_raw_gt0
    pairs_cr_call_yes = sum(1 for r in rows if r[7] == 1)
    pairs_cr_call_no = pairs_total - pairs_cr_call_yes
    delta_sum_cr_call_yes = sum(r[5] for r in rows if r[7] == 1)
    delta_sum_cr_call_no = delta_sum - delta_sum_cr_call_yes

    out_summary = out_dir / "STAR_M1_DELTA_VS_CR_SUMMARY.txt"
    with open(out_summary, "w") as fh:
        fh.write(f"generated_utc={datetime.now(timezone.utc).isoformat()}\n")
        fh.write(f"star_m1_dir={star_m1_dir}\n")
        fh.write(f"star_m0_dir={star_m0_dir}\n")
        fh.write(f"cr_mex_dir={cr_mex_dir}\n")
        fh.write(f"cr_calls_csv={cr_calls_csv}\n")
        fh.write(f"translation_file={translation_file}\n")
        fh.write(f"feature_type={args.feature_type}\n")
        fh.write(f"min_delta={args.min_delta}\n")
        fh.write(f"pairs_total={pairs_total}\n")
        fh.write(f"delta_sum={delta_sum}\n")
        fh.write(f"pairs_cr_raw_gt0={pairs_cr_raw_gt0}\n")
        fh.write(f"pairs_cr_raw_eq0={pairs_cr_raw_eq0}\n")
        fh.write(f"delta_sum_cr_raw_gt0={delta_sum_cr_raw_gt0}\n")
        fh.write(f"delta_sum_cr_raw_eq0={delta_sum_cr_raw_eq0}\n")
        fh.write(f"pairs_cr_call_yes={pairs_cr_call_yes}\n")
        fh.write(f"pairs_cr_call_no={pairs_cr_call_no}\n")
        fh.write(f"delta_sum_cr_call_yes={delta_sum_cr_call_yes}\n")
        fh.write(f"delta_sum_cr_call_no={delta_sum_cr_call_no}\n")

    print(out_tsv)
    print(out_summary)


if __name__ == "__main__":
    main()
