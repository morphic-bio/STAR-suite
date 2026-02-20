#!/usr/bin/env python3
"""
Compare feature MEX outputs between Cell Ranger and STAR assignBarcodes.

Supports optional barcode translation (e.g., CellPlex NXT -> TRU).
Computes per-barcode and per-feature counts for a selected feature type.
Supports integer and real-valued MatrixMarket files (for multimapper outputs).
"""

import argparse
import gzip
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, Tuple


def open_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def resolve_mex_file(mex_dir: Path, basename: str) -> Path:
    plain = mex_dir / basename
    gz = mex_dir / f"{basename}.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise FileNotFoundError(f"Missing {basename}(.gz) in {mex_dir}")


def read_features(features_path: Path) -> Tuple[list, set]:
    features = []
    with open_maybe_gz(str(features_path)) as f:
        for idx, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split(",")
            feature_id = parts[0]
            feature_name = parts[1] if len(parts) > 1 else feature_id
            feature_type = parts[2] if len(parts) > 2 else "Gene Expression"
            features.append({
                "row": idx,
                "id": feature_id,
                "name": feature_name,
                "type": feature_type,
            })
    return features, {f["row"] for f in features}


def read_barcodes(barcodes_path: Path, strip_suffix: bool) -> list:
    barcodes = []
    with open_maybe_gz(str(barcodes_path)) as f:
        for line in f:
            bc = line.strip()
            if not bc:
                continue
            if strip_suffix and "-" in bc:
                base, suffix = bc.rsplit("-", 1)
                if suffix.isdigit():
                    bc = base
            barcodes.append(bc)
    return barcodes


def iter_matrix_entries(matrix_path: Path) -> Iterable[Tuple[int, int, float]]:
    with open_maybe_gz(str(matrix_path)) as f:
        header = f.readline().strip()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"Invalid MatrixMarket header: {header}")
        line = f.readline().strip()
        while line.startswith("%"):
            line = f.readline().strip()
        if not line:
            raise ValueError("Missing matrix dimensions")
        for raw in f:
            raw = raw.strip()
            if not raw:
                continue
            parts = raw.split()
            if len(parts) < 3:
                continue
            row = int(parts[0])
            col = int(parts[1])
            val = float(parts[2])
            yield row, col, val


def load_mex_counts(
    mex_dir: Path,
    feature_types: set,
    strip_suffix: bool,
    matrix_basename: str = "matrix.mtx",
) -> Tuple[list, Dict[str, float], Dict[str, float], float]:
    features_path = None
    for name in ("features.tsv", "features.txt"):
        try:
            features_path = resolve_mex_file(mex_dir, name)
            break
        except FileNotFoundError:
            continue
    if features_path is None:
        raise FileNotFoundError(f"Missing features.tsv(.gz) or features.txt(.gz) in {mex_dir}")

    barcodes_path = None
    for name in ("barcodes.tsv", "barcodes.txt"):
        try:
            barcodes_path = resolve_mex_file(mex_dir, name)
            break
        except FileNotFoundError:
            continue
    if barcodes_path is None:
        raise FileNotFoundError(f"Missing barcodes.tsv(.gz) or barcodes.txt(.gz) in {mex_dir}")

    matrix_path = resolve_mex_file(mex_dir, matrix_basename)

    features, _ = read_features(features_path)
    if feature_types:
        keep_rows = {f["row"] for f in features if f["type"] in feature_types}
        kept_features = [f for f in features if f["row"] in keep_rows]
    else:
        keep_rows = {f["row"] for f in features}
        kept_features = features
    row_to_id = {f["row"]: f["id"] for f in kept_features}

    barcodes = read_barcodes(barcodes_path, strip_suffix)

    barcode_counts: Dict[str, float] = defaultdict(float)
    feature_counts: Dict[str, float] = defaultdict(float)
    total_counts = 0.0

    for row, col, val in iter_matrix_entries(matrix_path):
        if row not in keep_rows or val == 0.0:
            continue
        if col - 1 >= len(barcodes) or col <= 0:
            continue
        bc = barcodes[col - 1]
        barcode_counts[bc] += val
        feature_id = row_to_id.get(row)
        if feature_id:
            feature_counts[feature_id] += val
        total_counts += val

    return kept_features, barcode_counts, feature_counts, total_counts


def load_translation_map(path: Path, direction: str, wanted: set) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open_maybe_gz(str(path)) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            left, right = parts[0], parts[1]
            if direction == "right-to-left":
                src, dst = right, left
            else:
                src, dst = left, right
            if src in wanted:
                mapping[src] = dst
    return mapping


def apply_translation(counts: Dict[str, float], mapping: Dict[str, str]) -> Dict[str, float]:
    out: Dict[str, float] = defaultdict(float)
    for bc, count in counts.items():
        new_bc = mapping.get(bc, bc)
        out[new_bc] += count
    return out


def fmt_num(value: float) -> str:
    rounded = round(value)
    if abs(value - rounded) < 1e-9:
        return str(int(rounded))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def summarize(label: str, features: list, barcode_counts: Dict[str, float], total_counts: float):
    print(f"{label}:")
    print(f"  Features: {len(features)}")
    print(f"  Barcodes with counts: {len(barcode_counts)}")
    print(f"  Total counts: {fmt_num(total_counts)}")


def load_barcode_filter(path: Path, strip_suffix: bool) -> set:
    barcodes = set()
    with open_maybe_gz(str(path)) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # accept CSV like "GRCh38,AAAC...-1" or plain barcode
            if "," in line:
                line = line.split(",")[-1]
            if strip_suffix and "-" in line:
                base, suffix = line.rsplit("-", 1)
                if suffix.isdigit():
                    line = base
            if line:
                barcodes.add(line)
    return barcodes


def main():
    parser = argparse.ArgumentParser(
        description="Compare feature MEX outputs between Cell Ranger and STAR assignBarcodes."
    )
    parser.add_argument("cr_mex", help="Cell Ranger MEX directory")
    parser.add_argument("star_mex", help="STAR assignBarcodes MEX directory")
    parser.add_argument(
        "--feature-type",
        default="Multiplexing Capture",
        help="Feature type to compare (default: Multiplexing Capture)",
    )
    parser.add_argument(
        "--matrix-basename",
        default="matrix.mtx",
        help=(
            "Matrix basename to load from each MEX directory (default: matrix.mtx). "
            "Examples: matrix.mtx, UniqueAndMult-Rescue.mtx"
        ),
    )
    parser.add_argument(
        "--matrix-basename-a",
        default="",
        help=(
            "Optional matrix basename for first MEX directory only. "
            "If unset, --matrix-basename is used."
        ),
    )
    parser.add_argument(
        "--matrix-basename-b",
        default="",
        help=(
            "Optional matrix basename for second MEX directory only. "
            "If unset, --matrix-basename is used."
        ),
    )
    parser.add_argument(
        "--barcode-translation",
        default="",
        help="Two-column translation file (optional)",
    )
    parser.add_argument(
        "--translation-direction",
        default="left-to-right",
        choices=("left-to-right", "right-to-left"),
        help="Translation direction (default: left-to-right)",
    )
    parser.add_argument(
        "--translate",
        default="star",
        choices=("star", "cr", "both", "none"),
        help="Which side to translate (default: star)",
    )
    parser.add_argument(
        "--keep-suffix",
        action="store_true",
        help="Keep barcode suffixes like -1 (default: strip)",
    )
    parser.add_argument(
        "--barcode-filter",
        default="",
        help="Optional barcode list to filter counts (CSV or one barcode per line).",
    )
    parser.add_argument(
        "--barcode-filter-side",
        default="both",
        choices=("cr", "star", "both"),
        help="Which side to apply barcode filter to (default: both)",
    )
    args = parser.parse_args()

    cr_dir = Path(args.cr_mex)
    star_dir = Path(args.star_mex)
    if not cr_dir.exists():
        raise SystemExit(f"CR MEX dir not found: {cr_dir}")
    if not star_dir.exists():
        raise SystemExit(f"STAR MEX dir not found: {star_dir}")

    feature_types = set()
    if args.feature_type and args.feature_type.lower() not in ("all", "*", "any", "none"):
        feature_types = {args.feature_type}

    strip_suffix = not args.keep_suffix

    matrix_basename_a = args.matrix_basename_a if args.matrix_basename_a else args.matrix_basename
    matrix_basename_b = args.matrix_basename_b if args.matrix_basename_b else args.matrix_basename

    cr_features, cr_bc_counts, cr_feat_counts, cr_total = load_mex_counts(
        cr_dir, feature_types, strip_suffix, matrix_basename_a
    )
    star_features, star_bc_counts, star_feat_counts, star_total = load_mex_counts(
        star_dir, feature_types, strip_suffix, matrix_basename_b
    )

    # Apply translation mapping if provided
    if args.barcode_translation and args.translate != "none":
        translation_path = Path(args.barcode_translation)
        if not translation_path.exists():
            raise SystemExit(f"Translation file not found: {translation_path}")

        if args.translate in ("cr", "both"):
            wanted = set(cr_bc_counts.keys())
            mapping = load_translation_map(translation_path, args.translation_direction, wanted)
            cr_bc_counts = apply_translation(cr_bc_counts, mapping)

        if args.translate in ("star", "both"):
            wanted = set(star_bc_counts.keys())
            mapping = load_translation_map(translation_path, args.translation_direction, wanted)
            star_bc_counts = apply_translation(star_bc_counts, mapping)

    # Optional barcode filter (applied after translation)
    if args.barcode_filter:
        filter_path = Path(args.barcode_filter)
        if not filter_path.exists():
            raise SystemExit(f"Barcode filter file not found: {filter_path}")
        filter_set = load_barcode_filter(filter_path, strip_suffix)
        if args.barcode_filter_side in ("cr", "both"):
            cr_bc_counts = {bc: c for bc, c in cr_bc_counts.items() if bc in filter_set}
        if args.barcode_filter_side in ("star", "both"):
            star_bc_counts = {bc: c for bc, c in star_bc_counts.items() if bc in filter_set}

    # Recompute totals after translation/filtering
    cr_total = sum(cr_bc_counts.values())
    star_total = sum(star_bc_counts.values())

    print("=" * 70)
    print("Feature MEX Comparison")
    print("=" * 70)
    print()
    summarize("Cell Ranger", cr_features, cr_bc_counts, cr_total)
    summarize("STAR assignBarcodes", star_features, star_bc_counts, star_total)
    print()

    cr_bcs = set(cr_bc_counts.keys())
    star_bcs = set(star_bc_counts.keys())
    common = cr_bcs & star_bcs

    print("Barcode overlap:")
    print(f"  CR unique barcodes: {len(cr_bcs)}")
    print(f"  STAR unique barcodes: {len(star_bcs)}")
    print(f"  Common barcodes: {len(common)}")

    if common:
        cr_common_total = sum(cr_bc_counts[bc] for bc in common)
        star_common_total = sum(star_bc_counts[bc] for bc in common)
        print(f"  CR counts on common barcodes: {fmt_num(cr_common_total)}")
        print(f"  STAR counts on common barcodes: {fmt_num(star_common_total)}")
        print(f"  Delta (STAR-CR): {fmt_num(star_common_total - cr_common_total)}")

    print()
    cr_feat_ids = {f["id"] for f in cr_features}
    star_feat_ids = {f["id"] for f in star_features}
    common_feats = cr_feat_ids & star_feat_ids
    print("Feature overlap:")
    print(f"  CR features: {len(cr_feat_ids)}")
    print(f"  STAR features: {len(star_feat_ids)}")
    print(f"  Common features: {len(common_feats)}")


if __name__ == "__main__":
    main()
