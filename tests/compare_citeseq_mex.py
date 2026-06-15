#!/usr/bin/env python3
"""Compare STAR-suite ADT/HTO feature MEX against a Cell Ranger MEX oracle."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def open_maybe_gz(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def resolve_mex_file(mex_dir: Path, basename: str) -> Path:
    plain = mex_dir / basename
    if plain.exists():
        return plain
    gz = mex_dir / f"{basename}.gz"
    if gz.exists():
        return gz
    raise FileNotFoundError(f"Missing {basename}(.gz) in {mex_dir}")


def strip_barcode_suffix(barcode: str) -> str:
    if "-" not in barcode:
        return barcode
    base, suffix = barcode.rsplit("-", 1)
    if suffix.isdigit():
        return base
    return barcode


def translate_nxt_tru_barcode(barcode: str) -> str:
    """Apply STAR-suite's reversible NXT/TRU two-base barcode transform."""
    if len(barcode) < 9:
        return barcode
    bases = list(barcode)
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    bases[7] = bases[7].translate(complement)
    bases[8] = bases[8].translate(complement)
    return "".join(bases)


def transform_mex_barcodes(mex: dict, transform: str) -> dict:
    if transform == "none":
        return mex
    if transform != "nxt_tru":
        raise ValueError(f"Unsupported barcode transform: {transform}")

    barcode_totals: Dict[str, float] = defaultdict(float)
    entries: Dict[Tuple[str, str], float] = defaultdict(float)
    transformed_barcodes = []
    seen = set()

    for barcode in mex["barcodes"]:
        translated = translate_nxt_tru_barcode(barcode)
        if translated not in seen:
            transformed_barcodes.append(translated)
            seen.add(translated)

    for barcode, total in mex["barcode_totals"].items():
        barcode_totals[translate_nxt_tru_barcode(barcode)] += total

    for (feature_id, barcode), value in mex["entries"].items():
        entries[(feature_id, translate_nxt_tru_barcode(barcode))] += value

    mex = dict(mex)
    mex["barcodes"] = transformed_barcodes
    mex["barcode_totals"] = dict(barcode_totals)
    mex["entries"] = dict(entries)
    return mex


def read_features(mex_dir: Path) -> List[dict]:
    path = resolve_mex_file(mex_dir, "features.tsv")
    features = []
    with open_maybe_gz(path) as handle:
        for row, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split(",")
            feature_id = parts[0]
            feature_name = parts[1] if len(parts) > 1 else feature_id
            feature_type = parts[2] if len(parts) > 2 else "Gene Expression"
            features.append(
                {
                    "row": row,
                    "id": feature_id,
                    "name": feature_name,
                    "type": feature_type,
                }
            )
    return features


def read_barcodes(mex_dir: Path, keep_suffix: bool) -> List[str]:
    path = resolve_mex_file(mex_dir, "barcodes.tsv")
    barcodes = []
    with open_maybe_gz(path) as handle:
        for raw in handle:
            barcode = raw.strip()
            if not barcode:
                continue
            if not keep_suffix:
                barcode = strip_barcode_suffix(barcode)
            barcodes.append(barcode)
    return barcodes


def iter_matrix_entries(mex_dir: Path) -> Iterable[Tuple[int, int, float]]:
    path = resolve_mex_file(mex_dir, "matrix.mtx")
    with open_maybe_gz(path) as handle:
        header = handle.readline().strip()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"Invalid MatrixMarket header in {path}: {header}")
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            # First non-comment record after the header is the matrix shape.
            break
        else:
            raise ValueError(f"Missing MatrixMarket shape in {path}")

        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            yield int(parts[0]), int(parts[1]), float(parts[2])


def load_selected_mex(mex_dir: Path, feature_type: str, keep_suffix: bool) -> dict:
    features = read_features(mex_dir)
    barcodes = read_barcodes(mex_dir, keep_suffix)

    if feature_type.lower() in {"all", "*", "any"}:
        selected = features
    else:
        selected = [feature for feature in features if feature["type"] == feature_type]

    row_to_feature = {feature["row"]: feature["id"] for feature in selected}
    feature_ids = [feature["id"] for feature in selected]
    feature_names = {feature["id"]: feature["name"] for feature in selected}
    feature_totals: Dict[str, float] = {feature_id: 0.0 for feature_id in feature_ids}
    barcode_totals: Dict[str, float] = {barcode: 0.0 for barcode in barcodes}
    entries: Dict[Tuple[str, str], float] = defaultdict(float)

    for row, col, value in iter_matrix_entries(mex_dir):
        feature_id = row_to_feature.get(row)
        if feature_id is None or value == 0.0:
            continue
        if col <= 0 or col > len(barcodes):
            continue
        barcode = barcodes[col - 1]
        feature_totals[feature_id] = feature_totals.get(feature_id, 0.0) + value
        barcode_totals[barcode] = barcode_totals.get(barcode, 0.0) + value
        entries[(feature_id, barcode)] += value

    return {
        "mex_dir": str(mex_dir),
        "features": selected,
        "feature_ids": feature_ids,
        "feature_names": feature_names,
        "barcodes": barcodes,
        "feature_totals": feature_totals,
        "barcode_totals": barcode_totals,
        "entries": dict(entries),
        "total": sum(feature_totals.values()),
    }


def pearson(xs: List[float], ys: List[float]) -> float | None:
    if len(xs) != len(ys) or not xs:
        return None
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if den_x == 0.0 or den_y == 0.0:
        return 1.0 if all(abs(x - y) < 1e-9 for x, y in zip(xs, ys)) else 0.0
    return num / (den_x * den_y)


def rel_delta(a: float, b: float) -> float:
    denom = max(abs(a), 1.0)
    return abs(b - a) / denom


def write_feature_totals(path: Path, cr: dict, star: dict, common_features: List[str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("feature_id\tfeature_name\tcellranger_total\tstar_total\tdelta\n")
        for feature_id in common_features:
            cr_total = cr["feature_totals"].get(feature_id, 0.0)
            star_total = star["feature_totals"].get(feature_id, 0.0)
            name = cr["feature_names"].get(feature_id) or star["feature_names"].get(feature_id) or feature_id
            handle.write(f"{feature_id}\t{name}\t{cr_total:g}\t{star_total:g}\t{star_total - cr_total:g}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cr-mex", required=True, help="Cell Ranger MEX directory")
    parser.add_argument("--star-mex", required=True, help="STAR-suite MEX directory")
    parser.add_argument("--feature-type", default="Antibody Capture", help="Feature type to compare")
    parser.add_argument("--keep-suffix", action="store_true", help="Keep barcode suffixes like -1")
    parser.add_argument(
        "--translate-star-nxt-tru",
        action="store_true",
        help="Translate STAR barcodes between NXT/TRU namespaces before comparison",
    )
    parser.add_argument(
        "--translate-cr-nxt-tru",
        action="store_true",
        help="Translate Cell Ranger barcodes between NXT/TRU namespaces before comparison",
    )
    parser.add_argument("--report-json", default="", help="Optional JSON report path")
    parser.add_argument("--feature-totals-tsv", default="", help="Optional per-feature total TSV path")
    parser.add_argument("--min-feature-pearson", type=float, default=0.95)
    parser.add_argument("--min-barcode-pearson", type=float, default=0.90)
    parser.add_argument("--max-total-rel-delta", type=float, default=0.25)
    parser.add_argument("--max-missing-feature-frac", type=float, default=0.0)
    args = parser.parse_args()

    cr_transform = "nxt_tru" if args.translate_cr_nxt_tru else "none"
    star_transform = "nxt_tru" if args.translate_star_nxt_tru else "none"
    cr = transform_mex_barcodes(
        load_selected_mex(Path(args.cr_mex), args.feature_type, args.keep_suffix),
        cr_transform,
    )
    star = transform_mex_barcodes(
        load_selected_mex(Path(args.star_mex), args.feature_type, args.keep_suffix),
        star_transform,
    )

    cr_features = set(cr["feature_ids"])
    star_features = set(star["feature_ids"])
    common_features = sorted(cr_features & star_features)
    missing_from_star = sorted(cr_features - star_features)
    extra_in_star = sorted(star_features - cr_features)

    cr_barcodes = set(cr["barcodes"])
    star_barcodes = set(star["barcodes"])
    common_barcodes = sorted(cr_barcodes & star_barcodes)

    feature_corr = pearson(
        [cr["feature_totals"].get(feature_id, 0.0) for feature_id in common_features],
        [star["feature_totals"].get(feature_id, 0.0) for feature_id in common_features],
    )
    barcode_corr = pearson(
        [cr["barcode_totals"].get(barcode, 0.0) for barcode in common_barcodes],
        [star["barcode_totals"].get(barcode, 0.0) for barcode in common_barcodes],
    )

    common_entry_keys = {
        key
        for key in (set(cr["entries"]) | set(star["entries"]))
        if key[0] in common_features and key[1] in common_barcodes
    }
    abs_deltas = [
        abs(star["entries"].get(key, 0.0) - cr["entries"].get(key, 0.0))
        for key in common_entry_keys
    ]
    nonzero_entry_diffs = sum(1 for value in abs_deltas if value > 1e-9)

    total_rel_delta = rel_delta(cr["total"], star["total"])
    missing_feature_frac = (
        len(missing_from_star) / len(cr_features)
        if cr_features
        else 1.0
    )

    report = {
        "feature_type": args.feature_type,
        "cellranger_mex": cr["mex_dir"],
        "star_mex": star["mex_dir"],
        "cellranger_barcode_transform": cr_transform,
        "star_barcode_transform": star_transform,
        "cellranger_features": len(cr_features),
        "star_features": len(star_features),
        "common_features": len(common_features),
        "missing_features_from_star": missing_from_star,
        "extra_features_in_star": extra_in_star,
        "cellranger_barcodes": len(cr_barcodes),
        "star_barcodes": len(star_barcodes),
        "common_barcodes": len(common_barcodes),
        "cellranger_total": cr["total"],
        "star_total": star["total"],
        "total_rel_delta": total_rel_delta,
        "feature_total_pearson": feature_corr,
        "barcode_total_pearson": barcode_corr,
        "common_entry_count": len(common_entry_keys),
        "common_entry_abs_delta_sum": sum(abs_deltas),
        "common_entry_max_abs_delta": max(abs_deltas) if abs_deltas else 0.0,
        "common_entry_nonzero_differences": nonzero_entry_diffs,
        "thresholds": {
            "min_feature_pearson": args.min_feature_pearson,
            "min_barcode_pearson": args.min_barcode_pearson,
            "max_total_rel_delta": args.max_total_rel_delta,
            "max_missing_feature_frac": args.max_missing_feature_frac,
        },
    }

    failures = []
    if not cr_features:
        failures.append(f"Cell Ranger MEX has no {args.feature_type!r} features")
    if not star_features:
        failures.append(f"STAR MEX has no {args.feature_type!r} features")
    if not common_features:
        failures.append("No common selected features")
    if not common_barcodes:
        failures.append("No common barcodes")
    if missing_feature_frac > args.max_missing_feature_frac:
        failures.append(
            f"Missing feature fraction {missing_feature_frac:.4f} exceeds {args.max_missing_feature_frac:.4f}"
        )
    if feature_corr is None or feature_corr < args.min_feature_pearson:
        failures.append(
            f"Feature-total Pearson {feature_corr} below {args.min_feature_pearson:.4f}"
        )
    if barcode_corr is None or barcode_corr < args.min_barcode_pearson:
        failures.append(
            f"Barcode-total Pearson {barcode_corr} below {args.min_barcode_pearson:.4f}"
        )
    if total_rel_delta > args.max_total_rel_delta:
        failures.append(
            f"Total relative delta {total_rel_delta:.4f} exceeds {args.max_total_rel_delta:.4f}"
        )

    report["pass"] = not failures
    report["failures"] = failures

    print("CITE-seq feature MEX comparison")
    print(f"  feature_type: {args.feature_type}")
    print(f"  barcode_transforms: Cell Ranger={cr_transform} STAR={star_transform}")
    print(f"  features: Cell Ranger={len(cr_features)} STAR={len(star_features)} common={len(common_features)}")
    print(f"  barcodes: Cell Ranger={len(cr_barcodes)} STAR={len(star_barcodes)} common={len(common_barcodes)}")
    print(f"  total counts: Cell Ranger={cr['total']:.0f} STAR={star['total']:.0f} rel_delta={total_rel_delta:.4f}")
    print(f"  feature_total_pearson: {feature_corr}")
    print(f"  barcode_total_pearson: {barcode_corr}")
    print(f"  entry_abs_delta_sum: {report['common_entry_abs_delta_sum']:.0f}")

    if args.report_json:
        path = Path(args.report_json)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"  wrote: {path}")

    if args.feature_totals_tsv:
        path = Path(args.feature_totals_tsv)
        path.parent.mkdir(parents=True, exist_ok=True)
        write_feature_totals(path, cr, star, common_features)
        print(f"  wrote: {path}")

    if failures:
        print("FAIL:")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
