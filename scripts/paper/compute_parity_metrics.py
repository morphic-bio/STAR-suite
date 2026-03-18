#!/usr/bin/env python3
"""
Compute parity metrics between STAR and Cell Ranger filtered MEX outputs.

Metrics:
  - Barcode Jaccard index
  - Per-barcode total-UMI Pearson and Spearman correlations (per feature type)
  - Per-feature total-UMI Pearson and Spearman correlations (per feature type)
  - CRISPR call concordance (from protospacer_calls_per_cell.csv)

Usage:
  python compute_parity_metrics.py \
    --cr-mex /path/to/cr/filtered_feature_bc_matrix \
    --star-mex /path/to/star/filtered_feature_bc_matrix \
    [--cr-crispr /path/to/cr/protospacer_calls_per_cell.csv] \
    [--star-crispr /path/to/star/protospacer_calls_per_cell.csv] \
    [--label DATASET_NAME] \
    [--tsv-out /path/to/output.tsv]
"""

import argparse
import csv
import gzip
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def open_gz(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def resolve(mex_dir: Path, basename: str) -> Path:
    for candidate in [mex_dir / basename, mex_dir / f"{basename}.gz"]:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Missing {basename}(.gz) in {mex_dir}")


def load_barcodes(mex_dir: Path, strip_suffix: bool = True) -> List[str]:
    path = resolve(mex_dir, "barcodes.tsv")
    barcodes = []
    with open_gz(path) as fh:
        for line in fh:
            bc = line.strip()
            if not bc:
                continue
            if strip_suffix and "-" in bc:
                base, suf = bc.rsplit("-", 1)
                if suf.isdigit():
                    bc = base
            barcodes.append(bc)
    return barcodes


def load_features(mex_dir: Path) -> List[Tuple[int, str, str, str]]:
    """Returns list of (1-based-row, id, name, type)."""
    path = resolve(mex_dir, "features.tsv")
    features = []
    with open_gz(path) as fh:
        for idx, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            fid = parts[0]
            fname = parts[1] if len(parts) > 1 else fid
            ftype = parts[2] if len(parts) > 2 else "Gene Expression"
            features.append((idx, fid, fname, ftype))
    return features


def load_matrix(mex_dir: Path, basename: str = "matrix.mtx"):
    """Yields (row_1based, col_1based, value)."""
    path = resolve(mex_dir, basename)
    with open_gz(path) as fh:
        header = fh.readline().strip()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"Invalid MatrixMarket: {header}")
        line = fh.readline().strip()
        while line.startswith("%"):
            line = fh.readline().strip()
        for raw in fh:
            raw = raw.strip()
            if not raw:
                continue
            parts = raw.split()
            if len(parts) < 3:
                continue
            yield int(parts[0]), int(parts[1]), float(parts[2])


def rankdata(values: List[float]) -> List[float]:
    pairs = sorted((v, i) for i, v in enumerate(values))
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        while j < len(values) and pairs[j][0] == pairs[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[pairs[k][1]] = avg_rank
        i = j
    return ranks


def pearson(xs: List[float], ys: List[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if dx == 0 or dy == 0:
        return float("nan")
    return num / (dx * dy)


def spearman(xs: List[float], ys: List[float]) -> float:
    if len(xs) < 2:
        return float("nan")
    return pearson(rankdata(xs), rankdata(ys))


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return 1.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union > 0 else 0.0


def load_mex_per_barcode_per_feature(
    mex_dir: Path,
    feature_type_filter: Optional[str] = None,
) -> Tuple[Dict[str, float], Dict[str, float], float, int, int]:
    """
    Returns (barcode_counts, feature_counts, total, n_barcodes, n_features).
    Barcodes are stripped of -N suffix.
    """
    barcodes = load_barcodes(mex_dir)
    features = load_features(mex_dir)

    if feature_type_filter:
        keep_rows = {row for row, fid, fname, ftype in features if ftype == feature_type_filter}
        fid_map = {row: fid for row, fid, fname, ftype in features if ftype == feature_type_filter}
    else:
        keep_rows = {row for row, fid, fname, ftype in features}
        fid_map = {row: fid for row, fid, fname, ftype in features}

    bc_counts: Dict[str, float] = defaultdict(float)
    feat_counts: Dict[str, float] = defaultdict(float)
    total = 0.0

    for row, col, val in load_matrix(mex_dir):
        if row not in keep_rows or val == 0:
            continue
        if col < 1 or col - 1 >= len(barcodes):
            continue
        bc = barcodes[col - 1]
        bc_counts[bc] += val
        feat_counts[fid_map[row]] += val
        total += val

    return dict(bc_counts), dict(feat_counts), total, len(barcodes), len(fid_map)


def load_crispr_calls(path: Path, strip_suffix: bool = True) -> Dict[str, dict]:
    out = {}
    with open_gz(path) as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            return out
        key_col = "cell_barcode" if "cell_barcode" in reader.fieldnames else reader.fieldnames[0]
        for row in reader:
            bc = (row.get(key_col, "") or "").strip()
            if not bc:
                continue
            if strip_suffix and "-" in bc:
                base, suf = bc.rsplit("-", 1)
                if suf.isdigit():
                    bc = base
            call = (row.get("feature_call", "") or "").strip()
            if not call or call.lower() == "none":
                call = "None"
            num_features = int(row.get("num_features", "0") or "0")
            umis_raw = (row.get("num_umis", "") or "").strip()
            umis = 0
            if umis_raw:
                for tok in umis_raw.split("|"):
                    try:
                        umis += int(tok)
                    except ValueError:
                        pass
            out[bc] = {"call": call, "num_features": num_features, "umis": umis}
    return out


def fmt(v: float, digits: int = 6) -> str:
    if math.isnan(v):
        return "NA"
    return f"{v:.{digits}f}"


def compute_correlation_metrics(
    left: Dict[str, float], right: Dict[str, float]
) -> dict:
    common = sorted(set(left) & set(right))
    xs = [left[k] for k in common]
    ys = [right[k] for k in common]
    return {
        "n_common": len(common),
        "n_left": len(left),
        "n_right": len(right),
        "pearson": pearson(xs, ys),
        "spearman": spearman(xs, ys),
        "left_sum": sum(xs),
        "right_sum": sum(ys),
    }


def compute_crispr_concordance(
    cr_calls: Dict[str, dict], star_calls: Dict[str, dict]
) -> dict:
    cr_bcs = set(cr_calls)
    star_bcs = set(star_calls)
    common = sorted(cr_bcs & star_bcs)

    if not common:
        return {
            "cr_cells": len(cr_calls),
            "star_cells": len(star_calls),
            "common_cells": 0,
            "exact_match": 0,
            "exact_match_pct": float("nan"),
            "umi_pearson": float("nan"),
            "umi_spearman": float("nan"),
        }

    exact = sum(1 for bc in common if cr_calls[bc]["call"] == star_calls[bc]["call"])
    cr_umis = [cr_calls[bc]["umis"] for bc in common]
    star_umis = [star_calls[bc]["umis"] for bc in common]

    return {
        "cr_cells": len(cr_calls),
        "star_cells": len(star_calls),
        "common_cells": len(common),
        "exact_match": exact,
        "exact_match_pct": 100.0 * exact / len(common),
        "umi_pearson": pearson([float(x) for x in cr_umis], [float(x) for x in star_umis]),
        "umi_spearman": spearman([float(x) for x in cr_umis], [float(x) for x in star_umis]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cr-mex", required=True, help="CR filtered_feature_bc_matrix directory")
    parser.add_argument("--star-mex", required=True, help="STAR filtered_feature_bc_matrix directory")
    parser.add_argument("--cr-crispr", default="", help="CR protospacer_calls_per_cell.csv")
    parser.add_argument("--star-crispr", default="", help="STAR protospacer_calls_per_cell.csv")
    parser.add_argument("--label", default="dataset", help="Dataset label for output")
    parser.add_argument("--tsv-out", default="", help="Write machine-readable TSV summary")
    parser.add_argument("--feature-types", default="Gene Expression,CRISPR Guide Capture",
                        help="Comma-separated feature types to analyze independently (default: GEX,CRISPR)")
    args = parser.parse_args()

    cr_mex = Path(args.cr_mex)
    star_mex = Path(args.star_mex)

    feature_types = [ft.strip() for ft in args.feature_types.split(",") if ft.strip()]
    feature_types_plus_all = ["all"] + feature_types

    results = {}

    print(f"{'='*70}")
    print(f"Parity Metrics: {args.label}")
    print(f"{'='*70}")
    print(f"CR MEX:   {cr_mex}")
    print(f"STAR MEX: {star_mex}")
    print()

    for ft_label in feature_types_plus_all:
        ft_filter = None if ft_label == "all" else ft_label
        tag = ft_label.replace(" ", "_")

        print(f"--- {ft_label} ---")

        try:
            cr_bc, cr_feat, cr_total, cr_n_bc, cr_n_feat = load_mex_per_barcode_per_feature(cr_mex, ft_filter)
            star_bc, star_feat, star_total, star_n_bc, star_n_feat = load_mex_per_barcode_per_feature(star_mex, ft_filter)
        except FileNotFoundError as e:
            print(f"  SKIP: {e}")
            print()
            continue

        cr_bc_set = set(cr_bc)
        star_bc_set = set(star_bc)
        j = jaccard(cr_bc_set, star_bc_set)

        bc_corr = compute_correlation_metrics(cr_bc, star_bc)
        feat_corr = compute_correlation_metrics(cr_feat, star_feat)

        print(f"  CR cells:       {cr_n_bc}")
        print(f"  STAR cells:     {star_n_bc}")
        print(f"  CR features:    {cr_n_feat}")
        print(f"  STAR features:  {star_n_feat}")
        print(f"  CR total UMIs:  {cr_total:,.0f}")
        print(f"  STAR total UMIs:{star_total:,.0f}")
        print()
        print(f"  Barcode Jaccard:           {fmt(j, 4)}")
        print(f"  Common barcodes:           {bc_corr['n_common']}")
        print(f"  Per-barcode Pearson:       {fmt(bc_corr['pearson'])}")
        print(f"  Per-barcode Spearman:      {fmt(bc_corr['spearman'])}")
        print()
        print(f"  Common features:           {feat_corr['n_common']}")
        print(f"  Per-feature Pearson:       {fmt(feat_corr['pearson'])}")
        print(f"  Per-feature Spearman:      {fmt(feat_corr['spearman'])}")
        print()

        results[tag] = {
            "barcode_jaccard": j,
            "bc_pearson": bc_corr["pearson"],
            "bc_spearman": bc_corr["spearman"],
            "feat_pearson": feat_corr["pearson"],
            "feat_spearman": feat_corr["spearman"],
            "cr_cells": cr_n_bc,
            "star_cells": star_n_bc,
            "common_barcodes": bc_corr["n_common"],
            "common_features": feat_corr["n_common"],
            "cr_total": cr_total,
            "star_total": star_total,
        }

    crispr_result = None
    if args.cr_crispr and args.star_crispr:
        cr_calls = load_crispr_calls(Path(args.cr_crispr))
        star_calls = load_crispr_calls(Path(args.star_crispr))
        crispr_result = compute_crispr_concordance(cr_calls, star_calls)

        print(f"--- CRISPR Call Concordance ---")
        print(f"  CR cells with calls:       {crispr_result['cr_cells']}")
        print(f"  STAR cells with calls:     {crispr_result['star_cells']}")
        print(f"  Common cells:              {crispr_result['common_cells']}")
        print(f"  Exact call match:          {crispr_result['exact_match']} / {crispr_result['common_cells']}")
        print(f"  Exact match %:             {fmt(crispr_result['exact_match_pct'], 2)}")
        print(f"  UMI Pearson:               {fmt(crispr_result['umi_pearson'])}")
        print(f"  UMI Spearman:              {fmt(crispr_result['umi_spearman'])}")
        print()

    if args.tsv_out:
        tsv_path = Path(args.tsv_out)
        with open(tsv_path, "w") as fh:
            cols = [
                "label", "feature_type",
                "barcode_jaccard", "common_barcodes",
                "bc_pearson", "bc_spearman",
                "common_features", "feat_pearson", "feat_spearman",
                "cr_cells", "star_cells", "cr_total_umis", "star_total_umis",
                "crispr_exact_match", "crispr_exact_match_pct",
                "crispr_umi_pearson", "crispr_umi_spearman",
            ]
            fh.write("\t".join(cols) + "\n")
            for tag, r in results.items():
                row = [
                    args.label, tag,
                    fmt(r["barcode_jaccard"], 6), str(r["common_barcodes"]),
                    fmt(r["bc_pearson"]), fmt(r["bc_spearman"]),
                    str(r["common_features"]), fmt(r["feat_pearson"]), fmt(r["feat_spearman"]),
                    str(r["cr_cells"]), str(r["star_cells"]),
                    f"{r['cr_total']:.0f}", f"{r['star_total']:.0f}",
                    str(crispr_result["exact_match"]) if crispr_result else "NA",
                    fmt(crispr_result["exact_match_pct"], 2) if crispr_result else "NA",
                    fmt(crispr_result["umi_pearson"]) if crispr_result else "NA",
                    fmt(crispr_result["umi_spearman"]) if crispr_result else "NA",
                ]
                fh.write("\t".join(row) + "\n")
        print(f"TSV written to {tsv_path}")


if __name__ == "__main__":
    main()
