#!/usr/bin/env python3
"""
Evaluate a simple post-EmptyDrops QC filter applied only to accepted gray-zone cells.

This mirrors the downstream perturb/scRNA post-processing semantics used via
combineFilters.py:
- genes_detected >= min_genes
- genes_detected <= max_genes
- mt_pct <= mt_pct_cutoff

Optionally, max_genes can be made adaptive using the same median + N*MAD logic
used by scripts/compute_adaptive_qc_threshold.py.
"""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path
import statistics


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


def read_barcodes(path: Path) -> set[str]:
    out = set()
    with open_maybe_gz(path) as handle:
        for raw in handle:
            bc = normalize_barcode(raw)
            if bc:
                out.add(bc)
    return out


def summarize_median_mad(values: list[float]) -> tuple[float, float]:
    if len(values) < 1:
        raise SystemExit("Need at least one value to estimate median/MAD")
    median = statistics.median(values)
    abs_dev = [abs(v - median) for v in values]
    mad = statistics.median(abs_dev)
    return float(median), float(mad)


def load_features(features_path: Path) -> list[bool]:
    mito = []
    with open(features_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            parts = raw.rstrip("\n").split("\t")
            gene_name = parts[1] if len(parts) > 1 else parts[0]
            mito.append(gene_name.upper().startswith("MT-"))
    return mito


def compute_metrics(raw_dir: Path, target_barcodes: set[str]) -> dict[str, dict[str, float]]:
    barcodes_path = raw_dir / "barcodes.tsv"
    features_path = raw_dir / "features.tsv"
    matrix_path = raw_dir / "matrix.mtx"
    if not barcodes_path.exists() or not features_path.exists() or not matrix_path.exists():
        raise SystemExit(f"Missing raw MEX files under {raw_dir}")

    mito_feature = load_features(features_path)
    target_cols: dict[int, str] = {}
    metrics = {
        bc: {
            "total_umi": 0.0,
            "genes_detected": 0.0,
            "mt_umi": 0.0,
        }
        for bc in target_barcodes
    }

    with open(barcodes_path, "r", encoding="utf-8") as handle:
        for idx, raw in enumerate(handle, start=1):
            bc = normalize_barcode(raw)
            if bc in target_barcodes:
                target_cols[idx] = bc

    with open(matrix_path, "r", encoding="utf-8") as handle:
        dims_seen = False
        for raw in handle:
            if raw.startswith("%"):
                continue
            if not dims_seen:
                dims_seen = True
                continue
            parts = raw.split()
            if len(parts) != 3:
                continue
            row = int(parts[0]) - 1
            col = int(parts[1])
            bc = target_cols.get(col)
            if bc is None:
                continue
            val = int(float(parts[2]))
            rec = metrics[bc]
            rec["total_umi"] += val
            rec["genes_detected"] += 1.0
            if 0 <= row < len(mito_feature) and mito_feature[row]:
                rec["mt_umi"] += val

    for bc, rec in metrics.items():
        total = rec["total_umi"]
        rec["mt_frac"] = (rec["mt_umi"] / total) if total else 0.0
    return metrics


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-raw-dir", type=Path, required=True, help="Base run raw GeneFull MEX dir")
    ap.add_argument("--base-filtered-barcodes", type=Path, required=True, help="Base/current filtered barcodes.tsv")
    ap.add_argument("--gray-zone-barcodes", type=Path, required=True, help="Gray-zone barcodes.tsv")
    ap.add_argument("--cr-filtered-barcodes", type=Path, required=True, help="CR filtered barcodes.tsv(.gz)")
    ap.add_argument("--min-genes", type=int, default=200, help="Minimum genes cutoff")
    ap.add_argument("--max-genes", type=int, default=2500, help="Maximum genes cutoff when not using adaptive mode")
    ap.add_argument("--mt-pct-cutoff", type=float, default=5.0, help="Maximum mitochondrial percentage cutoff")
    ap.add_argument("--adaptive-max-genes", action="store_true", help="Use adaptive max_genes = median + N*MAD from accepted non-gray cells")
    ap.add_argument("--n-mad", type=float, default=3.0, help="MAD multiplier for adaptive max_genes")
    ap.add_argument("--outdir", type=Path, required=True, help="Output directory")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    base_filtered = read_barcodes(args.base_filtered_barcodes)
    gray_zone = read_barcodes(args.gray_zone_barcodes)
    cr_filtered = read_barcodes(args.cr_filtered_barcodes)

    accepted_gray = base_filtered & gray_zone
    accepted_non_gray = base_filtered - gray_zone
    if not accepted_gray:
        raise SystemExit("No accepted gray-zone cells found")
    if not accepted_non_gray:
        raise SystemExit("No accepted non-gray cells found")

    metrics = compute_metrics(args.base_raw_dir, base_filtered)

    genes_non_gray = [metrics[bc]["genes_detected"] for bc in accepted_non_gray]
    mt_non_gray = [metrics[bc]["mt_frac"] for bc in accepted_non_gray]
    genes_median, genes_mad = summarize_median_mad(genes_non_gray)
    if args.adaptive_max_genes:
        genes_max = max(int(genes_median + args.n_mad * genes_mad), int(args.min_genes))
    else:
        genes_max = int(args.max_genes)
    genes_min = int(args.min_genes)
    mt_max = float(args.mt_pct_cutoff) / 100.0

    removed = set()
    for bc in accepted_gray:
        rec = metrics[bc]
        if rec["genes_detected"] < genes_min or rec["genes_detected"] > genes_max or rec["mt_frac"] > mt_max:
            removed.add(bc)

    filtered_postqc = base_filtered - removed

    base_common = len(base_filtered & cr_filtered)
    base_union = len(base_filtered | cr_filtered)
    post_common = len(filtered_postqc & cr_filtered)
    post_union = len(filtered_postqc | cr_filtered)

    summary_path = args.outdir / "summary.tsv"
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write("metric\tvalue\n")
        out.write(f"accepted_gray\t{len(accepted_gray)}\n")
        out.write(f"accepted_non_gray\t{len(accepted_non_gray)}\n")
        out.write(f"adaptive_max_genes\t{1 if args.adaptive_max_genes else 0}\n")
        out.write(f"n_mad\t{args.n_mad:.6f}\n")
        out.write(f"genes_median_non_gray\t{genes_median:.6f}\n")
        out.write(f"genes_mad_non_gray\t{genes_mad:.6f}\n")
        out.write(f"genes_min_cutoff\t{genes_min}\n")
        out.write(f"genes_max_cutoff\t{genes_max}\n")
        out.write(f"mt_pct_cutoff\t{args.mt_pct_cutoff:.6f}\n")
        out.write(f"mt_max_frac_cutoff\t{mt_max:.6f}\n")
        out.write(f"non_gray_mt_mean_pct\t{statistics.fmean(mt_non_gray) * 100.0:.6f}\n")
        out.write(f"gray_removed\t{len(removed)}\n")
        out.write(f"gray_removed_in_cr\t{len(removed & cr_filtered)}\n")
        out.write(f"gray_removed_not_in_cr\t{len(removed - cr_filtered)}\n")
        out.write(f"base_filtered_cells\t{len(base_filtered)}\n")
        out.write(f"base_common_cr\t{base_common}\n")
        out.write(f"base_union_cr\t{base_union}\n")
        out.write(f"base_jaccard_cr\t{(base_common / base_union) if base_union else 1.0:.6f}\n")
        out.write(f"postqc_filtered_cells\t{len(filtered_postqc)}\n")
        out.write(f"postqc_common_cr\t{post_common}\n")
        out.write(f"postqc_union_cr\t{post_union}\n")
        out.write(f"postqc_jaccard_cr\t{(post_common / post_union) if post_union else 1.0:.6f}\n")

    with open(args.outdir / "removed_gray_barcodes.tsv", "w", encoding="utf-8") as out:
        for bc in sorted(removed):
            out.write(f"{bc}\n")

    with open(args.outdir / "removed_gray_in_cr.tsv", "w", encoding="utf-8") as out:
        for bc in sorted(removed & cr_filtered):
            out.write(f"{bc}\n")

    with open(args.outdir / "per_barcode_qc.tsv", "w", encoding="utf-8") as out:
        out.write("barcode\tis_gray\tin_cr\ttotal_umi\tgenes_detected\tmt_umi\tmt_frac\tremoved_by_postqc\n")
        for bc in sorted(base_filtered):
            rec = metrics[bc]
            out.write(
                f"{bc}\t{1 if bc in accepted_gray else 0}\t{1 if bc in cr_filtered else 0}\t"
                f"{rec['total_umi']:.0f}\t{rec['genes_detected']:.0f}\t{rec['mt_umi']:.0f}\t"
                f"{rec['mt_frac'] * 100.0:.6f}\t{1 if bc in removed else 0}\n"
            )

    print(summary_path.read_text(), end="")
    print(f"Wrote gray-zone post-QC evaluation to: {args.outdir}")


if __name__ == "__main__":
    main()
