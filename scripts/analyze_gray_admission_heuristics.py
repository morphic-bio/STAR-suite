#!/usr/bin/env python3
"""
Analyze cheap per-barcode heuristics for gray-zone admission decisions.

This consumes the UCSF gray-zone harness outputs and computes descriptive
statistics on the current/base raw matrix for gray-zone barcodes only.
"""

from __future__ import annotations

import argparse
import gzip
import math
import statistics
from pathlib import Path


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


def percentile(sorted_vals: list[float], p: float) -> float:
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    idx = (len(sorted_vals) - 1) * p
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return float(sorted_vals[lo])
    frac = idx - lo
    return float(sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac)


def summarize(values: list[float]) -> dict[str, str]:
    if not values:
        return {
            "n": "0",
            "mean": "NA",
            "median": "NA",
            "p25": "NA",
            "p75": "NA",
            "min": "NA",
            "max": "NA",
        }
    vals = sorted(values)
    return {
        "n": str(len(vals)),
        "mean": f"{statistics.fmean(vals):.6f}",
        "median": f"{statistics.median(vals):.6f}",
        "p25": f"{percentile(vals, 0.25):.6f}",
        "p75": f"{percentile(vals, 0.75):.6f}",
        "min": f"{vals[0]:.6f}",
        "max": f"{vals[-1]:.6f}",
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-raw-dir", type=Path, required=True, help="Base run raw GeneFull MEX dir")
    ap.add_argument("--gray-harness-outdir", type=Path, required=True, help="Output dir from evaluate_gray_em_admission.py")
    ap.add_argument("--base-filtered-barcodes", type=Path, required=True, help="Base/current filtered barcodes.tsv")
    ap.add_argument("--em-filtered-barcodes", type=Path, required=True, help="EM-like reference filtered barcodes.tsv")
    ap.add_argument("--cr-filtered-barcodes", type=Path, required=True, help="CR filtered barcodes.tsv(.gz)")
    ap.add_argument("--outdir", type=Path, required=True, help="Output directory")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    gray_zone = read_barcodes(args.gray_harness_outdir / "gray_zone_barcodes.tsv")
    base_filtered = read_barcodes(args.base_filtered_barcodes)
    em_filtered = read_barcodes(args.em_filtered_barcodes)
    cr_filtered = read_barcodes(args.cr_filtered_barcodes)

    base_gray = gray_zone & base_filtered
    em_gray = gray_zone & em_filtered

    groups = {
        "gray_base_and_em": base_gray & em_gray,
        "gray_base_only": base_gray - em_gray,
        "gray_em_only": em_gray - base_gray,
        "gray_neither": gray_zone - (base_gray | em_gray),
        "gray_base_only_in_cr": (base_gray - em_gray) & cr_filtered,
        "gray_base_only_not_in_cr": (base_gray - em_gray) - cr_filtered,
        "gray_em_only_in_cr": (em_gray - base_gray) & cr_filtered,
        "gray_em_only_not_in_cr": (em_gray - base_gray) - cr_filtered,
        "gray_base_and_em_in_cr": (base_gray & em_gray) & cr_filtered,
        "gray_base_and_em_not_in_cr": (base_gray & em_gray) - cr_filtered,
        "gray_neither_in_cr": (gray_zone - (base_gray | em_gray)) & cr_filtered,
        "gray_neither_not_in_cr": (gray_zone - (base_gray | em_gray)) - cr_filtered,
    }

    barcodes_path = args.base_raw_dir / "barcodes.tsv"
    matrix_path = args.base_raw_dir / "matrix.mtx"
    if not barcodes_path.exists() or not matrix_path.exists():
        raise SystemExit(f"Missing raw MEX under {args.base_raw_dir}")

    target = gray_zone
    metrics = {
        bc: {
            "total_umi": 0,
            "genes_detected": 0,
            "top_gene_count": 0,
            "sum_c_log_c": 0.0,
        }
        for bc in target
    }

    target_indices: dict[int, str] = {}
    with open(barcodes_path, "r", encoding="utf-8") as handle:
        for idx, raw in enumerate(handle, start=1):
            bc = normalize_barcode(raw)
            if bc in target:
                target_indices[idx] = bc

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
            col = int(parts[1])
            bc = target_indices.get(col)
            if bc is None:
                continue
            val = int(float(parts[2]))
            rec = metrics[bc]
            rec["total_umi"] += val
            rec["genes_detected"] += 1
            if val > rec["top_gene_count"]:
                rec["top_gene_count"] = val
            rec["sum_c_log_c"] += val * math.log(val)

    per_barcode_path = args.outdir / "per_barcode_metrics.tsv"
    with open(per_barcode_path, "w", encoding="utf-8") as out:
        out.write(
            "barcode\tgroup_base_em\tin_cr\ttotal_umi\tgenes_detected\tcomplexity_genes_per_umi\t"
            "top_gene_count\ttop_gene_frac\tentropy\n"
        )
        for bc in sorted(target):
            rec = metrics[bc]
            total = rec["total_umi"]
            genes = rec["genes_detected"]
            top = rec["top_gene_count"]
            complexity = (genes / total) if total else 0.0
            top_frac = (top / total) if total else 0.0
            entropy = (math.log(total) - (rec["sum_c_log_c"] / total)) if total else 0.0
            if bc in base_gray and bc in em_gray:
                label = "base_and_em"
            elif bc in base_gray:
                label = "base_only"
            elif bc in em_gray:
                label = "em_only"
            else:
                label = "neither"
            out.write(
                f"{bc}\t{label}\t{1 if bc in cr_filtered else 0}\t{total}\t{genes}\t"
                f"{complexity:.6f}\t{top}\t{top_frac:.6f}\t{entropy:.6f}\n"
            )

    feature_extractors = {
        "total_umi": lambda bc: float(metrics[bc]["total_umi"]),
        "genes_detected": lambda bc: float(metrics[bc]["genes_detected"]),
        "complexity_genes_per_umi": lambda bc: (
            float(metrics[bc]["genes_detected"]) / float(metrics[bc]["total_umi"])
            if metrics[bc]["total_umi"]
            else 0.0
        ),
        "top_gene_count": lambda bc: float(metrics[bc]["top_gene_count"]),
        "top_gene_frac": lambda bc: (
            float(metrics[bc]["top_gene_count"]) / float(metrics[bc]["total_umi"])
            if metrics[bc]["total_umi"]
            else 0.0
        ),
        "entropy": lambda bc: (
            math.log(metrics[bc]["total_umi"]) - (metrics[bc]["sum_c_log_c"] / metrics[bc]["total_umi"])
            if metrics[bc]["total_umi"]
            else 0.0
        ),
    }

    summary_path = args.outdir / "heuristic_summary.tsv"
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write("group\tfeature\tn\tmean\tmedian\tp25\tp75\tmin\tmax\n")
        for group_name, barcode_set in groups.items():
            for feature_name, fn in feature_extractors.items():
                vals = [fn(bc) for bc in barcode_set]
                row = summarize(vals)
                out.write(
                    f"{group_name}\t{feature_name}\t{row['n']}\t{row['mean']}\t{row['median']}\t"
                    f"{row['p25']}\t{row['p75']}\t{row['min']}\t{row['max']}\n"
                )

    print(summary_path.read_text(), end="")
    print(f"Wrote per-barcode metrics to: {per_barcode_path}")


if __name__ == "__main__":
    main()
