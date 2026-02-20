#!/usr/bin/env python3
"""Export publication-ready miss tables for STAR vs CR (0 and 1 Hamming analyses)."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Export combined H0/H1 miss tables for manuscript use.")
    parser.add_argument("--analysis-dir", required=True, help="Output directory from run_ucsf_star1h_cr_analysis.sh")
    parser.add_argument("--out-dir", required=True, help="Output directory for paper tables")
    parser.add_argument("--top-n", type=int, default=25, help="Top-N rows per mode in markdown preview")
    args = parser.parse_args()

    analysis_dir = Path(args.analysis_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    h1_path = analysis_dir / "STAR_M1_DELTA_CR_MISSING_WITH_PARTNER_COUNTS.tsv"
    h0_path = analysis_dir / "STAR_EXACT_VS_CR.tsv"

    if not h1_path.exists():
        raise FileNotFoundError(f"Missing file: {h1_path}")
    if not h0_path.exists():
        raise FileNotFoundError(f"Missing file: {h0_path}")

    h1_rows = read_tsv(h1_path)
    h0_all = read_tsv(h0_path)
    h0_rows = [r for r in h0_all if int(r["cr_raw_count"]) == 0]

    combined = []

    for r in h1_rows:
        combined.append(
            {
                "hamming_mode": "1",
                "feature": r["feature"],
                "barcode_tru": r["barcode_tru"],
                "barcode_nxt": r["barcode_nxt"],
                "star_metric_type": "star_m1_minus_m0",
                "star_metric": str(int(r["delta"])),
                "star_m0": r["star_m0"],
                "star_m1": r["star_m1"],
                "cr_raw_count": r["cr_self"],
                "in_cr_call_list": "0",
                "classification": r["classification"],
            }
        )

    for r in h0_rows:
        combined.append(
            {
                "hamming_mode": "0",
                "feature": r["feature"],
                "barcode_tru": r["barcode_tru"],
                "barcode_nxt": r["barcode_nxt"],
                "star_metric_type": "star_exact_count",
                "star_metric": str(int(r["star_exact_count"])),
                "star_m0": r["star_exact_count"],
                "star_m1": "",
                "cr_raw_count": r["cr_raw_count"],
                "in_cr_call_list": r["in_cr_call_list"],
                "classification": "exact_missing",
            }
        )

    combined.sort(key=lambda x: (x["hamming_mode"], -int(x["star_metric"]), x["feature"], x["barcode_tru"]))

    combined_out = out_dir / "UCSF_STAR_vs_CR_MISSES_H0_H1.tsv"
    with open(combined_out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "hamming_mode",
                "feature",
                "barcode_tru",
                "barcode_nxt",
                "star_metric_type",
                "star_metric",
                "star_m0",
                "star_m1",
                "cr_raw_count",
                "in_cr_call_list",
                "classification",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(combined)

    # Summary table for manuscript text
    summary_rows = []

    # H0 summary
    h0_count = len(h0_rows)
    h0_sum = sum(int(r["star_exact_count"]) for r in h0_rows)
    summary_rows.append({"hamming_mode": "0", "metric": "miss_pairs", "value": str(h0_count)})
    summary_rows.append({"hamming_mode": "0", "metric": "miss_star_sum", "value": str(h0_sum)})

    # H1 summary
    h1_count = len(h1_rows)
    h1_sum = sum(int(r["delta"]) for r in h1_rows)
    summary_rows.append({"hamming_mode": "1", "metric": "miss_pairs", "value": str(h1_count)})
    summary_rows.append({"hamming_mode": "1", "metric": "miss_star_sum", "value": str(h1_sum)})

    h1_class_count = Counter(r["classification"] for r in h1_rows)
    h1_class_sum: Dict[str, int] = defaultdict(int)
    for r in h1_rows:
        h1_class_sum[r["classification"]] += int(r["delta"])

    for klass in sorted(h1_class_count):
        summary_rows.append(
            {
                "hamming_mode": "1",
                "metric": f"{klass}_pairs",
                "value": str(h1_class_count[klass]),
            }
        )
        summary_rows.append(
            {
                "hamming_mode": "1",
                "metric": f"{klass}_star_sum",
                "value": str(h1_class_sum[klass]),
            }
        )

    summary_out = out_dir / "UCSF_STAR_vs_CR_MISSES_H0_H1_SUMMARY.tsv"
    with open(summary_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["hamming_mode", "metric", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    # Markdown preview
    md_out = out_dir / "UCSF_STAR_vs_CR_MISSES_H0_H1_PREVIEW.md"
    top_n = args.top_n
    h0_top = sorted(h0_rows, key=lambda r: int(r["star_exact_count"]), reverse=True)[:top_n]
    h1_top = sorted(h1_rows, key=lambda r: int(r["delta"]), reverse=True)[:top_n]

    with open(md_out, "w") as fh:
        fh.write("# UCSF STAR vs CR Missing Rows (H0 and H1)\n\n")
        fh.write("## H0 (Exact-match misses; STAR m=0, CR raw=0)\n\n")
        fh.write("| feature | barcode_tru | star_exact_count | in_cr_call_list |\n")
        fh.write("|---|---|---:|---:|\n")
        for r in h0_top:
            fh.write(
                f"| {r['feature']} | {r['barcode_tru']} | {int(r['star_exact_count'])} | {int(r['in_cr_call_list'])} |\n"
            )

        fh.write("\n## H1 (1-Hamming rescue misses; STAR delta, CR raw=0)\n\n")
        fh.write("| feature | barcode_tru | star_m1_minus_m0 | classification |\n")
        fh.write("|---|---|---:|---|\n")
        for r in h1_top:
            fh.write(
                f"| {r['feature']} | {r['barcode_tru']} | {int(r['delta'])} | {r['classification']} |\n"
            )

    print(combined_out)
    print(summary_out)
    print(md_out)


if __name__ == "__main__":
    main()
