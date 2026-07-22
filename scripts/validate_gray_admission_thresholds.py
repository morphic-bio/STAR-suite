#!/usr/bin/env python3
"""
5-fold validation for a simple gray-zone admission rule.

The rule is intentionally low-flexibility:
- gray zone is fixed upstream
- admission uses only monotone thresholds on
  genes_detected, entropy, top_gene_frac
- thresholds are derived from training-fold class medians

Two evaluation modes are supported:
- full_gray: predict keep/drop across the whole gray zone
- rescue_only: preserve current accepted gray-zone barcodes and only decide
  whether to rescue current rejected gray-zone barcodes
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import random
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


def midpoint(a: float, b: float) -> float:
    return (a + b) / 2.0


def safe_div(num: float, den: float) -> float:
    return num / den if den else float("nan")


def fmt(x: float) -> str:
    if math.isnan(x):
        return "NA"
    return f"{x:.6f}"


def load_rows(path: Path) -> list[dict[str, object]]:
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "barcode": row["barcode"],
                    "label": int(row["in_cr"]),
                    "genes_detected": float(row["genes_detected"]),
                    "top_gene_frac": float(row["top_gene_frac"]),
                    "entropy": float(row["entropy"]),
                }
            )
    return rows


def derive_thresholds(train_rows: list[dict[str, object]]) -> dict[str, float]:
    pos = [r for r in train_rows if r["label"] == 1]
    neg = [r for r in train_rows if r["label"] == 0]
    if not pos or not neg:
        raise SystemExit("Training fold lacks positive or negative examples")

    pos_genes = [float(r["genes_detected"]) for r in pos]
    neg_genes = [float(r["genes_detected"]) for r in neg]
    pos_entropy = [float(r["entropy"]) for r in pos]
    neg_entropy = [float(r["entropy"]) for r in neg]
    pos_topfrac = [float(r["top_gene_frac"]) for r in pos]
    neg_topfrac = [float(r["top_gene_frac"]) for r in neg]

    return {
        "genes_detected_min": midpoint(statistics.median(pos_genes), statistics.median(neg_genes)),
        "entropy_min": midpoint(statistics.median(pos_entropy), statistics.median(neg_entropy)),
        "top_gene_frac_max": midpoint(statistics.median(pos_topfrac), statistics.median(neg_topfrac)),
    }


def predict(row: dict[str, object], thr: dict[str, float]) -> int:
    return int(
        float(row["genes_detected"]) >= thr["genes_detected_min"]
        and float(row["entropy"]) >= thr["entropy_min"]
        and float(row["top_gene_frac"]) <= thr["top_gene_frac_max"]
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--per-barcode-metrics", type=Path, required=True, help="per_barcode_metrics.tsv from analyze_gray_admission_heuristics.py")
    ap.add_argument("--gray-zone-barcodes", type=Path, required=True, help="gray_zone_barcodes.tsv from gray admission harness")
    ap.add_argument("--base-filtered-barcodes", type=Path, required=True, help="Base/current filtered barcodes.tsv")
    ap.add_argument("--cr-filtered-barcodes", type=Path, required=True, help="CR filtered barcodes.tsv(.gz)")
    ap.add_argument(
        "--mode",
        choices=["full_gray", "rescue_only"],
        default="rescue_only",
        help="Whether to classify the whole gray zone or only rescue among current rejected gray-zone barcodes",
    )
    ap.add_argument("--folds", type=int, default=5, help="Number of CV folds")
    ap.add_argument("--seed", type=int, default=42, help="Shuffle seed")
    ap.add_argument("--outdir", type=Path, required=True, help="Output directory")
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(args.per_barcode_metrics)
    gray_zone = read_barcodes(args.gray_zone_barcodes)
    base_filtered = read_barcodes(args.base_filtered_barcodes)
    cr_filtered = read_barcodes(args.cr_filtered_barcodes)
    base_gray = gray_zone & base_filtered

    row_barcodes = {str(r["barcode"]) for r in rows}
    if row_barcodes != gray_zone:
        missing = len(gray_zone - row_barcodes)
        extra = len(row_barcodes - gray_zone)
        raise SystemExit(f"Gray-zone mismatch between metrics and gray list: missing={missing} extra={extra}")

    if args.mode == "rescue_only":
        rows = [r for r in rows if str(r["barcode"]) not in base_gray]
        if not rows:
            raise SystemExit("No current-rejected gray-zone barcodes remain after applying rescue_only mode")

    rng = random.Random(args.seed)
    rows_shuffled = rows[:]
    rng.shuffle(rows_shuffled)
    folds: list[list[dict[str, object]]] = [[] for _ in range(args.folds)]
    for i, row in enumerate(rows_shuffled):
        folds[i % args.folds].append(row)

    oof_pred: dict[str, int] = {}
    fold_summary_path = args.outdir / "fold_summary.tsv"
    with open(fold_summary_path, "w", encoding="utf-8") as out:
        out.write(
            "fold\tn_train\tn_test\tgenes_detected_min\tentropy_min\ttop_gene_frac_max\t"
            "tp\tfp\ttn\tfn\tprecision\trecall\tspecificity\tbalanced_accuracy\n"
        )
        for i in range(args.folds):
            test_rows = folds[i]
            train_rows = [r for j, fold in enumerate(folds) if j != i for r in fold]
            thr = derive_thresholds(train_rows)

            tp = fp = tn = fn = 0
            for row in test_rows:
                y = int(row["label"])
                yhat = predict(row, thr)
                oof_pred[str(row["barcode"])] = yhat
                if yhat == 1 and y == 1:
                    tp += 1
                elif yhat == 1 and y == 0:
                    fp += 1
                elif yhat == 0 and y == 0:
                    tn += 1
                else:
                    fn += 1

            precision = safe_div(tp, tp + fp)
            recall = safe_div(tp, tp + fn)
            specificity = safe_div(tn, tn + fp)
            bal_acc = safe_div((recall if not math.isnan(recall) else 0.0) + (specificity if not math.isnan(specificity) else 0.0), 2.0)
            out.write(
                f"{i+1}\t{len(train_rows)}\t{len(test_rows)}\t{thr['genes_detected_min']:.6f}\t"
                f"{thr['entropy_min']:.6f}\t{thr['top_gene_frac_max']:.6f}\t"
                f"{tp}\t{fp}\t{tn}\t{fn}\t{fmt(precision)}\t{fmt(recall)}\t{fmt(specificity)}\t{fmt(bal_acc)}\n"
            )

    tp = fp = tn = fn = 0
    predicted_positive = set()
    for row in rows:
        bc = str(row["barcode"])
        y = int(row["label"])
        yhat = oof_pred[bc]
        if yhat == 1:
            predicted_positive.add(bc)
        if yhat == 1 and y == 1:
            tp += 1
        elif yhat == 1 and y == 0:
            fp += 1
        elif yhat == 0 and y == 0:
            tn += 1
        else:
            fn += 1

    precision = safe_div(tp, tp + fp)
    recall = safe_div(tp, tp + fn)
    specificity = safe_div(tn, tn + fp)
    bal_acc = safe_div((recall if not math.isnan(recall) else 0.0) + (specificity if not math.isnan(specificity) else 0.0), 2.0)
    f1 = safe_div(2 * tp, 2 * tp + fp + fn)

    base_outside_gray = base_filtered - gray_zone
    if args.mode == "rescue_only":
        cv_filtered = base_outside_gray | base_gray | predicted_positive
    else:
        cv_filtered = base_outside_gray | predicted_positive
    base_common = len(base_filtered & cr_filtered)
    base_union = len(base_filtered | cr_filtered)
    cv_common = len(cv_filtered & cr_filtered)
    cv_union = len(cv_filtered | cr_filtered)

    summary_path = args.outdir / "summary.tsv"
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write("metric\tvalue\n")
        out.write(f"mode\t{args.mode}\n")
        out.write(f"folds\t{args.folds}\n")
        out.write(f"seed\t{args.seed}\n")
        out.write(f"gray_zone_barcodes\t{len(gray_zone)}\n")
        out.write(f"base_gray_barcodes\t{len(base_gray)}\n")
        out.write(f"eval_rows\t{len(rows)}\n")
        out.write(f"tp\t{tp}\n")
        out.write(f"fp\t{fp}\n")
        out.write(f"tn\t{tn}\n")
        out.write(f"fn\t{fn}\n")
        out.write(f"precision\t{fmt(precision)}\n")
        out.write(f"recall\t{fmt(recall)}\n")
        out.write(f"specificity\t{fmt(specificity)}\n")
        out.write(f"balanced_accuracy\t{fmt(bal_acc)}\n")
        out.write(f"f1\t{fmt(f1)}\n")
        out.write(f"predicted_positive_gray\t{len(predicted_positive)}\n")
        out.write(f"predicted_positive_gray_in_cr\t{len(predicted_positive & cr_filtered)}\n")
        out.write(f"predicted_positive_gray_not_in_cr\t{len(predicted_positive - cr_filtered)}\n")
        out.write(f"base_filtered_cells\t{len(base_filtered)}\n")
        out.write(f"base_common_cr\t{base_common}\n")
        out.write(f"base_union_cr\t{base_union}\n")
        out.write(f"base_jaccard_cr\t{fmt(safe_div(base_common, base_union))}\n")
        out.write(f"cv_filtered_cells\t{len(cv_filtered)}\n")
        out.write(f"cv_common_cr\t{cv_common}\n")
        out.write(f"cv_union_cr\t{cv_union}\n")
        out.write(f"cv_jaccard_cr\t{fmt(safe_div(cv_common, cv_union))}\n")
        out.write(f"cv_added_vs_base\t{len(cv_filtered - base_filtered)}\n")
        out.write(f"cv_removed_vs_base\t{len(base_filtered - cv_filtered)}\n")
        out.write(f"cv_added_in_cr\t{len((cv_filtered - base_filtered) & cr_filtered)}\n")
        out.write(f"cv_removed_in_cr\t{len((base_filtered - cv_filtered) & cr_filtered)}\n")

    for name, barcode_set in {
        "predicted_positive_gray.tsv": predicted_positive,
        "predicted_negative_gray.tsv": gray_zone - predicted_positive,
        "tp_barcodes.tsv": {bc for bc in predicted_positive if bc in cr_filtered},
        "fp_barcodes.tsv": {bc for bc in predicted_positive if bc not in cr_filtered},
        "fn_barcodes.tsv": {str(r["barcode"]) for r in rows if int(r["label"]) == 1 and oof_pred[str(r["barcode"])] == 0},
        "cv_filtered_barcodes.tsv": cv_filtered,
    }.items():
        with open(args.outdir / name, "w", encoding="utf-8") as out:
            for bc in sorted(barcode_set):
                out.write(f"{bc}\n")

    print(summary_path.read_text(), end="")
    print(f"Wrote fold summary to: {fold_summary_path}")


if __name__ == "__main__":
    main()
