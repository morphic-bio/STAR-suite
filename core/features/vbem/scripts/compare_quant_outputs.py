#!/usr/bin/env python3
"""
compare_quant_outputs.py - Compare two quant.sf-style files (TPM/NumReads)
and optionally generate a volcano-style plot.

Usage:
  python3 compare_quant_outputs.py --a path/to/quant.sf --b path/to/quant.sf \
    --out-prefix /path/to/output_prefix

Notes:
  - Spearman/Pearson computed for TPM and NumReads.
  - Volcano plot uses log2 fold-change on the chosen metric (default TPM) and
    -log10 p-values from a two-proportion z-test on NumReads.
"""

import argparse
import csv
import math
import sys
from typing import Dict, Tuple

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - environment-specific
    print("ERROR: numpy is required for this script.", file=sys.stderr)
    raise

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:  # pragma: no cover - environment-specific
    HAVE_MPL = False

try:
    from scipy.stats import spearmanr
    from scipy.stats import norm
    HAVE_SCIPY = True
except Exception:  # pragma: no cover - environment-specific
    HAVE_SCIPY = False


def load_quant(path: str) -> Dict[str, Tuple[float, float]]:
    data = {}
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Name", "TPM", "NumReads"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing required columns in {path}: {required}")
        for row in reader:
            name = row["Name"]
            data[name] = (float(row["TPM"]), float(row["NumReads"]))
    return data


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    xm = x.mean()
    ym = y.mean()
    cov = np.sum((x - xm) * (y - ym))
    denom = math.sqrt(np.sum((x - xm) ** 2) * np.sum((y - ym) ** 2))
    if denom == 0:
        return float("nan")
    return float(cov / denom)


def rankdata(a: np.ndarray) -> np.ndarray:
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    sorted_a = a[order]
    n = len(a)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_a[j + 1] == sorted_a[i]:
            j += 1
        # Average rank for ties (1-based ranks)
        avg_rank = (i + j) / 2.0 + 1.0
        ranks[order[i:j + 1]] = avg_rank
        i = j + 1
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if HAVE_SCIPY:
        return float(spearmanr(x, y).correlation)
    return pearson(rankdata(x), rankdata(y))


def two_prop_pvals(counts_a: np.ndarray, counts_b: np.ndarray) -> np.ndarray:
    total_a = counts_a.sum()
    total_b = counts_b.sum()
    if total_a == 0 or total_b == 0:
        return np.ones_like(counts_a)
    p1 = counts_a / total_a
    p2 = counts_b / total_b
    p = (counts_a + counts_b) / (total_a + total_b)
    se = np.sqrt(p * (1 - p) * (1 / total_a + 1 / total_b))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (p2 - p1) / se
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    if HAVE_SCIPY:
        pvals = 2.0 * norm.sf(np.abs(z))
    else:
        # Use erfc approximation without SciPy.
        pvals = np.vectorize(lambda v: math.erfc(abs(v) / math.sqrt(2.0)))(z)
    pvals = np.clip(pvals, 1e-300, 1.0)
    return pvals


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare two quant.sf-style files")
    parser.add_argument("--a", required=True, help="First quant.sf file")
    parser.add_argument("--b", required=True, help="Second quant.sf file")
    parser.add_argument("--label-a", default="A", help="Label for first file")
    parser.add_argument("--label-b", default="B", help="Label for second file")
    parser.add_argument("--out-prefix", required=True, help="Output prefix for plots")
    parser.add_argument("--metric", choices=["TPM", "NumReads"], default="TPM",
                        help="Metric for log2 fold-change in volcano plot")
    parser.add_argument("--eps", type=float, default=1e-3,
                        help="Pseudocount for log2 fold-change")
    parser.add_argument("--no-plot", action="store_true",
                        help="Skip plotting (print correlations only)")
    args = parser.parse_args()

    data_a = load_quant(args.a)
    data_b = load_quant(args.b)
    common = sorted(set(data_a) & set(data_b))
    if not common:
        print("ERROR: No common entries between files.", file=sys.stderr)
        return 2

    tpm_a = np.array([data_a[k][0] for k in common], dtype=float)
    tpm_b = np.array([data_b[k][0] for k in common], dtype=float)
    nr_a = np.array([data_a[k][1] for k in common], dtype=float)
    nr_b = np.array([data_b[k][1] for k in common], dtype=float)

    print("=== Correlations ===")
    print(f"Common entries: {len(common)}")
    print(f"TPM Pearson:   {pearson(tpm_a, tpm_b):.6f}")
    print(f"TPM Spearman:  {spearman(tpm_a, tpm_b):.6f}")
    print(f"NumReads Pearson:  {pearson(nr_a, nr_b):.6f}")
    print(f"NumReads Spearman: {spearman(nr_a, nr_b):.6f}")

    if args.no_plot:
        return 0

    if not HAVE_MPL:
        print("ERROR: matplotlib is required for plotting.", file=sys.stderr)
        return 2

    # Volcano plot: log2 fold-change vs -log10 p-value (two-proportion z-test).
    pvals = two_prop_pvals(nr_a, nr_b)
    neg_log_p = -np.log10(pvals)

    if args.metric == "TPM":
        x = np.log2((tpm_b + args.eps) / (tpm_a + args.eps))
        metric_label = "TPM"
    else:
        x = np.log2((nr_b + args.eps) / (nr_a + args.eps))
        metric_label = "NumReads"

    out_png = f"{args.out_prefix}.volcano.{args.metric.lower()}.png"

    plt.figure(figsize=(8, 6))
    plt.scatter(x, neg_log_p, s=4, alpha=0.3, linewidths=0)
    plt.xlabel(f"log2 fold-change ({args.label_b}/{args.label_a}) [{metric_label}]")
    plt.ylabel("-log10 p-value (two-proportion z-test on NumReads)")
    plt.title(f"Volcano plot ({metric_label})")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"Volcano plot written: {out_png}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
