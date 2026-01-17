#!/usr/bin/env python3
"""
Plot SNP mismatch-fraction histogram + Kneedle curve from a BAM file.

This mirrors the auto-estimation logic in SlamQuant::estimateSnpMismatchFrac.
"""
import argparse
from collections import defaultdict
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def iter_reads(bam):
    try:
        return bam.fetch()
    except ValueError:
        return bam.fetch(until_eof=True)


def collect_mismatch_fractions(bam_path, max_reads, min_cov):
    import pysam
    bam = pysam.AlignmentFile(bam_path, 'rb')

    cov = defaultdict(int)
    mis = defaultdict(int)
    reads = 0

    for read in iter_reads(bam):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        reads += 1
        if max_reads > 0 and reads > max_reads:
            break
        seq = read.query_sequence
        if not seq:
            continue

        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
            if read_pos is None or ref_pos is None or ref_base is None:
                continue
            key = (read.reference_id, ref_pos)
            cov[key] += 1
            if seq[read_pos].upper() != ref_base.upper():
                mis[key] += 1

    bam.close()

    fractions = []
    for key, c in cov.items():
        if c >= min_cov:
            m = mis.get(key, 0)
            fractions.append(m / c)

    return fractions, reads


def kneedle_threshold(fractions, bins, epsilon):
    hist = np.zeros(bins, dtype=np.uint64)
    for f in fractions:
        f = min(max(f, 0.0), 1.0)
        idx = int(f * (bins - 1))
        hist[idx] += 1

    cumulative = np.cumsum(hist[::-1])[::-1].astype(float)
    log_cum = np.log1p(cumulative)

    y_min = log_cum[-1]
    y_max = log_cum[0]
    if y_max <= y_min:
        return 0.0, 0, hist, cumulative, log_cum, None

    x_norm = np.linspace(0.0, 1.0, bins)
    y_norm = (log_cum - y_min) / (y_max - y_min)
    dist = y_norm + x_norm - 1.0
    knee_bin = int(np.argmax(dist))
    max_dist = dist[knee_bin]
    if max_dist < epsilon:
        return 0.0, knee_bin, hist, cumulative, log_cum, dist

    threshold = knee_bin / (bins - 1)
    return threshold, knee_bin, hist, cumulative, log_cum, dist


def main():
    parser = argparse.ArgumentParser(description="Plot SNP mismatch-fraction curve + Kneedle fit")
    parser.add_argument('--bam', required=True, help='BAM file path')
    parser.add_argument('--max-reads', type=int, default=100000, help='Max reads to scan')
    parser.add_argument('--min-coverage', type=int, default=10, help='Min coverage per site')
    parser.add_argument('--bins', type=int, default=100, help='Histogram bins')
    parser.add_argument('--epsilon', type=float, default=0.02, help='Kneedle epsilon')
    parser.add_argument('--output', default='docs/figures/snp_threshold_curve.png', help='Output PNG path')
    parser.add_argument('--stats-out', default='docs/figures/snp_threshold_curve.tsv', help='Output TSV path')
    args = parser.parse_args()

    fractions, reads = collect_mismatch_fractions(args.bam, args.max_reads, args.min_coverage)
    if not fractions:
        print("No eligible positions found; cannot plot.")
        return

    threshold, knee_bin, hist, cumulative, log_cum, dist = kneedle_threshold(
        fractions, args.bins, args.epsilon)

    print(f"Scanned reads: {reads}")
    print(f"Eligible sites: {len(fractions)} (min_cov={args.min_coverage})")
    print(f"Kneedle threshold: {threshold:.4f} (bin {knee_bin})")

    # Save TSV
    bin_centers = np.linspace(0.0, 1.0, args.bins)
    with open(args.stats_out, 'w') as out:
        out.write("\t".join([
            "bin", "bin_center", "count", "cumulative", "log1p_cum", "dist"
        ]) + "\n")
        for i in range(args.bins):
            out.write("\t".join(map(str, [
                i, bin_centers[i], int(hist[i]), float(cumulative[i]),
                float(log_cum[i]), float(dist[i]) if dist is not None else float('nan')
            ])) + "\n")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))

    # Histogram
    axes[0].bar(bin_centers, hist, width=1.0/args.bins, color='steelblue', alpha=0.8)
    axes[0].set_title('Mismatch Fraction Histogram')
    axes[0].set_xlabel('Mismatch Fraction')
    axes[0].set_ylabel('Count')
    axes[0].grid(True, alpha=0.2)

    # Cumulative + log
    axes[1].plot(bin_centers, cumulative, 'k-', label='Cumulative N(>=f)')
    axes[1].plot(bin_centers, log_cum, 'r--', label='log1p cumulative')
    if threshold > 0:
        axes[1].axvline(threshold, color='purple', linestyle='--', label=f'knee={threshold:.3f}')
    axes[1].set_title('Cumulative Curve')
    axes[1].set_xlabel('Mismatch Fraction')
    axes[1].set_ylabel('Count / log1p')
    axes[1].grid(True, alpha=0.2)
    axes[1].legend()

    # Kneedle distance
    if dist is not None:
        axes[2].plot(bin_centers, dist, 'g-', label='dist = y + x - 1')
        axes[2].axhline(0, color='gray', linestyle='--', alpha=0.5)
        if threshold > 0:
            axes[2].axvline(threshold, color='purple', linestyle='--', label=f'knee={threshold:.3f}')
        axes[2].set_title('Kneedle Distance')
        axes[2].set_xlabel('Mismatch Fraction')
        axes[2].set_ylabel('Distance')
        axes[2].grid(True, alpha=0.2)
        axes[2].legend()
    else:
        axes[2].axis('off')

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Plot saved to: {args.output}")
    print(f"Stats saved to: {args.stats_out}")


if __name__ == '__main__':
    main()
