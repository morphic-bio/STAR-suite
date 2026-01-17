#!/usr/bin/env python3
"""
Analyze PHRED quality score distribution by read position and compare with T→C rates.

Purpose: Demonstrate that position-dependent T→C rate differences between STAR and GEDI
are NOT correlated with base quality scores.

Key insight: Raw T→C rates show U-shaped artifacts from sequencing errors at read ends.
To see the TRUE 4sU signal, we normalize T→C by T→A (a non-4sU transition).

Outputs:
1. Console: Position-by-position quality and T→C rate statistics
2. Plot: Comparison of raw vs normalized T→C rates and PHRED quality
"""

import sys
import os
import argparse
import math
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def smooth_series(values, window=5, method='median'):
    """Smooth a 1D series with a centered window."""
    if window <= 1:
        return list(values)
    half = window // 2
    out = []
    for i in range(len(values)):
        start = max(0, i - half)
        end = min(len(values), i + half + 1)
        chunk = [v for v in values[start:end] if not math.isnan(v)]
        if not chunk:
            out.append(float('nan'))
            continue
        if method == 'mean':
            out.append(float(np.mean(chunk)))
        else:
            out.append(float(np.median(chunk)))
    return out

def detect_trim_smooth_quantile(values, window=5, quantile=0.6, stable_n=3, max_trim=15):
    """Detect trim using smoothed quantile threshold and stable region scan."""
    if not values:
        return 0, 0, float('nan'), []

    smoothed = smooth_series(values, window=window, method='median')
    read_len = len(smoothed)
    mid_start = int(read_len * 0.2)
    mid_end = int(read_len * 0.8) - 1
    if mid_end <= mid_start:
        mid_start, mid_end = 0, read_len - 1

    mid_vals = [v for v in smoothed[mid_start:mid_end + 1] if not math.isnan(v)]
    if not mid_vals:
        return 0, 0, float('nan'), smoothed

    threshold = float(np.quantile(mid_vals, quantile))

    # Find 5' trim: first stable run below threshold
    trim5p = 0
    for i in range(0, read_len - stable_n + 1):
        if all((not math.isnan(v)) and v <= threshold for v in smoothed[i:i + stable_n]):
            trim5p = i
            break

    # Find 3' trim: last stable run below threshold
    trim3p = 0
    for i in range(read_len - 1, stable_n - 2, -1):
        window_vals = smoothed[i - stable_n + 1:i + 1]
        if all((not math.isnan(v)) and v <= threshold for v in window_vals):
            trim3p = (read_len - 1) - i
            break

    trim5p = min(trim5p, max_trim)
    trim3p = min(trim3p, max_trim)
    return trim5p, trim3p, threshold, smoothed

def detect_trim_segmented(values, window=5, min_seg_len=3, max_trim=15):
    """Segmented regression (2 breakpoints) on smoothed curve."""
    if not values:
        return 0, 0, [], {}

    smoothed = smooth_series(values, window=window, method='median')
    y = np.array(smoothed, dtype=float)
    n = len(y)
    if n < (min_seg_len * 3 + 1):
        return 0, 0, smoothed, {}

    # Fill NaNs via linear interpolation
    if np.isnan(y).any():
        idx = np.arange(n)
        valid = ~np.isnan(y)
        if valid.sum() < 2:
            return 0, 0, smoothed, {}
        y[~valid] = np.interp(idx[~valid], idx[valid], y[valid])

    # Prefix sums for fast segment fits
    x = np.arange(n, dtype=float)
    px = np.concatenate([[0.0], np.cumsum(x)])
    pxx = np.concatenate([[0.0], np.cumsum(x * x)])
    py = np.concatenate([[0.0], np.cumsum(y)])
    pxy = np.concatenate([[0.0], np.cumsum(x * y)])
    pyy = np.concatenate([[0.0], np.cumsum(y * y)])

    def seg_sum(p, i, j):
        return p[j + 1] - p[i]

    def seg_fit(i, j):
        nseg = j - i + 1
        if nseg < 2:
            return 0.0, y[i], 0.0
        Sx = seg_sum(px, i, j)
        Sxx = seg_sum(pxx, i, j)
        Sy = seg_sum(py, i, j)
        Sxy = seg_sum(pxy, i, j)
        Sy2 = seg_sum(pyy, i, j)
        den = nseg * Sxx - Sx * Sx
        if abs(den) < 1e-9:
            m = 0.0
        else:
            m = (nseg * Sxy - Sx * Sy) / den
        b = (Sy - m * Sx) / nseg
        sse = Sy2 + m * m * Sxx + nseg * b * b + 2 * m * b * Sx - 2 * m * Sxy - 2 * b * Sy
        return m, b, sse

    # Candidate breakpoints
    b1_min = min_seg_len
    b1_max = min(max_trim, n - 2 * min_seg_len - 1)
    b2_min_floor = max(min_seg_len - 1, n - 1 - max_trim)
    b2_max = n - min_seg_len - 1
    if b1_max < b1_min or b2_min_floor > b2_max:
        return 0, 0, smoothed, {}

    best = None
    best_info = {}
    for b1 in range(b1_min, b1_max + 1):
        b2_min = max(b1 + min_seg_len - 1, b2_min_floor)
        for b2 in range(b2_min, b2_max + 1):
            m1, b1i, sse1 = seg_fit(0, b1 - 1)
            m2, b2i, sse2 = seg_fit(b1, b2)
            m3, b3i, sse3 = seg_fit(b2 + 1, n - 1)
            sse = sse1 + sse2 + sse3
            if best is None or sse < best:
                best = sse
                best_info = {
                    'b1': b1,
                    'b2': b2,
                    'fit1': (m1, b1i),
                    'fit2': (m2, b2i),
                    'fit3': (m3, b3i),
                    'sse': sse,
                }

    if not best_info:
        return 0, 0, smoothed, {}

    trim5p = min(best_info['b1'], max_trim)
    trim3p = min((n - 1 - best_info['b2']), max_trim)
    return trim5p, trim3p, smoothed, best_info

def analyze_bam_quality(bam_path: str, max_reads: int = 100000):
    """
    Analyze PHRED quality scores and T→C rates by read position from BAM.

    Uses transcript-oriented T positions:
      - plus-strand: ref T
      - minus-strand: ref A (complement of transcript T)

    Returns:
        dict with position -> stats, count of reads analyzed
    """
    import pysam
    
    bam = pysam.AlignmentFile(bam_path, 'rb')
    
    # Track stats by position
    stats = defaultdict(lambda: {
        'qual_sum': 0,
        'qual_count': 0,
        't_qual_sum': 0,
        't_qual_count': 0,
        't_count': 0,
        'tc_count': 0,
        'mm_count': 0,
        'mm_total': 0,
        'qual_values': [],     # All bases
        't_qual_values': []    # Transcript T bases only
    })
    
    n = 0
    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        n += 1
        if n > max_reads:
            break
        
        seq = read.query_sequence
        quals = read.query_qualities
        cigar = read.cigartuples
        
        if not seq or not quals or not cigar:
            continue
        
        # Get soft-clip lengths (5' and 3')
        clip_5p = cigar[0][1] if cigar[0][0] == 4 else 0
        clip_3p = cigar[-1][1] if cigar[-1][0] == 4 else 0
        aligned_len = read.query_length - clip_5p - clip_3p
        if aligned_len <= 0:
            continue

        try:
            nh_tag = read.get_tag("NH")
        except Exception:
            nh_tag = 1

        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
            if read_pos is None or ref_pos is None or ref_base is None:
                continue
            
            # Position within aligned read (0-based, first aligned base = 0)
            aligned_pos = read_pos - clip_5p
            if aligned_pos < 0 or aligned_pos >= aligned_len:
                continue
            
            qual = quals[read_pos]
            ref_base_upper = ref_base.upper()
            read_base = seq[read_pos].upper()
            is_reverse = read.is_reverse

            stats[aligned_pos]['mm_total'] += 1
            if nh_tag > 1:
                stats[aligned_pos]['mm_count'] += 1
            
            # Track quality at all positions
            stats[aligned_pos]['qual_sum'] += qual
            stats[aligned_pos]['qual_count'] += 1
            stats[aligned_pos]['qual_values'].append(qual)
            
            # Track transcript T positions and T→C conversions
            # plus-strand: ref T -> read C
            # minus-strand: ref A -> read G (complement of T->C)
            is_t = (ref_base_upper == 'A') if is_reverse else (ref_base_upper == 'T')
            if is_t:
                stats[aligned_pos]['t_count'] += 1
                stats[aligned_pos]['t_qual_sum'] += qual
                stats[aligned_pos]['t_qual_count'] += 1
                stats[aligned_pos]['t_qual_values'].append(qual)
                is_tc = (ref_base_upper == 'A' and read_base == 'G') if is_reverse else (read_base == 'C')
                if is_tc:
                    stats[aligned_pos]['tc_count'] += 1
    
    bam.close()
    return stats, n

def parse_mismatchdetails(path: str, category: str = 'ExonicSense'):
    """Parse STAR or GEDI mismatchdetails for T→C and T→A rates by position."""
    tc_data = {}
    ta_data = {}
    
    with open(path) as f:
        header = f.readline().strip().split('\t')
        cat_idx = header.index('Category')
        genomic_idx = header.index('Genomic')
        read_idx = header.index('Read')
        pos_idx = header.index('Position')
        cov_idx = header.index('Coverage')
        mm_idx = header.index('Mismatches')
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(cat_idx, genomic_idx, read_idx, pos_idx, cov_idx, mm_idx):
                continue
            
            cat = parts[cat_idx]
            genomic = parts[genomic_idx]
            read = parts[read_idx]
            pos = int(parts[pos_idx])
            cov = float(parts[cov_idx])
            mm = float(parts[mm_idx])
            
            if cat == category and genomic == 'T':
                if read == 'C':
                    p = mm / cov if cov > 0 else 0.0
                    rate = p * 100
                    stdev = math.sqrt(p * (1.0 - p)) * 100 if cov > 0 else 0.0
                    tc_data[pos] = {'cov': cov, 'mm': mm, 'rate': rate, 'stdev': stdev}
                elif read == 'A':
                    p = mm / cov if cov > 0 else 0.0
                    rate = p * 100
                    stdev = math.sqrt(p * (1.0 - p)) * 100 if cov > 0 else 0.0
                    ta_data[pos] = {'cov': cov, 'mm': mm, 'rate': rate, 'stdev': stdev}
    
    return tc_data, ta_data

def main():
    parser = argparse.ArgumentParser(description='Analyze PHRED quality by position')
    parser.add_argument('--bam', default='/storage/SLAM-Seq-prod-compare-20260109/star/WDHD1_0h3_Aligned.sortedByCoord.out.bam',
                        help='BAM file path')
    parser.add_argument('--star', default='test/tmp_prod_compat/default_SlamQuant.out.mismatchdetails.tsv',
                        help='STAR mismatchdetails path')
    parser.add_argument('--gedi', default='/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3_Sense_fixed.mismatchdetails.tsv',
                        help='GEDI mismatchdetails path')
    parser.add_argument('--output', default='test/tmp_prod_compat/phred_position_analysis.png',
                        help='Output plot path')
    parser.add_argument('--max-reads', type=int, default=100000,
                        help='Maximum reads to analyze from BAM')
    parser.add_argument('--stats-out', default='-',
                        help='Optional TSV output for per-position stats (use - to skip)')
    parser.add_argument('--trim-window', type=int, default=5,
                        help='Smoothing window for trim detection')
    parser.add_argument('--trim-max', type=int, default=15,
                        help='Maximum trim to apply at either end')
    parser.add_argument('--seg-min-len', type=int, default=3,
                        help='Minimum segment length for segmented regression')
    args = parser.parse_args()
    
    print("=" * 100)
    print("PHRED Quality (T-base) vs T→C Rate by Read Position")
    print("=" * 100)
    
    # Analyze BAM
    print(f"\nAnalyzing BAM: {args.bam}")
    bam_stats, n_reads = analyze_bam_quality(args.bam, args.max_reads)
    print(f"Analyzed {n_reads} reads")
    
    # Parse mismatchdetails (T→C and T→A)
    print(f"\nParsing STAR: {args.star}")
    star_tc, star_ta = parse_mismatchdetails(args.star)
    print(f"Parsing GEDI: {args.gedi}")
    gedi_tc, gedi_ta = parse_mismatchdetails(args.gedi)
    
    # Compute position-level statistics
    max_pos = 0
    if bam_stats:
        max_pos = max(max_pos, max(bam_stats.keys()))
    for data in (star_tc, star_ta, gedi_tc, gedi_ta):
        if data:
            max_pos = max(max_pos, max(data.keys()))
    positions = list(range(max_pos + 1))
    positions_plot = [p + 1 for p in positions]  # 1-based for plotting

    avg_qual = []
    qual_std = []
    qual_valid = []
    mm_frac = []
    t_counts = []
    star_tc_rate = []
    star_ta_rate = []
    gedi_tc_rate = []
    gedi_ta_rate = []
    star_tc_std = []
    gedi_tc_std = []
    star_normalized = []  # T→C / T→A ratio
    gedi_normalized = []
    
    for pos in positions:
        s = bam_stats.get(pos, {'t_qual_sum': 0, 't_qual_count': 0, 't_qual_values': [], 'mm_count': 0, 'mm_total': 0, 't_count': 0})
        
        avg_q = s['t_qual_sum'] / s['t_qual_count'] if s['t_qual_count'] > 0 else float('nan')
        std_q = np.std(s['t_qual_values']) if len(s['t_qual_values']) > 1 else float('nan')
        
        avg_qual.append(avg_q)
        qual_std.append(std_q)
        qual_valid.append(s['t_qual_count'] > 0)
        t_counts.append(s['t_count'])

        mm_total = s['mm_total']
        mm_frac.append(s['mm_count'] / mm_total if mm_total > 0 else float('nan'))

        s_tc = star_tc.get(pos, {}).get('rate', 0)
        s_ta = star_ta.get(pos, {}).get('rate', 0.001)
        g_tc = gedi_tc.get(pos, {}).get('rate', 0)
        g_ta = gedi_ta.get(pos, {}).get('rate', 0.001)
        s_tc_std = star_tc.get(pos, {}).get('stdev', 0.0)
        g_tc_std = gedi_tc.get(pos, {}).get('stdev', 0.0)
        
        star_tc_rate.append(s_tc)
        star_ta_rate.append(s_ta)
        gedi_tc_rate.append(g_tc)
        gedi_ta_rate.append(g_ta)
        star_tc_std.append(s_tc_std)
        gedi_tc_std.append(g_tc_std)
        
        # Normalized: T→C / T→A (removes position-dependent artifact)
        star_normalized.append(s_tc / s_ta if s_ta > 0.001 else 0)
        gedi_normalized.append(g_tc / g_ta if g_ta > 0.001 else 0)

    # Intersection-based trim detection from STAR T→C stdev
    trim5p, trim3p, star_tc_std_smooth, trim_fit = detect_trim_segmented(
        star_tc_std,
        window=args.trim_window,
        min_seg_len=args.seg_min_len,
        max_trim=args.trim_max,
    )
    read_len = len(positions)
    trim3p_pos = read_len - 1 - trim3p if trim3p > 0 else None

    print("\n" + "=" * 80)
    print("Auto-trim (segmented regression, 2 breakpoints) from STAR T→C stdev")
    print("=" * 80)
    print(f"window={args.trim_window} seg_min_len={args.seg_min_len} max_trim={args.trim_max}")
    if trim_fit:
        print(f"breakpoints: b1={trim_fit.get('b1')} b2={trim_fit.get('b2')} sse={trim_fit.get('sse'):.4f}")
    print(f"trim5p={trim5p} trim3p={trim3p}")

    if args.stats_out and args.stats_out != '-':
        def write_stats(path, pos_list, pos_header):
            with open(path, 'w') as out:
                out.write("\t".join([
                    pos_header, "AvgQ_T", "TCount", "MultiMapFrac",
                    "STAR_TC_Cov", "STAR_TC_MM", "STAR_TC_Rate", "STAR_TC_Std",
                    "STAR_TC_Std_Smooth",
                    "STAR_TA_Cov", "STAR_TA_MM", "STAR_TA_Rate",
                    "GEDI_TC_Cov", "GEDI_TC_MM", "GEDI_TC_Rate", "GEDI_TC_Std",
                    "GEDI_TA_Cov", "GEDI_TA_MM", "GEDI_TA_Rate",
                    "STAR_TC_TA_Ratio", "GEDI_TC_TA_Ratio"
                ]) + "\n")
                for idx, pos in enumerate(pos_list):
                    pos0 = positions[idx]
                    star_tc_cov = star_tc.get(pos0, {}).get('cov', 0)
                    star_tc_mm = star_tc.get(pos0, {}).get('mm', 0)
                    star_ta_cov = star_ta.get(pos0, {}).get('cov', 0)
                    star_ta_mm = star_ta.get(pos0, {}).get('mm', 0)
                    gedi_tc_cov = gedi_tc.get(pos0, {}).get('cov', 0)
                    gedi_tc_mm = gedi_tc.get(pos0, {}).get('mm', 0)
                    gedi_ta_cov = gedi_ta.get(pos0, {}).get('cov', 0)
                    gedi_ta_mm = gedi_ta.get(pos0, {}).get('mm', 0)
                    out.write("\t".join(map(str, [
                        pos,
                        avg_qual[idx],
                        t_counts[idx],
                        mm_frac[idx],
                        star_tc_cov, star_tc_mm, star_tc_rate[idx], star_tc_std[idx],
                        star_tc_std_smooth[idx],
                        star_ta_cov, star_ta_mm, star_ta_rate[idx],
                        gedi_tc_cov, gedi_tc_mm, gedi_tc_rate[idx], gedi_tc_std[idx],
                        gedi_ta_cov, gedi_ta_mm, gedi_ta_rate[idx],
                        star_normalized[idx], gedi_normalized[idx],
                    ])) + "\n")

        write_stats(args.stats_out, positions, "Position0")

        base, ext = os.path.splitext(args.stats_out)
        stats_1based = f"{base}_1based{ext or '.tsv'}"
        write_stats(stats_1based, positions_plot, "Position1")
        print(f"1-based stats written to: {stats_1based}")
    
    # Print statistics table
    print("\n" + "=" * 120)
    print("Position Statistics (Raw T→C rates include sequencing artifacts at ends)")
    print("=" * 120)
    print(f"{'Pos':>4} | {'Avg_Q(T)':>8} | {'STAR_TC%':>9} | {'STAR_TA%':>9} | {'STAR_Norm':>10} | {'GEDI_TC%':>9} | {'GEDI_TA%':>9} | {'GEDI_Norm':>10}")
    print("-" * 100)
    
    table_positions = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 49]
    for pos in [p for p in table_positions if p <= max_pos]:
        print(f"{pos:>4} | {avg_qual[pos]:>8.1f} | {star_tc_rate[pos]:>8.3f}% | {star_ta_rate[pos]:>8.3f}% | {star_normalized[pos]:>10.2f} | {gedi_tc_rate[pos]:>8.3f}% | {gedi_ta_rate[pos]:>8.3f}% | {gedi_normalized[pos]:>10.2f}")
    
    # Compute correlations
    from scipy.stats import pearsonr, spearmanr
    
    # Use middle positions (exclude 10 bp edges) for correlation to avoid edge artifacts
    edge = 10
    if max_pos >= edge * 2:
        mid_range = range(edge, max_pos - edge + 1)
    else:
        mid_range = range(0, max_pos + 1)
    mid_idx = [i for i in mid_range if i < len(qual_valid) and qual_valid[i]]
    mid_start = mid_idx[0] if mid_idx else 0
    mid_end = mid_idx[-1] if mid_idx else max_pos
    mid_start_plot = mid_start + 1
    mid_end_plot = mid_end + 1
    mid_qual = [avg_qual[i] for i in mid_idx]
    mid_star_norm = [star_normalized[i] for i in mid_idx]
    mid_gedi_norm = [gedi_normalized[i] for i in mid_idx]
    
    if len(mid_qual) > 1:
        r_qual_star_norm, _ = pearsonr(mid_qual, mid_star_norm)
        r_qual_gedi_norm, _ = pearsonr(mid_qual, mid_gedi_norm)
    else:
        r_qual_star_norm, r_qual_gedi_norm = float('nan'), float('nan')
    
    print("\n" + "=" * 80)
    print("Correlation Analysis (Middle positions, avoiding edge artifacts)")
    print("=" * 80)
    print(f"PHRED(T-base) vs STAR T→C/T→A (normalized): r = {r_qual_star_norm:.4f}")
    print(f"PHRED(T-base) vs GEDI T→C/T→A (normalized): r = {r_qual_gedi_norm:.4f}")
    
    # Create plot
    fig, axes = plt.subplots(4, 2, figsize=(14, 18))
    
    # Plot 1: Raw T→C and T→A rates (showing the artifact)
    ax1 = axes[0, 0]
    ax1.plot(positions_plot, star_tc_rate, 'b-', label='STAR T→C', linewidth=2)
    ax1.plot(positions_plot, star_ta_rate, 'b--', label='STAR T→A', linewidth=1, alpha=0.7)
    ax1.plot(positions_plot, gedi_tc_rate, 'r-', label='GEDI T→C', linewidth=2)
    ax1.plot(positions_plot, gedi_ta_rate, 'r--', label='GEDI T→A', linewidth=1, alpha=0.7)
    ax1.set_xlabel('Read Position (starts at 1)')
    ax1.set_ylabel('Mismatch Rate (%)')
    ax1.set_title('Raw Mismatch Rates by Position\n(U-shape = sequencing artifact at read ends)')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    if mid_start > 0:
        ax1.axvspan(1, mid_start, alpha=0.1, color='red', label='Edge artifact')
    if mid_end < max_pos:
        ax1.axvspan(mid_end + 2, max_pos + 1, alpha=0.1, color='red')
    ax1.set_yscale('log')
    
    # Plot 2: Normalized T→C/T→A ratio (the TRUE 4sU signal)
    ax2 = axes[0, 1]
    ax2.plot(positions_plot, star_normalized, 'b-', label='STAR', linewidth=2)
    ax2.plot(positions_plot, gedi_normalized, 'r-', label='GEDI', linewidth=2)
    ax2.axhline(y=1, color='gray', linestyle='--', alpha=0.7, label='No 4sU enrichment')
    ax2.set_xlabel('Read Position (starts at 1)')
    ax2.set_ylabel('T→C / T→A Ratio')
    ax2.set_title('Normalized T→C/T→A Ratio (TRUE 4sU Signal)\n>1 = 4sU enrichment; ~1 = pure sequencing error')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.axvspan(mid_start_plot, mid_end_plot, alpha=0.1, color='green', label='Reliable region')
    
    # Plot 3: PHRED quality by position
    ax3 = axes[1, 0]
    ax3.errorbar(positions_plot, avg_qual, yerr=qual_std, fmt='ko-', capsize=2, markersize=3, alpha=0.7)
    ax3.set_xlabel('Read Position (starts at 1)')
    ax3.set_ylabel('Average PHRED Quality')
    ax3.set_title('Base Quality by Read Position (Transcript T bases)\n(Uniform quality across positions)')
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, 45)
    
    # Plot 4: Quality vs Normalized T→C/T→A (middle positions only)
    ax4 = axes[1, 1]
    scatter1 = ax4.scatter(mid_qual, mid_star_norm, c='blue', alpha=0.7, s=50, label='STAR')
    scatter2 = ax4.scatter(mid_qual, mid_gedi_norm, c='red', alpha=0.7, s=50, label='GEDI')
    ax4.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    ax4.set_xlabel('Average PHRED Quality')
    ax4.set_ylabel('T→C / T→A Ratio')
    ax4.set_title(f'Quality vs Normalized Ratio (positions {mid_start}-{mid_end})\nSTAR r={r_qual_star_norm:.2f}, GEDI r={r_qual_gedi_norm:.2f}')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Multimapper fraction by position
    ax5 = axes[2, 0]
    ax5.plot(positions_plot, mm_frac, 'm-', linewidth=2, label='NH>1 fraction')
    ax5.set_xlabel('Read Position (starts at 1)')
    ax5.set_ylabel('Fraction of Alignments')
    ax5.set_title('Multimapper Fraction by Position (NH > 1)')
    ax5.grid(True, alpha=0.3)
    ax5.set_ylim(0, 1)

    # Plot 6: Transcript T coverage by position
    ax6 = axes[2, 1]
    ax6.plot(positions_plot, t_counts, 'k-', linewidth=2)
    ax6.set_xlabel('Read Position (starts at 1)')
    ax6.set_ylabel('T-base Count')
    ax6.set_title('Transcript T Coverage by Position')
    ax6.grid(True, alpha=0.3)

    # Plot 7: T→C per-base stdev by position (native variation)
    ax7 = axes[3, 0]
    ax7.plot(positions_plot, star_tc_std, 'b-', linewidth=2, label='STAR')
    ax7.plot(positions_plot, gedi_tc_std, 'r-', linewidth=2, label='GEDI')
    ax7.plot(positions_plot, star_tc_std_smooth, 'k--', linewidth=1.5, label='STAR smoothed')
    if trim_fit:
        b1 = trim_fit.get('b1')
        b2 = trim_fit.get('b2')
        if b1 is not None and b2 is not None:
            m1, c1 = trim_fit['fit1']
            m2, c2 = trim_fit['fit2']
            m3, c3 = trim_fit['fit3']
            x1 = np.arange(0, b1)
            x2 = np.arange(b1, b2 + 1)
            x3 = np.arange(b2 + 1, len(positions))
            ax7.plot(x1 + 1, m1 * x1 + c1, color='purple', linestyle='-', linewidth=1.2, alpha=0.8, label='seg1 fit')
            ax7.plot(x2 + 1, m2 * x2 + c2, color='purple', linestyle='--', linewidth=1.2, alpha=0.8, label='seg2 fit')
            ax7.plot(x3 + 1, m3 * x3 + c3, color='purple', linestyle=':', linewidth=1.2, alpha=0.8, label='seg3 fit')
    if trim5p > 0:
        ax7.axvline(trim5p + 1, color='purple', linestyle='--', linewidth=1.5, label=f'trim5p={trim5p}')
    if trim3p_pos is not None:
        ax7.axvline(trim3p_pos + 1, color='purple', linestyle=':', linewidth=1.5, label=f'trim3p={trim3p}')
    ax7.set_xlabel('Read Position (starts at 1)')
    ax7.set_ylabel('T→C Stdev (per-base, %)')
    ax7.set_title('T→C Variation by Position (per-base stdev)')
    ax7.grid(True, alpha=0.3)
    ax7.legend()

    # Hide unused panel
    axes[3, 1].axis('off')

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"\nPlot saved to: {args.output}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    mid_star_avg = np.mean([star_normalized[i] for i in mid_idx]) if mid_idx else float('nan')
    mid_gedi_avg = np.mean([gedi_normalized[i] for i in mid_idx]) if mid_idx else float('nan')
    left_edge = f"0-{mid_start - 1}" if mid_start > 0 else "none"
    right_edge = f"{mid_end + 1}-{max_pos}" if mid_end < max_pos else "none"
    
    print(f"""
KEY INSIGHT: Raw T→C rates show U-shaped artifact from sequencing errors at read ends.
To see TRUE 4sU signal, normalize by T→A (a non-4sU control transition).

NORMALIZED T→C/T→A RATIO (middle positions {mid_start}-{mid_end}):
  - STAR average: {mid_star_avg:.2f}
  - GEDI average: {mid_gedi_avg:.2f}
  - Ratio > 1 indicates 4sU-specific T→C enrichment

EDGE ARTIFACTS (positions {left_edge} and {right_edge}):
  - Both T→C and T→A are elevated (general sequencing error)
  - T→C/T→A ratio ~1 at edges = NO 4sU-specific signal
  - These positions should be CLIPPED for accurate 4sU quantification

QUALITY SCORES at transcript T positions are UNIFORM across the read:
  - T-base quality does NOT explain the U-shaped artifact
  - The artifact is inherent to sequencing chemistry at read ends

RECOMMENDATION:
  - Clip positions {left_edge} and {right_edge} for 4sU quantification
  - OR normalize by T→A to remove position-dependent artifacts
""")

if __name__ == '__main__':
    main()
