#!/usr/bin/env python3
"""
Analyze T→C mismatch position histograms by category (Exonic vs ExonicSense).
Removes intronic confounding.
"""

import sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def parse_mismatch_details(path: str, categories: list):
    """
    Parse mismatch details file for specific categories.
    
    Returns: dict[category] -> dict[position] -> (coverage, mismatches)
    """
    data = defaultdict(lambda: defaultdict(lambda: {'cov': 0, 'mm': 0.0}))
    
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
            
            # Only T→C mismatches
            if genomic == 'T' and read == 'C' and cat in categories:
                data[cat][pos]['cov'] = cov
                data[cat][pos]['mm'] = mm
    
    return data

def compute_rates(data, max_pos=50):
    """Compute T→C rate by position."""
    rates = {}
    for cat in data:
        rates[cat] = {}
        for pos in range(max_pos):
            if pos in data[cat]:
                cov = data[cat][pos]['cov']
                mm = data[cat][pos]['mm']
                rates[cat][pos] = mm / cov if cov > 0 else 0
            else:
                rates[cat][pos] = 0
    return rates

def main():
    star_path = 'test/tmp_prod_compat/default_SlamQuant.out.mismatchdetails.tsv'
    gedi_path = '/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3_Sense_fixed.mismatchdetails.tsv'
    
    categories = ['Exonic', 'ExonicSense']
    
    print("=" * 100)
    print("Position Histogram Analysis - Exonic/ExonicSense Only (No Intronic)")
    print("=" * 100)
    
    star_data = parse_mismatch_details(star_path, categories)
    gedi_data = parse_mismatch_details(gedi_path, categories)
    
    star_rates = compute_rates(star_data)
    gedi_rates = compute_rates(gedi_data)
    
    # Print comparison tables
    for cat in categories:
        print(f"\n{'='*80}")
        print(f"Category: {cat}")
        print("="*80)
        print(f"{'Pos':>4} | {'STAR_Cov':>10} | {'STAR_MM':>10} | {'STAR_Rate':>10} | {'GEDI_Cov':>10} | {'GEDI_MM':>10} | {'GEDI_Rate':>10} | {'Ratio':>8}")
        print("-"*90)
        
        for pos in range(50):
            s_cov = star_data[cat].get(pos, {}).get('cov', 0)
            s_mm = star_data[cat].get(pos, {}).get('mm', 0)
            g_cov = gedi_data[cat].get(pos, {}).get('cov', 0)
            g_mm = gedi_data[cat].get(pos, {}).get('mm', 0)
            
            s_rate = s_mm / s_cov * 100 if s_cov > 0 else 0
            g_rate = g_mm / g_cov * 100 if g_cov > 0 else 0
            ratio = s_rate / g_rate if g_rate > 0 else float('inf')
            
            # Only print significant positions
            if pos in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 49] or ratio > 3 or ratio < 0.5:
                print(f"{pos:>4} | {s_cov:>10.0f} | {s_mm:>10.1f} | {s_rate:>9.4f}% | {g_cov:>10.0f} | {g_mm:>10.1f} | {g_rate:>9.4f}% | {ratio:>8.2f}")
    
    # Create comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    for idx, cat in enumerate(categories):
        ax_rate = axes[idx, 0]
        ax_cov = axes[idx, 1]
        
        positions = list(range(50))
        star_r = [star_rates[cat].get(p, 0) * 100 for p in positions]
        gedi_r = [gedi_rates[cat].get(p, 0) * 100 for p in positions]
        
        star_c = [star_data[cat].get(p, {}).get('cov', 0) for p in positions]
        gedi_c = [gedi_data[cat].get(p, {}).get('cov', 0) for p in positions]
        
        # Rate plot
        ax_rate.plot(positions, star_r, 'b-', label='STAR', linewidth=2, alpha=0.7)
        ax_rate.plot(positions, gedi_r, 'r-', label='GEDI', linewidth=2, alpha=0.7)
        ax_rate.set_xlabel('Read Position')
        ax_rate.set_ylabel('T→C Rate (%)')
        ax_rate.set_title(f'{cat}: T→C Rate by Position')
        ax_rate.legend()
        ax_rate.grid(True, alpha=0.3)
        ax_rate.set_yscale('log')
        
        # Coverage plot
        ax_cov.plot(positions, star_c, 'b-', label='STAR', linewidth=2, alpha=0.7)
        ax_cov.plot(positions, gedi_c, 'r-', label='GEDI', linewidth=2, alpha=0.7)
        ax_cov.set_xlabel('Read Position')
        ax_cov.set_ylabel('T Coverage')
        ax_cov.set_title(f'{cat}: T Coverage by Position')
        ax_cov.legend()
        ax_cov.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = 'test/tmp_prod_compat/position_histogram_exonic_only.png'
    plt.savefig(output_file, dpi=150)
    print(f"\nPlot saved to: {output_file}")
    
    # Summary statistics
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    for cat in categories:
        print(f"\n{cat}:")
        
        # Early positions (0-10)
        s_early_mm = sum(star_data[cat].get(p, {}).get('mm', 0) for p in range(11))
        s_early_cov = sum(star_data[cat].get(p, {}).get('cov', 0) for p in range(11))
        g_early_mm = sum(gedi_data[cat].get(p, {}).get('mm', 0) for p in range(11))
        g_early_cov = sum(gedi_data[cat].get(p, {}).get('cov', 0) for p in range(11))
        
        s_early_rate = s_early_mm / s_early_cov * 100 if s_early_cov > 0 else 0
        g_early_rate = g_early_mm / g_early_cov * 100 if g_early_cov > 0 else 0
        
        # Late positions (40-49)
        s_late_mm = sum(star_data[cat].get(p, {}).get('mm', 0) for p in range(40, 50))
        s_late_cov = sum(star_data[cat].get(p, {}).get('cov', 0) for p in range(40, 50))
        g_late_mm = sum(gedi_data[cat].get(p, {}).get('mm', 0) for p in range(40, 50))
        g_late_cov = sum(gedi_data[cat].get(p, {}).get('cov', 0) for p in range(40, 50))
        
        s_late_rate = s_late_mm / s_late_cov * 100 if s_late_cov > 0 else 0
        g_late_rate = g_late_mm / g_late_cov * 100 if g_late_cov > 0 else 0
        
        print(f"  Early (0-10):  STAR={s_early_rate:.4f}%, GEDI={g_early_rate:.4f}%, Ratio={s_early_rate/g_early_rate if g_early_rate > 0 else 0:.2f}")
        print(f"  Late (40-49):  STAR={s_late_rate:.4f}%, GEDI={g_late_rate:.4f}%, Ratio={s_late_rate/g_late_rate if g_late_rate > 0 else 0:.2f}")

if __name__ == '__main__':
    main()
