#!/usr/bin/env python3
"""
Analyze mismatch position histograms to detect systematic position filtering differences.

Compares T→C mismatch counts by read position between STAR and GEDI.
"""

import sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np

def parse_mismatch_details(path: str, tool_name: str = "STAR"):
    """
    Parse mismatch details file.
    
    Expected format:
    Category	Genomic	Read	Position	Overlap	Opposite	Coverage	Mismatches
    Exonic	T	C	0	0	0	219629	16452.7
    
    Returns: dict mapping position -> (coverage, mismatches)
    """
    position_data = defaultdict(lambda: {'coverage': 0, 'mismatches': 0.0})
    
    with open(path) as f:
        header = f.readline().strip().split('\t')
        
        # Find column indices
        try:
            genomic_idx = header.index('Genomic')
            read_idx = header.index('Read')
            pos_idx = header.index('Position')
            cov_idx = header.index('Coverage')
            mismatch_idx = header.index('Mismatches')
        except ValueError as e:
            print(f"Error: Missing column in {tool_name} mismatch details: {e}")
            print(f"Header: {header}")
            sys.exit(1)
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(genomic_idx, read_idx, pos_idx, cov_idx, mismatch_idx):
                continue
            
            genomic = parts[genomic_idx]
            read = parts[read_idx]
            position = int(parts[pos_idx])
            coverage = float(parts[cov_idx])
            mismatches = float(parts[mismatch_idx])
            
            # Only count T→C mismatches
            if genomic == 'T' and read == 'C':
                position_data[position]['coverage'] += coverage
                position_data[position]['mismatches'] += mismatches
    
    return position_data

def compute_rate_by_position(position_data, max_pos=50):
    """Compute T→C rate by position."""
    rates = []
    positions = []
    coverages = []
    
    for pos in range(max_pos + 1):
        if pos in position_data:
            cov = position_data[pos]['coverage']
            mism = position_data[pos]['mismatches']
            rate = mism / cov if cov > 0 else 0
            rates.append(rate)
            positions.append(pos)
            coverages.append(cov)
        else:
            rates.append(0)
            positions.append(pos)
            coverages.append(0)
    
    return positions, rates, coverages

def main():
    star_path = 'test/tmp_prod_compat/default_SlamQuant.out.mismatchdetails.tsv'
    gedi_path = '/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3.mismatchdetails.tsv'
    
    print("=" * 80)
    print("Mismatch Position Histogram Analysis")
    print("=" * 80)
    
    print("\nParsing STAR mismatch details...")
    star_data = parse_mismatch_details(star_path, "STAR")
    
    print("Parsing GEDI mismatch details...")
    gedi_data = parse_mismatch_details(gedi_path, "GEDI")
    
    # Compute rates by position
    max_pos = max(max(star_data.keys()), max(gedi_data.keys())) if star_data and gedi_data else 50
    max_pos = min(max_pos, 50)  # Cap at 50 for readability
    
    star_pos, star_rates, star_covs = compute_rate_by_position(star_data, max_pos)
    gedi_pos, gedi_rates, gedi_covs = compute_rate_by_position(gedi_data, max_pos)
    
    # Summary statistics
    print(f"\nPosition range: 0-{max_pos}")
    print(f"\nSTAR T→C totals:")
    print(f"  Total coverage: {sum(star_covs):.0f}")
    print(f"  Total mismatches: {sum(star_data[p]['mismatches'] for p in star_data):.1f}")
    print(f"  Overall rate: {sum(star_data[p]['mismatches'] for p in star_data) / sum(star_covs) if sum(star_covs) > 0 else 0:.6f}")
    
    print(f"\nGEDI T→C totals:")
    print(f"  Total coverage: {sum(gedi_covs):.0f}")
    print(f"  Total mismatches: {sum(gedi_data[p]['mismatches'] for p in gedi_data):.1f}")
    print(f"  Overall rate: {sum(gedi_data[p]['mismatches'] for p in gedi_data) / sum(gedi_covs) if sum(gedi_covs) > 0 else 0:.6f}")
    
    # Compare early vs late positions
    early_positions = list(range(0, 5))
    late_positions = list(range(max_pos - 4, max_pos + 1))
    
    star_early_rate = sum(star_data[p]['mismatches'] for p in early_positions if p in star_data) / \
                      sum(star_data[p]['coverage'] for p in early_positions if p in star_data) if \
                      sum(star_data[p]['coverage'] for p in early_positions if p in star_data) > 0 else 0
    
    star_late_rate = sum(star_data[p]['mismatches'] for p in late_positions if p in star_data) / \
                     sum(star_data[p]['coverage'] for p in late_positions if p in star_data) if \
                     sum(star_data[p]['coverage'] for p in late_positions if p in star_data) > 0 else 0
    
    gedi_early_rate = sum(gedi_data[p]['mismatches'] for p in early_positions if p in gedi_data) / \
                      sum(gedi_data[p]['coverage'] for p in early_positions if p in gedi_data) if \
                      sum(gedi_data[p]['coverage'] for p in early_positions if p in gedi_data) > 0 else 0
    
    gedi_late_rate = sum(gedi_data[p]['mismatches'] for p in late_positions if p in gedi_data) / \
                     sum(gedi_data[p]['coverage'] for p in late_positions if p in gedi_data) if \
                     sum(gedi_data[p]['coverage'] for p in late_positions if p in gedi_data) > 0 else 0
    
    print(f"\nEarly positions (0-4):")
    print(f"  STAR rate: {star_early_rate:.6f}")
    print(f"  GEDI rate: {gedi_early_rate:.6f}")
    print(f"  Ratio (STAR/GEDI): {star_early_rate / gedi_early_rate if gedi_early_rate > 0 else 'inf':.3f}")
    
    print(f"\nLate positions ({max_pos-4}-{max_pos}):")
    print(f"  STAR rate: {star_late_rate:.6f}")
    print(f"  GEDI rate: {gedi_late_rate:.6f}")
    print(f"  Ratio (STAR/GEDI): {star_late_rate / gedi_late_rate if gedi_late_rate > 0 else 'inf':.3f}")
    
    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot 1: Rates by position
    ax1 = axes[0]
    ax1.plot(star_pos, star_rates, 'b-', label='STAR', linewidth=2, alpha=0.7)
    ax1.plot(gedi_pos, gedi_rates, 'r-', label='GEDI', linewidth=2, alpha=0.7)
    ax1.set_xlabel('Read Position')
    ax1.set_ylabel('T→C Rate')
    ax1.set_title('T→C Mismatch Rate by Read Position')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Coverage by position
    ax2 = axes[1]
    ax2.plot(star_pos, star_covs, 'b-', label='STAR', linewidth=2, alpha=0.7)
    ax2.plot(gedi_pos, gedi_covs, 'r-', label='GEDI', linewidth=2, alpha=0.7)
    ax2.set_xlabel('Read Position')
    ax2.set_ylabel('Coverage')
    ax2.set_title('Coverage by Read Position')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')
    
    plt.tight_layout()
    output_file = 'test/tmp_prod_compat/mismatch_position_histogram.png'
    plt.savefig(output_file, dpi=150)
    print(f"\nPlot saved to: {output_file}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    
    if abs(star_early_rate / gedi_early_rate - 1) > 0.2 if gedi_early_rate > 0 else False:
        print("⚠️  Early position rates differ significantly")
    if abs(star_late_rate / gedi_late_rate - 1) > 0.2 if gedi_late_rate > 0 else False:
        print("⚠️  Late position rates differ significantly")
    
    # Check for systematic trimming
    star_first_nonzero = next((i for i, r in enumerate(star_rates) if r > 0), None)
    gedi_first_nonzero = next((i for i, r in enumerate(gedi_rates) if r > 0), None)
    
    if star_first_nonzero != gedi_first_nonzero:
        print(f"⚠️  First non-zero position differs: STAR={star_first_nonzero}, GEDI={gedi_first_nonzero}")
    
    star_last_nonzero = next((i for i in range(len(star_rates)-1, -1, -1) if star_rates[i] > 0), None)
    gedi_last_nonzero = next((i for i in range(len(gedi_rates)-1, -1, -1) if gedi_rates[i] > 0), None)
    
    if star_last_nonzero != gedi_last_nonzero:
        print(f"⚠️  Last non-zero position differs: STAR={star_last_nonzero}, GEDI={gedi_last_nonzero}")

if __name__ == '__main__':
    main()
