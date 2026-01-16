#!/usr/bin/env python3
"""
Compare quantification results between Salmon and em_quant.

Usage:
    compare_quant.py <salmon_quant.sf> <em_quant.tsv> [--tolerance <tol>]

Compares NumReads (estimated counts) between the two files and reports
differences exceeding the tolerance threshold.
"""

import sys
import argparse

def parse_tsv(filename):
    """Parse TSV file and return dict mapping transcript name to NumReads."""
    results = {}
    with open(filename, 'r') as f:
        header = f.readline().strip().split('\t')
        
        # Find column indices
        name_idx = header.index('Name')
        numreads_idx = header.index('NumReads')
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) > max(name_idx, numreads_idx):
                name = fields[name_idx]
                numreads = float(fields[numreads_idx])
                results[name] = numreads
    return results

def compare_results(salmon_results, em_results, tolerance=1e-4, near_zero_threshold=1.0):
    """Compare results and return list of differences.
    
    Args:
        salmon_results: dict of transcript -> NumReads from Salmon
        em_results: dict of transcript -> NumReads from em_quant
        tolerance: relative tolerance for significant transcripts
        near_zero_threshold: values below this are considered "near-zero"
    
    Returns:
        tuple of (significant_diffs, near_zero_diffs)
    """
    significant_diffs = []
    near_zero_diffs = []
    all_transcripts = set(salmon_results.keys()) | set(em_results.keys())
    
    for txp in sorted(all_transcripts):
        salmon_val = salmon_results.get(txp, 0.0)
        em_val = em_results.get(txp, 0.0)
        
        # Skip if both are zero
        if salmon_val == 0.0 and em_val == 0.0:
            continue
        
        # Handle near-zero cases separately
        # If both values are below threshold, consider them "effectively zero"
        if max(salmon_val, em_val) < near_zero_threshold:
            if salmon_val == 0.0 and em_val > 0:
                near_zero_diffs.append({
                    'transcript': txp,
                    'salmon': salmon_val,
                    'em': em_val,
                    'abs_diff': abs(em_val - salmon_val)
                })
            continue
        
        # Compute relative difference for significant values
        if salmon_val != 0.0:
            rel_diff = abs(em_val - salmon_val) / abs(salmon_val)
        elif em_val != 0.0:
            rel_diff = abs(em_val - salmon_val) / abs(em_val)
        else:
            rel_diff = 0.0
        
        if rel_diff > tolerance:
            significant_diffs.append({
                'transcript': txp,
                'salmon': salmon_val,
                'em': em_val,
                'rel_diff': rel_diff
            })
    
    return significant_diffs, near_zero_diffs

def main():
    parser = argparse.ArgumentParser(description='Compare Salmon and em_quant results')
    parser.add_argument('salmon_file', help='Salmon quant.sf file')
    parser.add_argument('em_file', help='em_quant output TSV file')
    parser.add_argument('--tolerance', type=float, default=1e-4,
                       help='Relative tolerance threshold (default: 1e-4)')
    parser.add_argument('--near-zero', type=float, default=1.0,
                       help='Threshold for near-zero values (default: 1.0)')
    
    args = parser.parse_args()
    
    # Parse files
    print(f"Loading Salmon results from: {args.salmon_file}")
    salmon_results = parse_tsv(args.salmon_file)
    print(f"  Found {len(salmon_results)} transcripts")
    
    print(f"Loading em_quant results from: {args.em_file}")
    em_results = parse_tsv(args.em_file)
    print(f"  Found {len(em_results)} transcripts")
    
    # Compare
    print(f"\nComparing with tolerance: {args.tolerance} (near-zero threshold: {args.near_zero})")
    significant_diffs, near_zero_diffs = compare_results(
        salmon_results, em_results, args.tolerance, args.near_zero
    )
    
    # Report near-zero differences (informational, not failure)
    if near_zero_diffs:
        print(f"\nNote: {len(near_zero_diffs)} transcripts have near-zero differences (Salmon=0, em_quant<{args.near_zero})")
    
    if not significant_diffs:
        print("\n✓ PASS: All significant transcripts within tolerance")
        return 0
    else:
        print(f"\n✗ FAIL: {len(significant_diffs)} transcripts exceed tolerance")
        print("\nSignificant differences:")
        print(f"{'Transcript':<30} {'Salmon':>15} {'em_quant':>15} {'Rel Diff':>15}")
        print("-" * 75)
        for diff in significant_diffs[:20]:  # Show first 20
            print(f"{diff['transcript']:<30} {diff['salmon']:>15.6f} "
                  f"{diff['em']:>15.6f} {diff['rel_diff']:>15.6e}")
        if len(significant_diffs) > 20:
            print(f"... and {len(significant_diffs) - 20} more")
        return 1

if __name__ == '__main__':
    sys.exit(main())
