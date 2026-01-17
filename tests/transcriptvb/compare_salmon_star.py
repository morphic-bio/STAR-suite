#!/usr/bin/env python3
"""
compare_salmon_star.py - Compare TranscriptVB output with Salmon

Usage:
    python3 compare_salmon_star.py salmon_quant.sf star_quant.sf [--verbose]

Exit codes:
    0 - PASS (correlation thresholds met)
    1 - FAIL (correlation below threshold)
    2 - ERROR (file not found, etc.)
"""

import pandas as pd
import sys
import argparse
from scipy.stats import spearmanr, pearsonr
import numpy as np


def load_quant_sf(filepath):
    """Load a quant.sf file into a DataFrame."""
    try:
        df = pd.read_csv(filepath, sep='\t')
        required_cols = ['Name', 'Length', 'EffectiveLength', 'TPM', 'NumReads']
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Missing required column: {col}")
        return df
    except Exception as e:
        print(f"Error loading {filepath}: {e}", file=sys.stderr)
        sys.exit(2)


def compare_quantifications(salmon_df, star_df, verbose=False):
    """Compare two quantification results."""
    
    # Merge on transcript name
    merged = salmon_df.merge(star_df, on='Name', suffixes=('_salmon', '_star'))
    
    results = {
        'total_transcripts': len(merged),
        'salmon_expressed': (merged['NumReads_salmon'] > 0).sum(),
        'star_expressed': (merged['NumReads_star'] > 0).sum(),
    }
    
    # Jointly expressed
    expressed = merged[(merged['NumReads_salmon'] > 0) & (merged['NumReads_star'] > 0)]
    results['jointly_expressed'] = len(expressed)
    
    # Only in one
    results['salmon_only'] = ((merged['NumReads_salmon'] > 0) & (merged['NumReads_star'] == 0)).sum()
    results['star_only'] = ((merged['NumReads_salmon'] == 0) & (merged['NumReads_star'] > 0)).sum()
    
    # Correlations - all transcripts
    results['spearman_all'], _ = spearmanr(merged['NumReads_salmon'], merged['NumReads_star'])
    results['pearson_all'], _ = pearsonr(merged['NumReads_salmon'], merged['NumReads_star'])
    
    # Correlations - jointly expressed
    if len(expressed) > 2:
        results['spearman_expressed'], _ = spearmanr(expressed['NumReads_salmon'], expressed['NumReads_star'])
        results['pearson_expressed'], _ = pearsonr(expressed['NumReads_salmon'], expressed['NumReads_star'])
    else:
        results['spearman_expressed'] = np.nan
        results['pearson_expressed'] = np.nan
    
    # Total reads
    results['salmon_total_reads'] = merged['NumReads_salmon'].sum()
    results['star_total_reads'] = merged['NumReads_star'].sum()
    
    # Max difference
    merged['abs_diff'] = abs(merged['NumReads_salmon'] - merged['NumReads_star'])
    results['max_diff'] = merged['abs_diff'].max()
    results['mean_diff'] = merged['abs_diff'].mean()
    
    if verbose:
        # Top differences
        top_diff = merged.nlargest(10, 'abs_diff')[['Name', 'NumReads_salmon', 'NumReads_star', 'abs_diff']]
        results['top_differences'] = top_diff
    
    return results


def print_results(results, verbose=False):
    """Print comparison results."""
    print("=" * 50)
    print("Salmon vs TranscriptVB Comparison")
    print("=" * 50)
    print()
    
    print("Transcript Counts:")
    print(f"  Total transcripts:    {results['total_transcripts']}")
    print(f"  Salmon expressed:     {results['salmon_expressed']}")
    print(f"  STAR expressed:       {results['star_expressed']}")
    print(f"  Jointly expressed:    {results['jointly_expressed']}")
    print(f"  Salmon only:          {results['salmon_only']}")
    print(f"  STAR only:            {results['star_only']}")
    print()
    
    print("Correlations (All Transcripts):")
    print(f"  Spearman:  {results['spearman_all']:.6f}")
    print(f"  Pearson:   {results['pearson_all']:.6f}")
    print()
    
    print("Correlations (Jointly Expressed):")
    if not np.isnan(results['spearman_expressed']):
        print(f"  Spearman:  {results['spearman_expressed']:.6f}")
        print(f"  Pearson:   {results['pearson_expressed']:.6f}")
    else:
        print("  (insufficient data)")
    print()
    
    print("Read Counts:")
    print(f"  Salmon total:  {results['salmon_total_reads']:.1f}")
    print(f"  STAR total:    {results['star_total_reads']:.1f}")
    print(f"  Difference:    {abs(results['salmon_total_reads'] - results['star_total_reads']):.1f}")
    print()
    
    print("Differences:")
    print(f"  Max difference:   {results['max_diff']:.3f}")
    print(f"  Mean difference:  {results['mean_diff']:.3f}")
    
    if verbose and 'top_differences' in results:
        print()
        print("Top 10 Differences:")
        print(results['top_differences'].to_string(index=False))
    
    print()


def check_thresholds(results):
    """Check if results meet quality thresholds."""
    thresholds = {
        'spearman_all': 0.95,
        'spearman_expressed': 0.99,
        'pearson_expressed': 0.99,
    }
    
    passed = True
    print("Threshold Checks:")
    
    for metric, threshold in thresholds.items():
        value = results.get(metric, np.nan)
        if np.isnan(value):
            status = "SKIP"
        elif value >= threshold:
            status = "PASS"
        else:
            status = "FAIL"
            passed = False
        print(f"  {metric}: {value:.4f} >= {threshold} [{status}]")
    
    print()
    return passed


def main():
    parser = argparse.ArgumentParser(description='Compare Salmon and TranscriptVB quantification')
    parser.add_argument('salmon_file', help='Salmon quant.sf file')
    parser.add_argument('star_file', help='STAR TranscriptVB quant.sf file')
    parser.add_argument('--verbose', '-v', action='store_true', help='Show detailed output')
    parser.add_argument('--gene', action='store_true', help='Compare gene-level files (quant.genes.sf)')
    args = parser.parse_args()
    
    # Load files
    salmon_df = load_quant_sf(args.salmon_file)
    star_df = load_quant_sf(args.star_file)
    
    # Compare
    results = compare_quantifications(salmon_df, star_df, verbose=args.verbose)
    
    # Print results
    print_results(results, verbose=args.verbose)
    
    # Check thresholds
    passed = check_thresholds(results)
    
    if passed:
        print("✓ OVERALL: PASS")
        return 0
    else:
        print("✗ OVERALL: FAIL")
        return 1


if __name__ == "__main__":
    sys.exit(main())

