#!/usr/bin/env python3
"""
Compare our expected GC distribution with Salmon's.

Salmon stores expected GC in binary GCFragModel format:
- Conditional bins (rows) × GC bins (columns) matrix
- We need to extract the marginal distribution (sum over rows)
"""

import struct
import gzip
import sys
import numpy as np

def parse_salmon_gc_binary(filename):
    """Parse Salmon's binary GCFragModel format."""
    with gzip.open(filename, 'rb') as f:
        # Read distribution space (0=LOG, 1=LINEAR)
        dtype_bytes = f.read(4)
        dtype = struct.unpack('i', dtype_bytes)[0]
        
        # Read matrix dimensions
        rows_bytes = f.read(8)  # Eigen::Index is typically int64_t
        cols_bytes = f.read(8)
        rows = struct.unpack('q', rows_bytes)[0]
        cols = struct.unpack('q', cols_bytes)[0]
        
        # Read row totals (modelTotals_)
        totals = []
        for i in range(rows):
            total_bytes = f.read(8)  # double
            total = struct.unpack('d', total_bytes)[0]
            totals.append(total)
        
        # Read matrix (counts_)
        matrix = []
        for r in range(rows):
            row = []
            for c in range(cols):
                val_bytes = f.read(8)  # double
                val = struct.unpack('d', val_bytes)[0]
                row.append(val)
            matrix.append(row)
        
        return {
            'dtype': dtype,  # 0=LOG, 1=LINEAR
            'rows': rows,     # Conditional bins
            'cols': cols,     # GC bins (should be 101)
            'totals': totals,
            'matrix': matrix
        }

def extract_marginal_distribution(gc_model):
    """Extract marginal GC distribution (sum over conditional bins)."""
    matrix = np.array(gc_model['matrix'])
    
    # Sum over rows (conditional bins) to get marginal GC distribution
    marginal = matrix.sum(axis=0)
    
    # Normalize
    total = marginal.sum()
    if total > 0:
        marginal = marginal / total
    
    return marginal

def load_our_expected_gc(filename):
    """Load our TSV format expected GC."""
    gc_dist = np.zeros(101)
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                gc_pct = int(parts[0])
                prob = float(parts[1])
                if 0 <= gc_pct < 101:
                    gc_dist[gc_pct] = prob
    return gc_dist

def compare_distributions(salmon_dist, our_dist):
    """Compare two GC distributions."""
    print("=" * 70)
    print("GC Distribution Comparison")
    print("=" * 70)
    
    salmon_bins = len(salmon_dist)
    our_bins = len(our_dist)
    
    print(f"\nSalmon distribution:")
    print(f"  Bins: {salmon_bins}")
    print(f"  Sum: {salmon_dist.sum():.10f}")
    print(f"  Non-zero bins: {(salmon_dist > 0).sum()}")
    salmon_gc_range = np.arange(salmon_bins) * (100.0 / (salmon_bins - 1)) if salmon_bins > 1 else np.array([0])
    print(f"  Mean GC%: {(salmon_dist * salmon_gc_range).sum():.2f}")
    
    print(f"\nOur distribution:")
    print(f"  Bins: {our_bins}")
    print(f"  Sum: {our_dist.sum():.10f}")
    print(f"  Non-zero bins: {(our_dist > 0).sum()}")
    our_gc_range = np.arange(our_bins)
    print(f"  Mean GC%: {(our_dist * our_gc_range).sum():.2f}")
    
    # Interpolate Salmon's distribution to 101 bins for comparison
    if salmon_bins != our_bins:
        print(f"\nNote: Salmon uses {salmon_bins} bins, we use {our_bins} bins")
        print("Interpolating Salmon's distribution to {our_bins} bins for comparison...")
        
        # Create interpolation: map Salmon bins to our bins
        salmon_gc_values = np.linspace(0, 100, salmon_bins)
        our_gc_values = np.arange(our_bins)
        
        # Interpolate
        from scipy.interpolate import interp1d
        interp_func = interp1d(salmon_gc_values, salmon_dist, kind='linear', 
                              bounds_error=False, fill_value=0.0)
        salmon_dist_interp = interp_func(our_gc_values)
        salmon_dist_interp = salmon_dist_interp / salmon_dist_interp.sum()  # Renormalize
        
        salmon_dist = salmon_dist_interp
        salmon_bins = our_bins
    
    # Compute differences
    diff = np.abs(salmon_dist - our_dist)
    max_diff = diff.max()
    mean_diff = diff.mean()
    
    print(f"\nDifferences:")
    print(f"  Max absolute difference: {max_diff:.10f}")
    print(f"  Mean absolute difference: {mean_diff:.10f}")
    
    # Find bins with largest differences
    max_diff_idx = diff.argmax()
    print(f"\nLargest difference at GC%={max_diff_idx}:")
    print(f"  Salmon: {salmon_dist[max_diff_idx]:.10f}")
    print(f"  Ours:   {our_dist[max_diff_idx]:.10f}")
    print(f"  Diff:   {diff[max_diff_idx]:.10f}")
    
    # Check if distributions are close
    if max_diff < 0.01 and mean_diff < 0.001:
        print("\n✓ Distributions match closely!")
        return True
    elif max_diff < 0.1 and mean_diff < 0.01:
        print("\n⚠ Distributions are similar but have some differences")
        return False
    else:
        print("\n✗ Distributions differ significantly")
        return False

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: compare_with_salmon.py <salmon_exp_gc.gz> <our_expected_gc.tsv>")
        sys.exit(1)
    
    salmon_file = sys.argv[1]
    our_file = sys.argv[2]
    
    print("Parsing Salmon's expected GC...")
    gc_model = parse_salmon_gc_binary(salmon_file)
    print(f"  Format: {'LOG' if gc_model['dtype'] == 0 else 'LINEAR'}")
    print(f"  Dimensions: {gc_model['rows']} × {gc_model['cols']}")
    
    if gc_model['dtype'] == 0:
        # Convert from log space to linear
        matrix = np.array(gc_model['matrix'])
        matrix = np.exp(np.clip(matrix, -700, 700))  # Avoid overflow
        gc_model['matrix'] = matrix.tolist()
    
    salmon_dist = extract_marginal_distribution(gc_model)
    
    print("\nLoading our expected GC...")
    our_dist = load_our_expected_gc(our_file)
    
    match = compare_distributions(salmon_dist, our_dist)
    sys.exit(0 if match else 1)
