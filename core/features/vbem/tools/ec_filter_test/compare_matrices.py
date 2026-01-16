#!/usr/bin/env python3
"""
Compare transition matrix dumps from our CLI and Salmon.
Identifies cells with largest divergence.
"""

import argparse
import sys
import numpy as np

def parse_matrix_file(filename):
    """Parse TSV matrix file."""
    matrix = []
    header_info = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                # Parse header info
                if 'bin=' in line:
                    parts = line.split(',')
                    for part in parts:
                        if '=' in part:
                            key, value = part.split('=')
                            header_info[key.strip()] = value.strip()
                continue
            
            # Parse matrix row
            row = [float(x) for x in line.split('\t')]
            matrix.append(row)
    
    return np.array(matrix), header_info

def compare_matrices(ours_file, salmon_file, tolerance=1e-10):
    """Compare two matrix files and report differences."""
    ours_matrix, ours_info = parse_matrix_file(ours_file)
    salmon_matrix, salmon_info = parse_matrix_file(salmon_file)
    
    print(f"Our matrix shape: {ours_matrix.shape}")
    print(f"Salmon matrix shape: {salmon_matrix.shape}")
    print()
    
    if ours_matrix.shape != salmon_matrix.shape:
        print(f"ERROR: Matrix shapes don't match!")
        return None
    
    # Compute differences
    diff = ours_matrix - salmon_matrix
    abs_diff = np.abs(diff)
    
    print(f"Matrix comparison:")
    print(f"  Mean absolute difference: {np.mean(abs_diff):.10e}")
    print(f"  Max absolute difference: {np.max(abs_diff):.10e}")
    print(f"  Min absolute difference: {np.min(abs_diff):.10e}")
    print(f"  Cells within tolerance ({tolerance}): {np.sum(abs_diff < tolerance)}/{abs_diff.size}")
    print(f"  Cells exceeding tolerance: {np.sum(abs_diff >= tolerance)}")
    print()
    
    # Find cells with largest differences
    if np.max(abs_diff) > tolerance:
        print("Top 10 cells with largest differences:")
        flat_indices = np.argsort(abs_diff.flatten())[::-1]
        for i in range(min(10, len(flat_indices))):
            idx = flat_indices[i]
            row = idx // ours_matrix.shape[1]
            col = idx % ours_matrix.shape[1]
            ours_val = ours_matrix[row, col]
            salmon_val = salmon_matrix[row, col]
            diff_val = abs_diff[row, col]
            print(f"  [{row:3d}, {col:3d}]: ours={ours_val:15.10f}, salmon={salmon_val:15.10f}, diff={diff_val:.10e}")
        print()
    
    # Summary statistics
    matching_cells = np.sum(abs_diff < tolerance)
    total_cells = abs_diff.size
    match_percent = 100.0 * matching_cells / total_cells
    
    print(f"Summary:")
    print(f"  Matching cells: {matching_cells}/{total_cells} ({match_percent:.2f}%)")
    print(f"  Diverging cells: {total_cells - matching_cells}/{total_cells} ({100-match_percent:.2f}%)")
    
    return {
        'mean_diff': np.mean(abs_diff),
        'max_diff': np.max(abs_diff),
        'matching_cells': matching_cells,
        'total_cells': total_cells,
        'match_percent': match_percent
    }

def main():
    parser = argparse.ArgumentParser(description='Compare transition matrix dumps')
    parser.add_argument('--ours', required=True, help='Our CLI matrix file')
    parser.add_argument('--salmon', required=True, help='Salmon matrix file')
    parser.add_argument('--tolerance', type=float, default=1e-10, help='Tolerance for comparison')
    parser.add_argument('--output', help='Output comparison report')
    
    args = parser.parse_args()
    
    results = compare_matrices(args.ours, args.salmon, args.tolerance)
    
    if args.output and results:
        with open(args.output, 'w') as f:
            f.write("Transition Matrix Comparison Report\n")
            f.write("="*50 + "\n\n")
            f.write(f"Mean absolute difference: {results['mean_diff']:.10e}\n")
            f.write(f"Max absolute difference: {results['max_diff']:.10e}\n")
            f.write(f"Matching cells: {results['matching_cells']}/{results['total_cells']} ({results['match_percent']:.2f}%)\n")
        print(f"\nReport written to {args.output}")

if __name__ == '__main__':
    main()
