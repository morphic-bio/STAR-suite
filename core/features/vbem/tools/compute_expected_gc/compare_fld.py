#!/usr/bin/env python3
"""
Compare fragment length distributions between Salmon and our sample_fld tool.
"""

import sys
import gzip
import struct
import numpy as np

def read_salmon_fld(filename):
    """Read Salmon's binary FLD format.
    
    Salmon's FLD format stores raw counts as 4-byte integers, one per length.
    File size / 4 = number of length bins.
    """
    with gzip.open(filename, 'rb') as f:
        data = f.read()
    
    # FLD is stored as array of 4-byte unsigned integers (counts)
    num_entries = len(data) // 4
    counts = []
    for i in range(num_entries):
        val = struct.unpack('<I', data[i*4:(i+1)*4])[0]
        counts.append(val)
    
    counts = np.array(counts, dtype=float)
    
    # Normalize to PMF
    total = counts.sum()
    if total > 0:
        pmf = counts / total
    else:
        pmf = counts
    
    return pmf, num_entries - 1

def read_our_fld(filename):
    """Read our TSV FLD format."""
    fld = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                try:
                    length = int(parts[0])
                    # Use probability if available (3rd column), otherwise count
                    if len(parts) >= 3:
                        prob = float(parts[2])
                    else:
                        count = int(parts[1])
                        prob = float(count)  # Will normalize later
                    fld[length] = prob
                except ValueError:
                    continue  # Skip header lines
    
    # Find max length
    max_len = max(fld.keys()) if fld else 0
    
    # Convert to array (use size to match Salmon's 1001 bins)
    array_size = max(max_len + 1, 1001)
    pmf = np.zeros(array_size)
    for length, prob in fld.items():
        if length < len(pmf):
            pmf[length] = prob
    
    # Normalize if needed (probabilities should already sum to 1)
    total = pmf.sum()
    if total > 0 and abs(total - 1.0) > 0.01:  # Only normalize if not already normalized
        pmf /= total
    
    return pmf, max_len

def compare_fld(salmon_file, our_file, tolerance=1e-6):
    """Compare two FLD files."""
    print(f"Reading Salmon FLD: {salmon_file}")
    salmon_pmf, salmon_max = read_salmon_fld(salmon_file)
    
    print(f"Reading our FLD: {our_file}")
    our_pmf, our_max = read_our_fld(our_file)
    
    # Align lengths
    max_len = min(len(salmon_pmf), len(our_pmf))
    salmon_pmf = salmon_pmf[:max_len]
    our_pmf = our_pmf[:max_len]
    
    # Compare
    diff = np.abs(salmon_pmf - our_pmf)
    max_diff = np.max(diff)
    mean_diff = np.mean(diff)
    
    # Find non-zero differences
    nonzero_diff = diff[diff > tolerance]
    
    print(f"\nComparison Results:")
    print(f"  Max length: Salmon={len(salmon_pmf)-1}, Ours={len(our_pmf)-1}")
    print(f"  Max absolute difference: {max_diff:.10f}")
    print(f"  Mean absolute difference: {mean_diff:.10f}")
    print(f"  Non-zero differences (> {tolerance}): {len(nonzero_diff)}")
    
    if len(nonzero_diff) > 0:
        print(f"  Largest differences:")
        idx = np.argsort(diff)[-10:][::-1]
        for i in idx:
            if diff[i] > tolerance:
                print(f"    Length {i}: Salmon={salmon_pmf[i]:.10f}, Ours={our_pmf[i]:.10f}, Diff={diff[i]:.10f}")
    
    # Determine pass/fail
    if max_diff < tolerance:
        print(f"\n✅ PASS: Distributions match within tolerance ({tolerance})")
        return 0
    else:
        print(f"\n⚠️  WARN: Distributions differ (max diff: {max_diff:.10f} > {tolerance})")
        return 1

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: compare_fld.py <salmon_fld.gz> <our_fld.tsv> [tolerance]")
        sys.exit(1)
    
    salmon_file = sys.argv[1]
    our_file = sys.argv[2]
    tolerance = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
    
    sys.exit(compare_fld(salmon_file, our_file, tolerance))
