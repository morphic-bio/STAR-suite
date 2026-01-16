#!/usr/bin/env python3
"""
Compare debug traces between em_quant and Salmon to find first divergence.

Usage:
    compare_traces.py <em_trace.txt> <salmon_trace.txt> [--transcript <txp_id>]

Compares per-iteration and per-EC values to identify where divergence occurs.
"""

import sys
import argparse
from collections import defaultdict

def parse_trace(filename):
    """Parse trace file into structured data."""
    per_iter = {}  # (iter, transcript) -> {alpha, logNorm, expTheta, expected_count}
    ec_details = []  # List of EC entries
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if fields[0] == 'EC':
                # EC-level trace: EC	iter	ec_id	transcript	denom	expTheta	aux	contribution
                if len(fields) >= 8:
                    ec_details.append({
                        'iter': int(fields[1]),
                        'ec_id': int(fields[2]),
                        'transcript': fields[3],
                        'denom': float(fields[4]) if fields[4] != 'SINGLE' else 'SINGLE',
                        'expTheta': float(fields[5]) if fields[5] != 'SINGLE' else None,
                        'aux': float(fields[6]) if fields[6] != 'SINGLE' else None,
                        'contribution': float(fields[7]) if fields[7] != 'SINGLE' else float(fields[7])
                    })
            elif fields[0].isdigit() or fields[0] == '-1':
                # Per-iteration: iter	transcript	alpha	logNorm	expTheta	expected_count
                if len(fields) >= 6:
                    iter_num = int(fields[0])
                    transcript = fields[1]
                    per_iter[(iter_num, transcript)] = {
                        'alpha': float(fields[2]),
                        'logNorm': float(fields[3]),
                        'expTheta': float(fields[4]),
                        'expected_count': float(fields[5])
                    }
    
    return per_iter, ec_details

def compare_per_iter(em_per_iter, salmon_per_iter, transcript=None, tolerance=1e-10):
    """Compare per-iteration values."""
    print("=" * 80)
    print("Per-Iteration Comparison")
    print("=" * 80)
    
    # Find common (iter, transcript) pairs
    em_keys = set(em_per_iter.keys())
    salmon_keys = set(salmon_per_iter.keys())
    common_keys = em_keys & salmon_keys
    
    if transcript:
        common_keys = {k for k in common_keys if k[1] == transcript}
    
    if not common_keys:
        print("No common iterations/transcripts found!")
        return
    
    # Sort by iteration, then transcript
    sorted_keys = sorted(common_keys)
    
    first_diff = None
    diff_count = 0
    
    for key in sorted_keys:
        iter_num, txp = key
        em_vals = em_per_iter[key]
        salmon_vals = salmon_per_iter[key]
        
        diffs = {}
        for field in ['alpha', 'logNorm', 'expTheta', 'expected_count']:
            em_val = em_vals[field]
            salmon_val = salmon_vals[field]
            diff = abs(em_val - salmon_val)
            rel_diff = diff / abs(salmon_val) if salmon_val != 0 else diff
            diffs[field] = (diff, rel_diff)
            
            if diff > tolerance:
                if first_diff is None:
                    first_diff = (iter_num, txp, field, em_val, salmon_val, diff, rel_diff)
                diff_count += 1
        
        # Print if any differences found
        if any(d[0] > tolerance for d in diffs.values()):
            print(f"\nIteration {iter_num}, Transcript {txp}:")
            for field, (abs_diff, rel_diff) in diffs.items():
                if abs_diff > tolerance:
                    print(f"  {field:15s}: em={em_vals[field]:15.10f}, salmon={salmon_vals[field]:15.10f}, "
                          f"diff={abs_diff:.2e} ({rel_diff*100:.6f}%)")
    
    if first_diff:
        iter_num, txp, field, em_val, salmon_val, diff, rel_diff = first_diff
        print(f"\n{'='*80}")
        print(f"FIRST DIVERGENCE:")
        print(f"  Iteration: {iter_num}")
        print(f"  Transcript: {txp}")
        print(f"  Field: {field}")
        print(f"  em_quant:   {em_val:.10f}")
        print(f"  Salmon:     {salmon_val:.10f}")
        print(f"  Difference: {diff:.2e} ({rel_diff*100:.6f}%)")
        print(f"{'='*80}")
    else:
        print(f"\n✓ All per-iteration values match within tolerance ({tolerance:.2e})")
    
    print(f"\nTotal differences found: {diff_count}")

def compare_ec_details(em_ecs, salmon_ecs, transcript=None, iter_num=None, tolerance=1e-10):
    """Compare EC-level details."""
    print("\n" + "=" * 80)
    print("EC-Level Comparison")
    print("=" * 80)
    
    # Filter by transcript and iteration if specified
    em_filtered = [ec for ec in em_ecs 
                   if (transcript is None or ec['transcript'] == transcript) and
                      (iter_num is None or ec['iter'] == iter_num)]
    salmon_filtered = [ec for ec in salmon_ecs 
                       if (transcript is None or ec['transcript'] == transcript) and
                          (iter_num is None or ec['iter'] == iter_num)]
    
    # Group by (iter, ec_id, transcript)
    em_by_key = {(ec['iter'], ec['ec_id'], ec['transcript']): ec for ec in em_filtered}
    salmon_by_key = {(ec['iter'], ec['ec_id'], ec['transcript']): ec for ec in salmon_filtered}
    
    common_keys = em_by_key.keys() & salmon_by_key.keys()
    
    if not common_keys:
        print("No common ECs found!")
        return
    
    first_diff = None
    diff_count = 0
    
    for key in sorted(common_keys):
        iter_num, ec_id, txp = key
        em_ec = em_by_key[key]
        salmon_ec = salmon_by_key[key]
        
        # Skip single-transcript ECs for now (they should match exactly)
        if em_ec['denom'] == 'SINGLE':
            continue
        
        diffs = {}
        for field in ['denom', 'expTheta', 'aux', 'contribution']:
            if field in em_ec and field in salmon_ec:
                em_val = em_ec[field]
                salmon_val = salmon_ec[field]
                if em_val is None or salmon_val is None:
                    continue
                diff = abs(em_val - salmon_val)
                rel_diff = diff / abs(salmon_val) if salmon_val != 0 else diff
                diffs[field] = (diff, rel_diff)
                
                if diff > tolerance:
                    if first_diff is None:
                        first_diff = (iter_num, ec_id, txp, field, em_val, salmon_val, diff, rel_diff)
                    diff_count += 1
        
        # Print if any differences found
        if any(d[0] > tolerance for d in diffs.values()):
            print(f"\nIteration {iter_num}, EC {ec_id}, Transcript {txp}:")
            for field, (abs_diff, rel_diff) in diffs.items():
                if abs_diff > tolerance:
                    print(f"  {field:15s}: em={em_ec[field]:15.10f}, salmon={salmon_ec[field]:15.10f}, "
                          f"diff={abs_diff:.2e} ({rel_diff*100:.6f}%)")
    
    if first_diff:
        iter_num, ec_id, txp, field, em_val, salmon_val, diff, rel_diff = first_diff
        print(f"\n{'='*80}")
        print(f"FIRST EC DIVERGENCE:")
        print(f"  Iteration: {iter_num}")
        print(f"  EC ID: {ec_id}")
        print(f"  Transcript: {txp}")
        print(f"  Field: {field}")
        print(f"  em_quant:   {em_val:.10f}")
        print(f"  Salmon:     {salmon_val:.10f}")
        print(f"  Difference: {diff:.2e} ({rel_diff*100:.6f}%)")
        print(f"{'='*80}")
    else:
        print(f"\n✓ All EC-level values match within tolerance ({tolerance:.2e})")
    
    print(f"\nTotal EC differences found: {diff_count}")

def main():
    parser = argparse.ArgumentParser(description='Compare em_quant and Salmon debug traces')
    parser.add_argument('em_trace', help='em_quant trace file')
    parser.add_argument('salmon_trace', help='Salmon trace file')
    parser.add_argument('--transcript', help='Focus on specific transcript')
    parser.add_argument('--iter', type=int, help='Focus on specific iteration')
    parser.add_argument('--tolerance', type=float, default=1e-10, help='Tolerance for comparisons')
    
    args = parser.parse_args()
    
    print(f"Loading em_quant trace: {args.em_trace}")
    em_per_iter, em_ecs = parse_trace(args.em_trace)
    print(f"  Found {len(em_per_iter)} per-iteration entries")
    print(f"  Found {len(em_ecs)} EC-level entries")
    
    print(f"\nLoading Salmon trace: {args.salmon_trace}")
    salmon_per_iter, salmon_ecs = parse_trace(args.salmon_trace)
    print(f"  Found {len(salmon_per_iter)} per-iteration entries")
    print(f"  Found {len(salmon_ecs)} EC-level entries")
    
    # Compare per-iteration values
    compare_per_iter(em_per_iter, salmon_per_iter, args.transcript, args.tolerance)
    
    # Compare EC-level details
    compare_ec_details(em_ecs, salmon_ecs, args.transcript, args.iter, args.tolerance)

if __name__ == '__main__':
    main()
