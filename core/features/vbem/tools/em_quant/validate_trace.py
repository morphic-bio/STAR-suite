#!/usr/bin/env python3
"""
Validate trace file format and check for inconsistencies.

Usage:
    validate_trace.py <trace_file>
"""

import sys
import argparse
from collections import defaultdict

def validate_trace(filename):
    """Validate trace file format."""
    errors = []
    warnings = []
    
    per_iter = defaultdict(dict)  # transcript -> {iter: values}
    ec_by_iter = defaultdict(list)  # (iter, transcript) -> [ecs]
    
    expected_header1 = "iter\ttranscript\talpha\tlogNorm\texpTheta\texpected_count"
    expected_header2 = "# EC-level trace: EC\titer\tec_id\ttranscript\tdenom\texpTheta\taux\tcontribution"
    
    header1_found = False
    header2_found = False
    
    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('#'):
                if expected_header2 in line:
                    header2_found = True
                continue
            
            if line.startswith('iter\t'):
                if line == expected_header1:
                    header1_found = True
                else:
                    errors.append(f"Line {line_num}: Unexpected header format: {line}")
                continue
            
            fields = line.split('\t')
            
            if fields[0] == 'EC':
                # EC-level trace
                if len(fields) < 8:
                    errors.append(f"Line {line_num}: EC line has {len(fields)} fields, expected 8")
                    continue
                
                try:
                    iter_num = int(fields[1])
                    ec_id = int(fields[2])
                    transcript = fields[3]
                    
                    if fields[4] != 'SINGLE':
                        denom = float(fields[4])
                        expTheta = float(fields[5])
                        aux = float(fields[6])
                        contribution = float(fields[7])
                    else:
                        denom = 'SINGLE'
                        expTheta = None
                        aux = None
                        contribution = float(fields[7])
                    
                    ec_by_iter[(iter_num, transcript)].append({
                        'ec_id': ec_id,
                        'denom': denom,
                        'expTheta': expTheta,
                        'aux': aux,
                        'contribution': contribution
                    })
                except (ValueError, IndexError) as e:
                    errors.append(f"Line {line_num}: Error parsing EC line: {e}")
            
            elif fields[0].isdigit() or fields[0] == '-1':
                # Per-iteration line
                if len(fields) < 6:
                    errors.append(f"Line {line_num}: Per-iteration line has {len(fields)} fields, expected 6")
                    continue
                
                try:
                    iter_num = int(fields[0])
                    transcript = fields[1]
                    alpha = float(fields[2])
                    logNorm = float(fields[3])
                    expTheta = float(fields[4])
                    expected_count = float(fields[5])
                    
                    per_iter[transcript][iter_num] = {
                        'alpha': alpha,
                        'logNorm': logNorm,
                        'expTheta': expTheta,
                        'expected_count': expected_count
                    }
                except (ValueError, IndexError) as e:
                    errors.append(f"Line {line_num}: Error parsing per-iteration line: {e}")
    
    # Validation checks
    if not header1_found:
        warnings.append("Per-iteration header not found")
    if not header2_found:
        warnings.append("EC-level header comment not found")
    
    # Check for consistency
    for transcript, iters in per_iter.items():
        if -1 not in iters:
            warnings.append(f"Transcript {transcript}: No initial state (iter -1) found")
        
        # Check if logNorm is consistent across transcripts in same iteration
        logNorms_by_iter = defaultdict(set)
        for iter_num, vals in iters.items():
            logNorms_by_iter[iter_num].add(vals['logNorm'])
        
        for iter_num, logNorms in logNorms_by_iter.items():
            if len(logNorms) > 1:
                errors.append(f"Transcript {transcript}, iter {iter_num}: Inconsistent logNorm values: {logNorms}")
    
    # Check EC contributions sum to expected_count
    for (iter_num, transcript), ecs in ec_by_iter.items():
        if iter_num == -1:
            continue
        
        if transcript in per_iter and iter_num in per_iter[transcript]:
            expected = per_iter[transcript][iter_num]['expected_count']
            total_contrib = sum(ec['contribution'] for ec in ecs if ec['contribution'] is not None)
            
            if abs(total_contrib - expected) > 1e-6:
                warnings.append(f"Transcript {transcript}, iter {iter_num}: "
                              f"EC contributions sum to {total_contrib:.10f}, "
                              f"but expected_count is {expected:.10f} "
                              f"(diff: {abs(total_contrib - expected):.2e})")
    
    # Print results
    print(f"Trace file: {filename}")
    print(f"  Transcripts: {len(per_iter)}")
    print(f"  Total per-iteration entries: {sum(len(iters) for iters in per_iter.values())}")
    print(f"  Total EC entries: {sum(len(ecs) for ecs in ec_by_iter.values())}")
    print()
    
    if warnings:
        print("WARNINGS:")
        for w in warnings:
            print(f"  - {w}")
        print()
    
    if errors:
        print("ERRORS:")
        for e in errors:
            print(f"  - {e}")
        print()
        return False
    else:
        print("✓ Trace file format is valid")
        return True

def main():
    parser = argparse.ArgumentParser(description='Validate trace file format')
    parser.add_argument('trace_file', help='Trace file to validate')
    
    args = parser.parse_args()
    
    if not validate_trace(args.trace_file):
        sys.exit(1)

if __name__ == '__main__':
    main()
