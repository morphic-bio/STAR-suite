#!/usr/bin/env python3
"""
Compare trace files from our CLI and Salmon to pinpoint divergence.

Usage:
    python3 compare_traces.py our_trace.txt salmon_trace.txt
"""

import sys
import re
from collections import defaultdict

def parse_trace_line(line):
    """Parse a structured trace line into a dictionary."""
    parts = line.strip().split('\t')
    if len(parts) < 2:
        return None
    
    qname = parts[0]
    data = {}
    
    # Parse semicolon-separated fields
    for field in parts[1].split(';'):
        if '=' not in field:
            continue
        key, value = field.split('=', 1)
        
        if key == 'txpIDs':
            data['txpIDs'] = [int(x) for x in value.split(',') if x]
        elif key == 'as':
            data['as'] = [int(x) for x in value.split(',') if x]
        elif key == 'bestAS':
            data['bestAS'] = int(value)
        elif key == 'logFragProb':
            data['logFragProb'] = [float(x) for x in value.split(',') if x]
        elif key == 'errLike':
            data['errLike'] = [float(x) for x in value.split(',') if x]
        elif key == 'logCompat':
            data['logCompat'] = [float(x) for x in value.split(',') if x]
        elif key == 'droppedIncompat':
            data['droppedIncompat'] = int(value)
        elif key == 'orphan':
            data['orphan'] = [int(x) for x in value.split(',') if x]
        elif key == 'auxProb':
            data['auxProb'] = [float(x) for x in value.split(',') if x]
        elif key == 'weight':
            data['weight'] = [float(x) for x in value.split(',') if x]
        elif key == 'ec_label':
            data['ec_label'] = value
    
    return {'qname': qname, **data}

def load_traces(filename):
    """Load trace file into a dictionary keyed by qname."""
    traces = {}
    with open(filename, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parsed = parse_trace_line(line)
            if parsed:
                traces[parsed['qname']] = parsed
    return traces

def compare_reads(our_read, salmon_read, tolerance=1e-6):
    """Compare two reads and return differences."""
    diffs = {
        'txpIDs_match': False,
        'weights_match': False,
        'component_diffs': {},
        'weight_deltas': []
    }
    
    # Check transcript IDs
    our_txpids = our_read.get('txpIDs', [])
    salmon_txpids = salmon_read.get('txpIDs', [])
    diffs['txpIDs_match'] = (our_txpids == salmon_txpids)
    
    # Check weights
    our_weights = our_read.get('weight', [])
    salmon_weights = salmon_read.get('weight', [])
    
    if len(our_weights) == len(salmon_weights):
        weight_match = True
        for i, (w1, w2) in enumerate(zip(our_weights, salmon_weights)):
            delta = abs(w1 - w2)
            diffs['weight_deltas'].append(delta)
            if delta > tolerance:
                weight_match = False
        diffs['weights_match'] = weight_match
    else:
        diffs['weights_match'] = False
    
    # Compare components
    our_errlike = our_read.get('errLike', [])
    salmon_errlike = salmon_read.get('errLike', [])
    if len(our_errlike) == len(salmon_errlike):
        errlike_diffs = [abs(a - b) for a, b in zip(our_errlike, salmon_errlike)]
        diffs['component_diffs']['errLike'] = max(errlike_diffs) if errlike_diffs else 0.0
    else:
        diffs['component_diffs']['errLike'] = float('inf')
    
    our_logfrag = our_read.get('logFragProb', [])
    salmon_logfrag = salmon_read.get('logFragProb', [])
    if len(our_logfrag) == len(salmon_logfrag):
        logfrag_diffs = [abs(a - b) for a, b in zip(our_logfrag, salmon_logfrag)]
        diffs['component_diffs']['logFragProb'] = max(logfrag_diffs) if logfrag_diffs else 0.0
    else:
        diffs['component_diffs']['logFragProb'] = float('inf')
    
    our_logcompat = our_read.get('logCompat', [])
    salmon_logcompat = salmon_read.get('logCompat', [])
    if len(our_logcompat) == len(salmon_logcompat):
        logcompat_diffs = [abs(a - b) for a, b in zip(our_logcompat, salmon_logcompat)]
        diffs['component_diffs']['logCompat'] = max(logcompat_diffs) if logcompat_diffs else 0.0
    else:
        diffs['component_diffs']['logCompat'] = float('inf')
    
    return diffs

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 compare_traces.py our_trace.txt salmon_trace.txt [tolerance]")
        sys.exit(1)
    
    our_file = sys.argv[1]
    salmon_file = sys.argv[2]
    tolerance = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6
    
    print("=" * 60)
    print("TRACE COMPARISON REPORT")
    print("=" * 60)
    print(f"Our trace: {our_file}")
    print(f"Salmon trace: {salmon_file}")
    print(f"Tolerance: {tolerance}")
    print()
    
    # Load traces
    print("Loading traces...")
    our_traces = load_traces(our_file)
    salmon_traces = load_traces(salmon_file)
    
    print(f"Loaded {len(our_traces)} reads from our trace")
    print(f"Loaded {len(salmon_traces)} reads from Salmon trace")
    print()
    
    # Find common reads
    common_qnames = set(our_traces.keys()) & set(salmon_traces.keys())
    our_only = set(our_traces.keys()) - set(salmon_traces.keys())
    salmon_only = set(salmon_traces.keys()) - set(our_traces.keys())
    
    print(f"Common reads: {len(common_qnames)}")
    print(f"Our-only reads: {len(our_only)}")
    print(f"Salmon-only reads: {len(salmon_only)}")
    print()
    
    # Compare common reads
    mismatches = []
    component_diff_counts = defaultdict(int)
    weight_deltas_all = []
    
    for qname in sorted(common_qnames):
        our_read = our_traces[qname]
        salmon_read = salmon_traces[qname]
        
        diffs = compare_reads(our_read, salmon_read, tolerance)
        
        if not diffs['txpIDs_match'] or not diffs['weights_match']:
            mismatches.append((qname, diffs))
            
            # Track which component differs most
            for comp, diff_val in diffs['component_diffs'].items():
                if diff_val > tolerance:
                    component_diff_counts[comp] += 1
        
        weight_deltas_all.extend(diffs['weight_deltas'])
    
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total common reads: {len(common_qnames)}")
    print(f"Mismatching reads: {len(mismatches)}")
    print(f"Match rate: {(len(common_qnames) - len(mismatches)) / len(common_qnames) * 100:.2f}%")
    print()
    
    if component_diff_counts:
        print("Component differences (reads with differences):")
        for comp, count in sorted(component_diff_counts.items(), key=lambda x: -x[1]):
            print(f"  {comp}: {count} reads")
        print()
    
    if weight_deltas_all:
        import statistics
        print("Weight delta statistics:")
        print(f"  Mean: {statistics.mean(weight_deltas_all):.6e}")
        print(f"  Median: {statistics.median(weight_deltas_all):.6e}")
        print(f"  Max: {max(weight_deltas_all):.6e}")
        print(f"  Min: {min(weight_deltas_all):.6e}")
        print()
    
    # Show first mismatch in detail
    if mismatches:
        print("=" * 60)
        print("FIRST MISMATCHING READ (detailed)")
        print("=" * 60)
        qname, diffs = mismatches[0]
        our_read = our_traces[qname]
        salmon_read = salmon_traces[qname]
        
        print(f"QNAME: {qname}")
        print()
        print("Transcript IDs:")
        print(f"  Ours:   {our_read.get('txpIDs', [])}")
        print(f"  Salmon: {salmon_read.get('txpIDs', [])}")
        print(f"  Match: {diffs['txpIDs_match']}")
        print()
        
        print("Weights:")
        print(f"  Ours:   {our_read.get('weight', [])}")
        print(f"  Salmon: {salmon_read.get('weight', [])}")
        print(f"  Match: {diffs['weights_match']}")
        if diffs['weight_deltas']:
            print(f"  Deltas: {diffs['weight_deltas']}")
        print()
        
        print("Components:")
        print(f"  errLike:")
        print(f"    Ours:   {our_read.get('errLike', [])}")
        print(f"    Salmon: {salmon_read.get('errLike', [])}")
        print(f"    Max diff: {diffs['component_diffs'].get('errLike', 0.0):.6e}")
        print()
        
        print(f"  logFragProb:")
        print(f"    Ours:   {our_read.get('logFragProb', [])}")
        print(f"    Salmon: {salmon_read.get('logFragProb', [])}")
        print(f"    Max diff: {diffs['component_diffs'].get('logFragProb', 0.0):.6e}")
        print()
        
        print(f"  logCompat:")
        print(f"    Ours:   {our_read.get('logCompat', [])}")
        print(f"    Salmon: {salmon_read.get('logCompat', [])}")
        print(f"    Max diff: {diffs['component_diffs'].get('logCompat', 0.0):.6e}")
        print()
        
        print(f"  AS tags:")
        print(f"    Ours:   {our_read.get('as', [])}")
        print(f"    Salmon: {salmon_read.get('as', [])}")
        print()
        
        print(f"  Best AS:")
        print(f"    Ours:   {our_read.get('bestAS', 'N/A')}")
        print(f"    Salmon: {salmon_read.get('bestAS', 'N/A')}")
        print()
        
        print(f"  Orphan flags:")
        print(f"    Ours:   {our_read.get('orphan', [])}")
        print(f"    Salmon: {salmon_read.get('orphan', [])}")
        print()
        
        print(f"  Dropped incompatible:")
        print(f"    Ours:   {our_read.get('droppedIncompat', 0)}")
        print(f"    Salmon: {salmon_read.get('droppedIncompat', 0)}")
        print()
        
        print(f"  EC label:")
        print(f"    Ours:   {our_read.get('ec_label', 'N/A')}")
        print(f"    Salmon: {salmon_read.get('ec_label', 'N/A')}")
        print()

if __name__ == '__main__':
    main()
