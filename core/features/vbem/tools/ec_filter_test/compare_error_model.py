#!/usr/bin/env python3
"""
Compare error model trace files from our CLI and Salmon.
Parses trace files and identifies divergences in errLike computation.
"""

import argparse
import sys
from collections import defaultdict
import re

def parse_trace_file(trace_file):
    """Parse trace file into structured data."""
    reads = defaultdict(lambda: {'alignments': [], 'read_info': None})
    
    if not trace_file:
        return reads
    
    try:
        with open(trace_file, 'r') as f:
            current_qname = None
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('READ '):
                    # READ <qname> numAlns=<N> errLikeSum=<X> modelUsed=<0/1> modelUpdated=<0/1>
                    parts = line.split()
                    qname = parts[1]
                    current_qname = qname
                    
                    read_info = {}
                    for part in parts[2:]:
                        if '=' in part:
                            key, value = part.split('=')
                            read_info[key] = value
                    
                    reads[qname]['read_info'] = read_info
                    reads[qname]['qname'] = qname
                
                elif line.startswith('ALN '):
                    # ALN <qname> tid=<X> fg=<Y> bg=<Z> errLike=<W> bin=<B>
                    if not current_qname:
                        continue
                    
                    parts = line.split()
                    aln_info = {'qname': parts[1]}
                    for part in parts[2:]:
                        if '=' in part:
                            key, value = part.split('=')
                            try:
                                aln_info[key] = float(value)
                            except ValueError:
                                aln_info[key] = value
                    
                    reads[current_qname]['alignments'].append(aln_info)
                
                elif line.startswith('TRANS '):
                    # TRANS <qname> tid=<X> readPos=<P> bin=<B> prevState=<S1> curState=<S2> logProb=<X>
                    if not current_qname:
                        continue
                    # For now, we'll focus on ALN level comparison
                    pass
    except FileNotFoundError:
        print(f"Warning: Trace file not found: {trace_file}", file=sys.stderr)
    except Exception as e:
        print(f"Error parsing trace file {trace_file}: {e}", file=sys.stderr)
    
    return reads

def compare_traces(ours_trace, salmon_trace, tolerance=1e-6):
    """Compare two trace files and report differences."""
    ours_reads = parse_trace_file(ours_trace)
    salmon_reads = parse_trace_file(salmon_trace)
    
    # Find common reads
    common_qnames = set(ours_reads.keys()) & set(salmon_reads.keys())
    ours_only = set(ours_reads.keys()) - set(salmon_reads.keys())
    salmon_only = set(salmon_reads.keys()) - set(ours_reads.keys())
    
    print(f"Total reads in ours: {len(ours_reads)}")
    print(f"Total reads in Salmon: {len(salmon_reads)}")
    print(f"Common reads: {len(common_qnames)}")
    print(f"Ours-only reads: {len(ours_only)}")
    print(f"Salmon-only reads: {len(salmon_only)}")
    print()
    
    # Compare common reads
    matching_reads = 0
    diverging_reads = []
    errLike_deltas = []
    
    for qname in sorted(common_qnames):
        ours_data = ours_reads[qname]
        salmon_data = salmon_reads[qname]
        
        # Compare number of alignments
        ours_alns = ours_data['alignments']
        salmon_alns = salmon_data['alignments']
        
        if len(ours_alns) != len(salmon_alns):
            diverging_reads.append({
                'qname': qname,
                'reason': f'Alignment count mismatch: ours={len(ours_alns)}, salmon={len(salmon_alns)}',
                'ours_alns': len(ours_alns),
                'salmon_alns': len(salmon_alns)
            })
            continue
        
        # Compare errLike values for each alignment
        read_matches = True
        for i, (ours_aln, salmon_aln) in enumerate(zip(ours_alns, salmon_alns)):
            ours_errLike = ours_aln.get('errLike', 0.0)
            salmon_errLike = salmon_aln.get('errLike', 0.0)
            
            delta = abs(ours_errLike - salmon_errLike)
            errLike_deltas.append(delta)
            
            if delta > tolerance:
                read_matches = False
                diverging_reads.append({
                    'qname': qname,
                    'reason': f'errLike mismatch at alignment {i}: ours={ours_errLike:.10f}, salmon={salmon_errLike:.10f}, delta={delta:.10f}',
                    'alignment_idx': i,
                    'ours_errLike': ours_errLike,
                    'salmon_errLike': salmon_errLike,
                    'delta': delta
                })
                break
        
        if read_matches:
            matching_reads += 1
    
    print(f"Matching reads (within tolerance {tolerance}): {matching_reads}/{len(common_qnames)}")
    print(f"Diverging reads: {len(diverging_reads)}")
    print()
    
    if errLike_deltas:
        import statistics
        print(f"errLike delta statistics:")
        print(f"  Mean: {statistics.mean(errLike_deltas):.10f}")
        print(f"  Median: {statistics.median(errLike_deltas):.10f}")
        print(f"  Max: {max(errLike_deltas):.10f}")
        print(f"  Min: {min(errLike_deltas):.10f}")
        print()
    
    # Report first divergence
    if diverging_reads:
        print("First divergence:")
        first = diverging_reads[0]
        print(f"  Read: {first['qname']}")
        print(f"  Reason: {first['reason']}")
        if 'alignment_idx' in first:
            print(f"  Alignment index: {first['alignment_idx']}")
            print(f"  Our errLike: {first['ours_errLike']:.10f}")
            print(f"  Salmon errLike: {first['salmon_errLike']:.10f}")
            print(f"  Delta: {first['delta']:.10f}")
        print()
        
        # Show first 10 divergences
        print("First 10 divergences:")
        for i, div in enumerate(diverging_reads[:10]):
            print(f"  {i+1}. {div['qname']}: {div['reason']}")
    
    return {
        'common_reads': len(common_qnames),
        'matching_reads': matching_reads,
        'diverging_reads': len(diverging_reads),
        'errLike_deltas': errLike_deltas
    }

def main():
    parser = argparse.ArgumentParser(description='Compare error model trace files')
    parser.add_argument('--ours', required=True, help='Our CLI trace file')
    parser.add_argument('--salmon', required=True, help='Salmon trace file')
    parser.add_argument('--tolerance', type=float, default=1e-6, help='Tolerance for errLike comparison')
    parser.add_argument('--report', help='Output report file')
    
    args = parser.parse_args()
    
    results = compare_traces(args.ours, args.salmon, args.tolerance)
    
    if args.report:
        with open(args.report, 'w') as f:
            f.write(f"Error Model Trace Comparison Report\n")
            f.write(f"{'='*50}\n\n")
            f.write(f"Common reads: {results['common_reads']}\n")
            f.write(f"Matching reads: {results['matching_reads']}\n")
            f.write(f"Diverging reads: {results['diverging_reads']}\n")
            if results['errLike_deltas']:
                import statistics
                f.write(f"\nErrLike delta statistics:\n")
                f.write(f"  Mean: {statistics.mean(results['errLike_deltas']):.10f}\n")
                f.write(f"  Median: {statistics.median(results['errLike_deltas']):.10f}\n")
                f.write(f"  Max: {max(results['errLike_deltas']):.10f}\n")
                f.write(f"  Min: {min(results['errLike_deltas']):.10f}\n")
        print(f"Report written to {args.report}")

if __name__ == '__main__':
    main()
