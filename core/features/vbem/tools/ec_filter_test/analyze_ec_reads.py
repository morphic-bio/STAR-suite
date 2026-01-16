#!/usr/bin/env python3
"""
Analyze EC-level read differences between Salmon and our CLI.

Parses trace files from both tools, groups reads by EC label, and identifies
which reads cause count differences.

Usage:
    python3 analyze_ec_reads.py --salmon-trace salmon_trace.txt --our-trace our_trace.txt --ec-file salmon_eq_classes.txt --our-ec-file our_eq_classes.txt --output report.txt
"""

import argparse
import gzip
import sys
from collections import defaultdict
from typing import Dict, Set, List, Tuple

def parse_trace_file(trace_file: str) -> Dict[str, List[str]]:
    """
    Parse trace file and group reads by EC label.
    
    Returns: dict mapping EC label -> list of read names
    """
    ec_to_reads = defaultdict(list)
    
    open_func = gzip.open if trace_file.endswith('.gz') else open
    mode = 'rt' if trace_file.endswith('.gz') else 'r'
    
    with open_func(trace_file, mode) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            
            qname = parts[0]
            
            # Extract EC label from fields
            ec_label = None
            for field in parts[1].split(';'):
                if field.startswith('ec_label='):
                    ec_label = field.split('=', 1)[1]
                    break
            
            if ec_label:
                ec_to_reads[ec_label].append(qname)
            else:
                # If no EC label, this read was filtered (zero-prob fragment)
                ec_to_reads['__FILTERED__'].append(qname)
    
    return ec_to_reads

def parse_ec_file(filepath: str) -> Dict[str, float]:
    """
    Parse EC file and extract EC labels and counts.
    
    Returns: dict mapping EC label (as tuple) -> count
    """
    open_func = gzip.open if filepath.endswith('.gz') else open
    mode = 'rt' if filepath.endswith('.gz') else 'r'
    
    with open_func(filepath, mode) as f:
        lines = [line.strip() for line in f if line.strip()]
    
    num_transcripts = int(lines[0])
    transcript_names = lines[2:2+num_transcripts]
    ec_lines = lines[2+num_transcripts:]
    
    ec_counts = {}
    for line in ec_lines:
        parts = line.split()
        if len(parts) < 2:
            continue
        
        k = int(parts[0])
        transcript_ids = [int(parts[i]) for i in range(1, k+1)]
        
        # Determine if weighted format
        if len(parts) == 2 * k + 2:
            count = float(parts[-1])
        elif len(parts) == k + 2:
            count = float(parts[-1])
        else:
            continue
        
        # Create EC label as comma-separated string (matching trace format)
        ec_label = ','.join(str(tid) for tid in transcript_ids)
        ec_counts[ec_label] = count
    
    return ec_counts

def normalize_ec_label(ec_label: str) -> str:
    """
    Normalize EC label for comparison (sort transcript IDs).
    """
    if ec_label == '__FILTERED__':
        return ec_label
    
    txp_ids = [int(x) for x in ec_label.split(',') if x]
    txp_ids.sort()
    return ','.join(str(tid) for tid in txp_ids)

def compare_ec_reads(
    salmon_trace: str,
    our_trace: str,
    salmon_ec_file: str,
    our_ec_file: str,
    output_file: str
):
    """
    Compare EC-to-reads mappings and identify differences.
    """
    print("Parsing trace files...")
    
    # Try to parse Salmon trace if available, otherwise use empty dict
    salmon_ec_reads = {}
    if salmon_trace and salmon_trace.strip():
        try:
            salmon_ec_reads = parse_trace_file(salmon_trace)
            print(f"  Parsed Salmon trace: {len(salmon_ec_reads)} ECs")
        except FileNotFoundError:
            print(f"  Warning: Salmon trace file not found: {salmon_trace}")
            print("  Will analyze based on EC file counts only")
    else:
        print("  No Salmon trace file provided, using EC counts only")
    
    our_ec_reads = parse_trace_file(our_trace)
    print(f"  Parsed our trace: {len(our_ec_reads)} ECs")
    
    print("Parsing EC files...")
    salmon_ec_counts = parse_ec_file(salmon_ec_file)
    our_ec_counts = parse_ec_file(our_ec_file)
    
    # Normalize EC labels (sort transcript IDs)
    salmon_ec_reads_norm = {}
    if salmon_ec_reads:
        for ec_label, reads in salmon_ec_reads.items():
            norm_label = normalize_ec_label(ec_label)
            if norm_label not in salmon_ec_reads_norm:
                salmon_ec_reads_norm[norm_label] = []
            salmon_ec_reads_norm[norm_label].extend(reads)
    
    our_ec_reads_norm = {}
    for ec_label, reads in our_ec_reads.items():
        norm_label = normalize_ec_label(ec_label)
        if norm_label not in our_ec_reads_norm:
            our_ec_reads_norm[norm_label] = []
        our_ec_reads_norm[norm_label].extend(reads)
    
    # Normalize EC counts labels too
    salmon_ec_counts_norm = {}
    for ec_label, count in salmon_ec_counts.items():
        norm_label = normalize_ec_label(ec_label)
        salmon_ec_counts_norm[norm_label] = count
    
    our_ec_counts_norm = {}
    for ec_label, count in our_ec_counts.items():
        norm_label = normalize_ec_label(ec_label)
        our_ec_counts_norm[norm_label] = count
    
    # Find ECs with count differences
    all_ec_labels = set(salmon_ec_counts_norm.keys()) | set(our_ec_counts_norm.keys())
    
    differing_ecs = []
    for ec_label in all_ec_labels:
        salmon_count = salmon_ec_counts_norm.get(ec_label, 0.0)
        our_count = our_ec_counts_norm.get(ec_label, 0.0)
        
        if abs(salmon_count - our_count) > 0.5:  # Allow small floating point differences
            differing_ecs.append((ec_label, salmon_count, our_count))
    
    differing_ecs.sort(key=lambda x: abs(x[2] - x[1]), reverse=True)
    
    # Write report
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("EC-LEVEL READ DIFFERENCE ANALYSIS\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total ECs with count differences: {len(differing_ecs)}\n")
        f.write(f"Salmon total reads: {sum(salmon_ec_counts_norm.values())}\n")
        f.write(f"Our total reads: {sum(our_ec_counts_norm.values())}\n")
        f.write(f"Difference: {sum(our_ec_counts_norm.values()) - sum(salmon_ec_counts_norm.values())}\n\n")
        
        # Analyze filtered reads
        salmon_filtered = salmon_ec_reads_norm.get('__FILTERED__', [])
        our_filtered = our_ec_reads_norm.get('__FILTERED__', [])
        
        f.write(f"Filtered reads (zero-prob fragments):\n")
        f.write(f"  Salmon: {len(salmon_filtered)}\n")
        f.write(f"  Ours: {len(our_filtered)}\n")
        f.write(f"  Difference: {len(our_filtered) - len(salmon_filtered)}\n\n")
        
        # For each differing EC, show read differences
        f.write("=" * 80 + "\n")
        f.write("DETAILED EC DIFFERENCES\n")
        f.write("=" * 80 + "\n\n")
        
        total_ours_only = 0
        total_salmon_only = 0
        
        for ec_label, salmon_count, our_count in differing_ecs[:20]:  # Top 20 differences
            if ec_label == '__FILTERED__':
                continue
            
            our_reads = set(our_ec_reads_norm.get(ec_label, []))
            
            if salmon_ec_reads_norm:
                salmon_reads = set(salmon_ec_reads_norm.get(ec_label, []))
                ours_only = our_reads - salmon_reads
                salmon_only = salmon_reads - our_reads
            else:
                # No Salmon trace, so we can't compare read-by-read
                # Just show our reads and the count difference
                salmon_reads = set()
                ours_only = our_reads  # All our reads are "ours-only" if we can't compare
                salmon_only = set()
            
            total_ours_only += len(ours_only)
            total_salmon_only += len(salmon_only)
            
            diff = our_count - salmon_count
            
            f.write(f"EC ({ec_label}):\n")
            f.write(f"  Salmon count: {salmon_count:.0f}\n")
            f.write(f"  Our count: {our_count:.0f}\n")
            f.write(f"  Difference: {diff:.0f}\n")
            f.write(f"  Reads in Salmon but not ours: {len(salmon_only)}\n")
            f.write(f"  Reads in ours but not Salmon: {len(ours_only)}\n")
            
            if ours_only:
                f.write(f"  Ours-only reads (sample of {min(10, len(ours_only))}):\n")
                for qname in sorted(ours_only)[:10]:
                    f.write(f"    {qname}\n")
                if len(ours_only) > 10:
                    f.write(f"    ... and {len(ours_only) - 10} more\n")
            
            if salmon_only:
                f.write(f"  Salmon-only reads (sample of {min(10, len(salmon_only))}):\n")
                for qname in sorted(salmon_only)[:10]:
                    f.write(f"    {qname}\n")
                if len(salmon_only) > 10:
                    f.write(f"    ... and {len(salmon_only) - 10} more\n")
            
            f.write("\n")
        
        f.write("=" * 80 + "\n")
        f.write("SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Total reads in ours but not Salmon: {total_ours_only}\n")
        f.write(f"Total reads in Salmon but not ours: {total_salmon_only}\n")
        f.write(f"Net difference: {total_ours_only - total_salmon_only}\n")
    
    print(f"\nAnalysis complete. Report written to: {output_file}")
    print(f"Found {len(differing_ecs)} ECs with count differences")
    print(f"Total ours-only reads: {total_ours_only}")
    print(f"Total salmon-only reads: {total_salmon_only}")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze EC-level read differences between Salmon and our CLI'
    )
    parser.add_argument('--salmon-trace', required=False, default=None,
                       help='Salmon trace file (from SALMON_TRACE_FILE, optional)')
    parser.add_argument('--our-trace', required=True,
                       help='Our CLI trace file (from --trace-reads)')
    parser.add_argument('--salmon-ec-file', required=True,
                       help='Salmon eq_classes.txt file')
    parser.add_argument('--our-ec-file', required=True,
                       help='Our eq_classes.txt file')
    parser.add_argument('--output', required=True,
                       help='Output report file')
    
    args = parser.parse_args()
    
    # Handle empty string as None
    salmon_trace = args.salmon_trace if args.salmon_trace and args.salmon_trace.strip() else None
    
    compare_ec_reads(
        salmon_trace,
        args.our_trace,
        args.salmon_ec_file,
        args.our_ec_file,
        args.output
    )

if __name__ == '__main__':
    main()
