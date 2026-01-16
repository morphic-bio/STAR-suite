#!/usr/bin/env python3
"""
Identify reads that likely differ between Salmon and our CLI.

Since we can't get Salmon's trace file, we analyze our reads to identify
potential patterns that might cause filtering differences.
"""

import argparse
import subprocess
import sys
from collections import defaultdict

def parse_our_trace(trace_file: str) -> dict:
    """Parse our trace file and extract read information."""
    reads = {}
    
    with open(trace_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            
            qname = parts[0]
            fields = {}
            for field in parts[1].split(';'):
                if '=' in field:
                    k, v = field.split('=', 1)
                    fields[k] = v
            
            reads[qname] = {
                'ec_label': fields.get('ec_label', ''),
                'txpIDs': fields.get('txpIDs', '').split(',') if fields.get('txpIDs') else [],
                'as': fields.get('as', '').split(',') if fields.get('as') else [],
                'orphan': fields.get('orphan', '0') == '1',
                'dropped_incompat': fields.get('droppedIncompat', '0') == '1',
            }
    
    return reads

def get_bam_characteristics(bam_file: str, read_names: list) -> dict:
    """Get BAM characteristics for a set of read names."""
    read_chars = defaultdict(lambda: {
        'flags': [],
        'mapqs': [],
        'nhs': [],
        'has_as': False,
        'as_values': [],
    })
    
    cmd = ['samtools', 'view', bam_file]
    read_set = set(read_names)
    
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in proc.stdout:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            
            qname = parts[0]
            if qname not in read_set:
                continue
            
            flag = int(parts[1])
            mapq = int(parts[4])
            
            # Extract tags
            nh = 0
            as_val = None
            for tag_field in parts[11:]:
                if ':' in tag_field:
                    tag_parts = tag_field.split(':')
                    if len(tag_parts) >= 3:
                        tag = tag_parts[0]
                        tag_type = tag_parts[1]
                        tag_value = tag_parts[2]
                        if tag == 'NH' and tag_type == 'i':
                            nh = int(tag_value)
                        elif tag == 'AS' and tag_type == 'i':
                            as_val = int(tag_value)
            
            read_chars[qname]['flags'].append(flag)
            read_chars[qname]['mapqs'].append(mapq)
            read_chars[qname]['nhs'].append(nh)
            if as_val is not None:
                read_chars[qname]['has_as'] = True
                read_chars[qname]['as_values'].append(as_val)
        
        proc.terminate()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
    
    return read_chars

def analyze_ec_differences(analysis_file: str, our_trace: str, bam_file: str, output_file: str):
    """Analyze reads from ECs with differences."""
    
    # Parse our trace
    print("Parsing our trace file...")
    our_reads = parse_our_trace(our_trace)
    print(f"  Found {len(our_reads)} reads")
    
    # Extract ECs with differences from analysis file
    differing_ecs = []
    current_ec = None
    ec_reads = defaultdict(list)
    
    with open(analysis_file, 'r') as f:
        for line in f:
            if line.startswith('EC ('):
                ec_label = line.split('(')[1].split(')')[0]
                current_ec = ec_label
            elif 'Difference:' in line and current_ec:
                try:
                    diff = float(line.split('Difference:')[1].strip().split()[0])
                    if abs(diff) > 0.5:
                        differing_ecs.append((current_ec, diff))
                except:
                    pass
            elif current_ec and line.strip().startswith('SRR'):
                read_name = line.strip().split()[0]
                ec_reads[current_ec].append(read_name)
    
    print(f"Found {len(differing_ecs)} ECs with differences")
    
    # Analyze top differing ECs
    differing_ecs.sort(key=lambda x: abs(x[1]), reverse=True)
    
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("READ CHARACTERISTICS ANALYSIS FOR DIFFERING ECs\n")
        f.write("=" * 80 + "\n\n")
        
        for ec_label, diff in differing_ecs[:10]:  # Top 10
            f.write(f"\nEC ({ec_label}): Difference = {diff:.0f}\n")
            f.write("-" * 80 + "\n")
            
            # Get reads for this EC
            sample_reads = ec_reads.get(ec_label, [])[:50]  # Sample 50 reads
            
            if not sample_reads:
                f.write("  No sample reads found in analysis file\n")
                continue
            
            # Get BAM characteristics
            bam_chars = get_bam_characteristics(bam_file, sample_reads)
            
            # Analyze characteristics
            total_flags = []
            total_mapqs = []
            total_nhs = []
            reads_with_as = 0
            reads_without_as = 0
            proper_pairs = 0
            orphans = 0
            
            for read_name in sample_reads:
                if read_name not in bam_chars:
                    continue
                
                chars = bam_chars[read_name]
                total_flags.extend(chars['flags'])
                total_mapqs.extend(chars['mapqs'])
                total_nhs.extend(chars['nhs'])
                
                if chars['has_as']:
                    reads_with_as += 1
                else:
                    reads_without_as += 1
                
                # Check if proper pair
                if chars['flags']:
                    flag = chars['flags'][0]
                    if flag & 0x1 and flag & 0x2:  # Paired and proper
                        proper_pairs += 1
                    elif flag & 0x1 and not (flag & 0x2):  # Paired but not proper
                        orphans += 1
            
            f.write(f"  Sample size: {len(sample_reads)} reads\n")
            if total_mapqs:
                f.write(f"  Average MAPQ: {sum(total_mapqs)/len(total_mapqs):.1f}\n")
            if total_nhs:
                f.write(f"  Average NH: {sum(total_nhs)/len(total_nhs):.1f}\n")
            f.write(f"  Reads with AS tag: {reads_with_as}/{len(sample_reads)}\n")
            f.write(f"  Reads without AS tag: {reads_without_as}/{len(sample_reads)}\n")
            f.write(f"  Proper pairs: {proper_pairs}\n")
            f.write(f"  Orphans: {orphans}\n")
            
            # Flag distribution
            flag_counts = defaultdict(int)
            for flag in total_flags:
                flag_counts[flag] += 1
            f.write(f"  Flag distribution: {dict(flag_counts)}\n")
    
    print(f"Analysis written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Identify characteristics of reads in differing ECs'
    )
    parser.add_argument('--analysis-file', required=True,
                       help='Output from analyze_ec_reads.py')
    parser.add_argument('--our-trace', required=True,
                       help='Our CLI trace file')
    parser.add_argument('--bam-file', required=True,
                       help='BAM file')
    parser.add_argument('--output', required=True,
                       help='Output file')
    
    args = parser.parse_args()
    
    analyze_ec_differences(
        args.analysis_file,
        args.our_trace,
        args.bam_file,
        args.output
    )

if __name__ == '__main__':
    main()
