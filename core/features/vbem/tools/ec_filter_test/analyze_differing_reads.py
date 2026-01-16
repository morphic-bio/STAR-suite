#!/usr/bin/env python3
"""
Analyze characteristics of reads that differ between Salmon and our CLI.

Extracts BAM characteristics for reads in ECs with count differences.
"""

import argparse
import subprocess
import sys
from collections import defaultdict

def extract_reads_from_analysis(analysis_file: str, max_reads: int = 50) -> dict:
    """Extract read names from the analysis file for ECs with differences."""
    ec_to_reads = defaultdict(list)
    current_ec = None
    
    with open(analysis_file, 'r') as f:
        for line in f:
            if line.startswith('EC ('):
                # Extract EC label
                ec_label = line.split('(')[1].split(')')[0]
                current_ec = ec_label
            elif current_ec and 'Ours-only reads' in line:
                # Skip header line
                continue
            elif current_ec and line.strip().startswith('SRR'):
                # Read name
                read_name = line.strip().split()[0]
                ec_to_reads[current_ec].append(read_name)
                if len(ec_to_reads[current_ec]) >= max_reads:
                    break
    
    return ec_to_reads

def analyze_bam_reads(bam_file: str, read_names: list) -> list:
    """Extract BAM characteristics for given read names."""
    results = []
    
    # Use samtools to extract reads
    read_set = set(read_names)
    cmd = ['samtools', 'view', bam_file]
    
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        for line in proc.stdout:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            
            qname = parts[0]
            if qname in read_set:
                flag = int(parts[1])
                tid = parts[2]
                pos = int(parts[3])
                mapq = int(parts[4])
                cigar = parts[5]
                
                # Extract tags (format: TAG:TYPE:VALUE)
                tags = {}
                for tag_field in parts[11:]:
                    if ':' in tag_field:
                        tag_parts = tag_field.split(':')
                        if len(tag_parts) >= 3:
                            tag = tag_parts[0]
                            tag_type = tag_parts[1]
                            tag_value = tag_parts[2]
                            if tag_type == 'i':  # Integer
                                tags[tag] = int(tag_value)
                            elif tag_type == 'f':  # Float
                                tags[tag] = float(tag_value)
                            else:
                                tags[tag] = tag_value
                
                results.append({
                    'qname': qname,
                    'flag': flag,
                    'tid': tid,
                    'pos': pos,
                    'mapq': mapq,
                    'cigar': cigar,
                    'NH': tags.get('NH', 0),
                    'HI': tags.get('HI', 0),
                    'AS': tags.get('AS'),
                })
                
                if len(results) >= len(read_names) * 2:  # Account for paired reads
                    break
        
        proc.terminate()
    except Exception as e:
        print(f"Error analyzing BAM: {e}", file=sys.stderr)
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description='Analyze characteristics of differing reads'
    )
    parser.add_argument('--analysis-file', required=True,
                       help='Output from analyze_ec_reads.py')
    parser.add_argument('--bam-file', required=True,
                       help='BAM file to analyze')
    parser.add_argument('--ec-label', required=True,
                       help='EC label to analyze (e.g., "93" or "38,39")')
    parser.add_argument('--max-reads', type=int, default=20,
                       help='Maximum reads to analyze per EC')
    parser.add_argument('--output', required=True,
                       help='Output file')
    
    args = parser.parse_args()
    
    # Extract reads for the specified EC
    ec_to_reads = extract_reads_from_analysis(args.analysis_file, args.max_reads)
    
    if args.ec_label not in ec_to_reads:
        print(f"Error: EC label '{args.ec_label}' not found in analysis file", file=sys.stderr)
        print(f"Available ECs: {list(ec_to_reads.keys())[:10]}...", file=sys.stderr)
        sys.exit(1)
    
    read_names = ec_to_reads[args.ec_label]
    print(f"Analyzing {len(read_names)} reads from EC ({args.ec_label})...")
    
    # Analyze BAM characteristics
    bam_data = analyze_bam_reads(args.bam_file, read_names)
    
    # Group by read name (paired reads)
    reads_by_name = defaultdict(list)
    for data in bam_data:
        reads_by_name[data['qname']].append(data)
    
    # Write analysis
    with open(args.output, 'w') as f:
        f.write(f"Analysis of reads from EC ({args.ec_label})\n")
        f.write(f"Total reads analyzed: {len(reads_by_name)}\n")
        f.write("=" * 80 + "\n\n")
        
        if not bam_data:
            f.write("No reads found in BAM file!\n")
            return
        
        # Summary statistics
        flags = [d['flag'] for d in bam_data]
        mapqs = [d['mapq'] for d in bam_data]
        nhs = [d['NH'] for d in bam_data]
        has_as = sum(1 for d in bam_data if d['AS'] is not None)
        
        f.write("Summary Statistics:\n")
        f.write(f"  Unique read names: {len(reads_by_name)}\n")
        f.write(f"  Total alignments: {len(bam_data)}\n")
        if mapqs:
            f.write(f"  Average MAPQ: {sum(mapqs)/len(mapqs):.1f}\n")
        f.write(f"  Reads with AS tag: {has_as}/{len(bam_data)}\n")
        if nhs:
            f.write(f"  Average NH (number of hits): {sum(nhs)/len(nhs):.1f}\n")
        f.write("\n")
        
        # Flag analysis
        flag_counts = defaultdict(int)
        for flag in flags:
            flag_counts[flag] += 1
        
        f.write("BAM Flag Distribution:\n")
        for flag, count in sorted(flag_counts.items()):
            paired = "paired" if flag & 0x1 else "single"
            proper = "proper" if flag & 0x2 else "not_proper"
            mapped = "mapped" if not (flag & 0x4) else "unmapped"
            f.write(f"  Flag {flag} ({paired}, {proper}, {mapped}): {count}\n")
        f.write("\n")
        
        # Detailed read information
        f.write("Detailed Read Information:\n")
        f.write("=" * 80 + "\n")
        for qname in sorted(reads_by_name.keys())[:args.max_reads]:
            reads = reads_by_name[qname]
            f.write(f"\n{qname}:\n")
            for read in reads:
                f.write(f"  Flag: {read['flag']}, TID: {read['tid']}, Pos: {read['pos']}, ")
                f.write(f"MAPQ: {read['mapq']}, NH: {read['NH']}, HI: {read['HI']}")
                if read['AS'] is not None:
                    f.write(f", AS: {read['AS']}")
                f.write(f"\n    CIGAR: {read['cigar']}\n")
    
    print(f"Analysis written to: {args.output}")

if __name__ == '__main__':
    main()
