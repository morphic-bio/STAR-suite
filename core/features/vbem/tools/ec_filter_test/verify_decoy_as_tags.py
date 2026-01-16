#!/usr/bin/env python3
"""
Verify decoy list parity and alignment score extraction (AS tags).

This script checks:
1. Decoy list parity between Salmon and CLI (same transcripts marked as decoys)
2. AS tag presence for all alignments in a small subset of reads

Usage:
    python3 verify_decoy_as_tags.py --bam <input.bam> --salmon-decoys <salmon_decoys.txt> --cli-decoys <cli_decoys.txt> --sample-size <N> --output <report.txt>
"""

import argparse
import sys
import pysam
from collections import defaultdict
from typing import Set, Dict, List

def load_decoy_list(decoy_file: str, name_to_tid: Dict[str, int] = None) -> Set[int]:
    """Load decoy transcript IDs or names from file (one per line).
    
    If name_to_tid is provided, names are mapped to IDs. Otherwise, numeric IDs are used.
    """
    decoys = set()
    with open(decoy_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Try as numeric ID first
                try:
                    tid = int(line)
                    decoys.add(tid)
                except ValueError:
                    # Try as name if name_to_tid mapping provided
                    if name_to_tid and line in name_to_tid:
                        decoys.add(name_to_tid[line])
                    # Otherwise skip (could be invalid)
                    continue
    return decoys

def verify_decoy_parity(salmon_decoys: Set[int], cli_decoys: Set[int]) -> Dict:
    """Compare decoy lists and return differences."""
    salmon_only = salmon_decoys - cli_decoys
    cli_only = cli_decoys - salmon_decoys
    common = salmon_decoys & cli_decoys
    
    return {
        'salmon_count': len(salmon_decoys),
        'cli_count': len(cli_decoys),
        'common_count': len(common),
        'salmon_only': sorted(salmon_only),
        'cli_only': sorted(cli_only),
        'match': len(salmon_only) == 0 and len(cli_only) == 0
    }

def verify_as_tags(bam_file: str, sample_size: int = 1000) -> Dict:
    """Verify AS tags are present for all alignments in a sample of reads."""
    bam = pysam.AlignmentFile(bam_file, 'rb')
    
    reads_without_as = []
    reads_with_as = []
    total_alignments = 0
    alignments_with_as = 0
    alignments_without_as = 0
    
    read_count = 0
    
    for read in bam:
        if read_count >= sample_size:
            break
        
        if read.is_unmapped:
            continue
        
        read_count += 1
        qname = read.query_name
        
        # Check if AS tag is present
        as_tag = read.get_tag('AS') if read.has_tag('AS') else None
        
        if as_tag is None:
            reads_without_as.append({
                'qname': qname,
                'tid': read.reference_id,
                'pos': read.reference_start
            })
            alignments_without_as += 1
        else:
            reads_with_as.append({
                'qname': qname,
                'tid': read.reference_id,
                'as': as_tag
            })
            alignments_with_as += 1
        
        total_alignments += 1
    
    bam.close()
    
    return {
        'total_reads': read_count,
        'total_alignments': total_alignments,
        'alignments_with_as': alignments_with_as,
        'alignments_without_as': alignments_without_as,
        'reads_without_as': reads_without_as[:20],  # Limit to first 20
        'as_tag_present_rate': alignments_with_as / total_alignments if total_alignments > 0 else 0.0
    }

def write_report(decoy_result: Dict, as_result: Dict, output_file: str):
    """Write verification report to file."""
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Decoy List Parity and AS Tag Verification Report\n")
        f.write("=" * 80 + "\n\n")
        
        # Decoy parity section
        f.write("DECOY LIST PARITY:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Salmon decoys: {decoy_result['salmon_count']}\n")
        f.write(f"CLI decoys: {decoy_result['cli_count']}\n")
        f.write(f"Common decoys: {decoy_result['common_count']}\n")
        
        if decoy_result['match']:
            f.write("✅ Decoy lists match perfectly!\n")
        else:
            f.write("⚠️  Decoy lists differ:\n")
            if decoy_result['salmon_only']:
                f.write(f"  Salmon-only decoys ({len(decoy_result['salmon_only'])}): {decoy_result['salmon_only'][:20]}\n")
            if decoy_result['cli_only']:
                f.write(f"  CLI-only decoys ({len(decoy_result['cli_only'])}): {decoy_result['cli_only'][:20]}\n")
        f.write("\n")
        
        # AS tag verification section
        f.write("AS TAG VERIFICATION:\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total reads sampled: {as_result['total_reads']}\n")
        f.write(f"Total alignments: {as_result['total_alignments']}\n")
        f.write(f"Alignments with AS tag: {as_result['alignments_with_as']}\n")
        f.write(f"Alignments without AS tag: {as_result['alignments_without_as']}\n")
        f.write(f"AS tag presence rate: {as_result['as_tag_present_rate']:.2%}\n")
        
        if as_result['alignments_without_as'] == 0:
            f.write("✅ All alignments have AS tags!\n")
        else:
            f.write(f"⚠️  {as_result['alignments_without_as']} alignments missing AS tags\n")
            f.write("First few reads without AS tags:\n")
            for read_info in as_result['reads_without_as']:
                f.write(f"  {read_info['qname']}: tid={read_info['tid']}, pos={read_info['pos']}\n")
        f.write("\n")
        
        # Summary
        f.write("=" * 80 + "\n")
        f.write("SUMMARY:\n")
        f.write("=" * 80 + "\n")
        if decoy_result['match'] and as_result['alignments_without_as'] == 0:
            f.write("✅ All checks passed!\n")
        else:
            f.write("⚠️  Issues found:\n")
            if not decoy_result['match']:
                f.write("  - Decoy lists differ\n")
            if as_result['alignments_without_as'] > 0:
                f.write("  - Some alignments missing AS tags\n")

def load_transcript_names(bam_file: str) -> Dict[str, int]:
    """Load transcript name to ID mapping from BAM header."""
    import pysam
    bam = pysam.AlignmentFile(bam_file, 'rb')
    name_to_tid = {}
    for i, name in enumerate(bam.references):
        name_to_tid[name] = i
    bam.close()
    return name_to_tid

def main():
    parser = argparse.ArgumentParser(description='Verify decoy list parity and AS tag extraction')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--salmon-decoys', required=True, help='Salmon decoy list file (IDs or names)')
    parser.add_argument('--cli-decoys', required=True, help='CLI decoy list file (IDs or names)')
    parser.add_argument('--transcripts', help='Transcriptome FASTA file (optional, for name mapping)')
    parser.add_argument('--sample-size', type=int, default=1000, help='Number of reads to sample for AS tag check')
    parser.add_argument('--output', required=True, help='Output report file')
    
    args = parser.parse_args()
    
    # Load transcript name to ID mapping from BAM header
    print(f"Loading transcript names from BAM header: {args.bam}")
    name_to_tid = load_transcript_names(args.bam)
    print(f"  Found {len(name_to_tid)} transcripts")
    
    print(f"Loading Salmon decoy list: {args.salmon_decoys}")
    salmon_decoys = load_decoy_list(args.salmon_decoys, name_to_tid)
    print(f"  Found {len(salmon_decoys)} decoys")
    
    print(f"Loading CLI decoy list: {args.cli_decoys}")
    cli_decoys = load_decoy_list(args.cli_decoys, name_to_tid)
    print(f"  Found {len(cli_decoys)} decoys")
    
    print("Comparing decoy lists...")
    decoy_result = verify_decoy_parity(salmon_decoys, cli_decoys)
    
    print(f"Verifying AS tags in BAM file: {args.bam}")
    print(f"  Sampling {args.sample_size} reads...")
    as_result = verify_as_tags(args.bam, args.sample_size)
    
    print(f"Writing report to: {args.output}")
    write_report(decoy_result, as_result, args.output)
    
    if decoy_result['match'] and as_result['alignments_without_as'] == 0:
        print("\n✅ All checks passed!")
        return 0
    else:
        print("\n⚠️  Issues found - see report for details")
        return 1

if __name__ == '__main__':
    sys.exit(main())
