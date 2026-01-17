#!/usr/bin/env python3
"""
Compare read position handling between STAR and GEDI by analyzing:
1. The same BAM file used by both tools
2. T→C conversion positions in specific reads
3. How read positions map to the position histograms

Goal: Determine if mapToRead() is causing the position discrepancy.
"""

import sys
import pysam
from collections import defaultdict

def analyze_single_read(read, verbose=False):
    """
    Analyze a single read's T and T→C positions.
    
    Returns dict with position details for comparison with GEDI.
    """
    if read.is_unmapped:
        return None
    
    result = {
        'read_name': read.query_name,
        'chrom': read.reference_name,
        'start': read.reference_start,
        'end': read.reference_end,
        'strand': '-' if read.is_reverse else '+',
        'cigar': read.cigarstring,
        'read_len': read.query_length,
        't_positions': [],  # List of (read_pos, genomic_pos, is_conversion)
        'tc_count': 0,
        't_count': 0
    }
    
    for read_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
        if read_pos is None or ref_pos is None or ref_base is None:
            continue
        
        ref_base_upper = ref_base.upper()
        read_base = read.query_sequence[read_pos].upper() if read.query_sequence else None
        
        if read_base is None:
            continue
        
        # For sense strand, look for T→C
        # For antisense, the genomic T is actually an A in the transcript
        is_reverse = read.is_reverse
        
        if is_reverse:
            # Minus strand: genomic A corresponds to transcript T
            # Conversion would be A→G in the read (complement of T→C)
            if ref_base_upper == 'A':
                result['t_count'] += 1
                is_conversion = (read_base == 'G')
                if is_conversion:
                    result['tc_count'] += 1
                result['t_positions'].append({
                    'read_pos': read_pos,
                    'genomic_pos': ref_pos,
                    'ref_base': ref_base_upper,
                    'read_base': read_base,
                    'is_conversion': is_conversion
                })
        else:
            # Plus strand: look for T→C
            if ref_base_upper == 'T':
                result['t_count'] += 1
                is_conversion = (read_base == 'C')
                if is_conversion:
                    result['tc_count'] += 1
                result['t_positions'].append({
                    'read_pos': read_pos,
                    'genomic_pos': ref_pos,
                    'ref_base': ref_base_upper,
                    'read_base': read_base,
                    'is_conversion': is_conversion
                })
    
    return result

def aggregate_positions(reads_analysis):
    """
    Aggregate T and T→C counts by read position.
    
    This should match STAR's mismatchdetails if STAR uses direct read positions.
    """
    by_position = defaultdict(lambda: {'t_count': 0, 'tc_count': 0})
    
    for read_data in reads_analysis:
        if read_data is None:
            continue
        for t_pos in read_data['t_positions']:
            pos = t_pos['read_pos']
            by_position[pos]['t_count'] += 1
            if t_pos['is_conversion']:
                by_position[pos]['tc_count'] += 1
    
    return by_position

def main():
    bam_path = '/storage/SLAM-Seq-prod-compare-20260109/star/WDHD1_0h3_Aligned.sortedByCoord.out.bam'
    
    print("="*100)
    print("Direct Read Position Analysis from BAM")
    print("="*100)
    print(f"\nAnalyzing: {bam_path}")
    
    bam = pysam.AlignmentFile(bam_path, 'rb')
    
    # Analyze first N reads
    max_reads = 50000
    reads_analysis = []
    
    for i, read in enumerate(bam.fetch()):
        if i >= max_reads:
            break
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        
        result = analyze_single_read(read)
        if result:
            reads_analysis.append(result)
    
    bam.close()
    
    print(f"\nAnalyzed {len(reads_analysis)} primary alignments")
    
    # Show some example reads with conversions
    reads_with_conv = [r for r in reads_analysis if r['tc_count'] > 0]
    print(f"Reads with T→C conversions: {len(reads_with_conv)}")
    
    print("\n" + "="*80)
    print("Example Reads with Conversions")
    print("="*80)
    
    for r in reads_with_conv[:5]:
        print(f"\n{r['read_name']}")
        print(f"  {r['chrom']}:{r['start']}-{r['end']} ({r['strand']})")
        print(f"  CIGAR: {r['cigar']}")
        print(f"  T positions: {r['t_count']}, T→C conversions: {r['tc_count']}")
        print(f"  Conversion positions:")
        for t in r['t_positions']:
            if t['is_conversion']:
                print(f"    read_pos={t['read_pos']}, genomic_pos={t['genomic_pos']}, {t['ref_base']}→{t['read_base']}")
    
    # Aggregate by position
    by_pos = aggregate_positions(reads_analysis)
    
    print("\n" + "="*80)
    print("Position Aggregation (from BAM analysis)")
    print("="*80)
    print(f"{'Pos':>4} | {'T_Count':>10} | {'TC_Count':>10} | {'Rate':>10}")
    print("-"*50)
    
    for pos in sorted(by_pos.keys()):
        t = by_pos[pos]['t_count']
        tc = by_pos[pos]['tc_count']
        rate = tc / t * 100 if t > 0 else 0
        if pos in [0, 5, 10, 20, 30, 40, 45, 49] or rate > 5:
            print(f"{pos:>4} | {t:>10} | {tc:>10} | {rate:>9.4f}%")
    
    # Summary comparison
    print("\n" + "="*80)
    print("Position Rate Summary")
    print("="*80)
    
    # Early positions (0-10)
    early_t = sum(by_pos[p]['t_count'] for p in range(11))
    early_tc = sum(by_pos[p]['tc_count'] for p in range(11))
    
    # Late positions (40-49)
    late_t = sum(by_pos[p]['t_count'] for p in range(40, 50))
    late_tc = sum(by_pos[p]['tc_count'] for p in range(40, 50))
    
    print(f"\nEarly (0-10): T={early_t}, T→C={early_tc}, Rate={early_tc/early_t*100 if early_t > 0 else 0:.4f}%")
    print(f"Late (40-49): T={late_t}, T→C={late_tc}, Rate={late_tc/late_t*100 if late_t > 0 else 0:.4f}%")
    
    if early_t > 0 and late_t > 0:
        early_rate = early_tc / early_t
        late_rate = late_tc / late_t
        print(f"\nLate/Early ratio: {late_rate/early_rate if early_rate > 0 else 'inf':.2f}")
    
    print("\n" + "="*80)
    print("Comparison with STAR mismatchdetails")
    print("="*80)
    print("""
If this BAM analysis shows similar late-position rates as STAR mismatchdetails,
then the issue is NOT in STAR's position handling - it's in the BAM data itself.

If this analysis shows LOW late-position rates (like GEDI), then STAR is doing
something different in its position aggregation.
""")

if __name__ == '__main__':
    main()
