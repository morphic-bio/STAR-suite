#!/usr/bin/env python3
"""
Analyze soft-clip distribution and its impact on position counting.

STAR may be counting soft-clipped positions differently than GEDI.
"""

import sys
import pysam
from collections import defaultdict

def analyze_softclips(bam_path, max_reads=50000):
    """Analyze soft-clip patterns in the BAM."""
    
    bam = pysam.AlignmentFile(bam_path, 'rb')
    
    stats = {
        'total_reads': 0,
        'reads_with_5p_clip': 0,
        'reads_with_3p_clip': 0,
        'reads_with_both_clips': 0,
        'reads_no_clips': 0,
        '5p_clip_lengths': defaultdict(int),
        '3p_clip_lengths': defaultdict(int),
    }
    
    # Position-based analysis
    # mapped_pos: positions in the aligned portion (not soft-clipped)
    # all_pos: all positions including soft-clipped
    mapped_by_pos = defaultdict(lambda: {'t': 0, 'tc': 0})
    all_by_pos = defaultdict(lambda: {'t': 0, 'tc': 0})
    
    for i, read in enumerate(bam.fetch()):
        if i >= max_reads:
            break
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        
        stats['total_reads'] += 1
        
        cigar = read.cigartuples
        if not cigar:
            continue
        
        # Check for soft-clips
        has_5p_clip = cigar[0][0] == 4  # S at start
        has_3p_clip = cigar[-1][0] == 4  # S at end
        
        clip_5p_len = cigar[0][1] if has_5p_clip else 0
        clip_3p_len = cigar[-1][1] if has_3p_clip else 0
        
        if has_5p_clip and has_3p_clip:
            stats['reads_with_both_clips'] += 1
        elif has_5p_clip:
            stats['reads_with_5p_clip'] += 1
        elif has_3p_clip:
            stats['reads_with_3p_clip'] += 1
        else:
            stats['reads_no_clips'] += 1
        
        if clip_5p_len > 0:
            stats['5p_clip_lengths'][clip_5p_len] += 1
        if clip_3p_len > 0:
            stats['3p_clip_lengths'][clip_3p_len] += 1
        
        # Get sequence and quality
        seq = read.query_sequence
        if not seq:
            continue
        
        read_len = len(seq)
        is_reverse = read.is_reverse
        
        # For mapped positions, use aligned pairs
        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
            if read_pos is None or ref_pos is None or ref_base is None:
                continue
            
            ref_base_upper = ref_base.upper()
            read_base = seq[read_pos].upper()
            
            # Track T positions in mapped region
            if is_reverse:
                if ref_base_upper == 'A':
                    mapped_by_pos[read_pos]['t'] += 1
                    if read_base == 'G':
                        mapped_by_pos[read_pos]['tc'] += 1
            else:
                if ref_base_upper == 'T':
                    mapped_by_pos[read_pos]['t'] += 1
                    if read_base == 'C':
                        mapped_by_pos[read_pos]['tc'] += 1
        
        # For ALL positions (including soft-clips), we can only look at read bases
        # We don't have reference for soft-clipped regions
        # But we can count by position in the read
        for read_pos in range(read_len):
            read_base = seq[read_pos].upper()
            # For all positions, just count T bases (we can't know if they're conversions in soft-clips)
            if read_base == 'T':
                all_by_pos[read_pos]['t'] += 1
            elif read_base == 'C':
                # Can't tell if C is a conversion without reference
                pass
    
    bam.close()
    
    return stats, mapped_by_pos, all_by_pos

def main():
    bam_path = '/storage/SLAM-Seq-prod-compare-20260109/star/WDHD1_0h3_Aligned.sortedByCoord.out.bam'
    
    print("="*100)
    print("Soft-Clip Analysis")
    print("="*100)
    
    stats, mapped_by_pos, all_by_pos = analyze_softclips(bam_path, max_reads=50000)
    
    print(f"\nTotal reads analyzed: {stats['total_reads']}")
    print(f"\nSoft-clip distribution:")
    print(f"  5' clip only:  {stats['reads_with_5p_clip']} ({stats['reads_with_5p_clip']/stats['total_reads']*100:.1f}%)")
    print(f"  3' clip only:  {stats['reads_with_3p_clip']} ({stats['reads_with_3p_clip']/stats['total_reads']*100:.1f}%)")
    print(f"  Both clips:    {stats['reads_with_both_clips']} ({stats['reads_with_both_clips']/stats['total_reads']*100:.1f}%)")
    print(f"  No clips:      {stats['reads_no_clips']} ({stats['reads_no_clips']/stats['total_reads']*100:.1f}%)")
    
    print(f"\n5' clip length distribution:")
    for length in sorted(stats['5p_clip_lengths'].keys())[:10]:
        count = stats['5p_clip_lengths'][length]
        print(f"  {length}bp: {count} reads")
    
    print(f"\n3' clip length distribution:")
    for length in sorted(stats['3p_clip_lengths'].keys())[:10]:
        count = stats['3p_clip_lengths'][length]
        print(f"  {length}bp: {count} reads")
    
    print("\n" + "="*80)
    print("T→C Rates by Position (Mapped Regions Only)")
    print("="*80)
    print(f"{'Pos':>4} | {'T_Count':>10} | {'TC_Count':>10} | {'Rate':>10}")
    print("-"*50)
    
    for pos in sorted(mapped_by_pos.keys()):
        t = mapped_by_pos[pos]['t']
        tc = mapped_by_pos[pos]['tc']
        rate = tc / t * 100 if t > 0 else 0
        if pos in [0, 5, 10, 20, 30, 40, 45, 49]:
            print(f"{pos:>4} | {t:>10} | {tc:>10} | {rate:>9.4f}%")
    
    print("\n" + "="*80)
    print("KEY INSIGHT")
    print("="*80)
    print("""
Soft-clips affect position counting:
- If a read has a 10bp 5' soft-clip, position 0 in the aligned region
  is actually position 10 in the full read sequence.
- STAR and GEDI may differ in how they handle soft-clipped positions.
- This could explain the position-dependent rate differences.

If GEDI uses mapToRead() which accounts for soft-clips, and STAR uses
direct alignment positions, the position histograms would differ.
""")

if __name__ == '__main__':
    main()
