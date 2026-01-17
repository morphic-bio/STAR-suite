#!/usr/bin/env python3
"""
Analyze how soft-clips affect position mapping.

Key question: Does GEDI include soft-clipped bases in its read position?
If so, position 49 in the aligned region might be position 49+clip_len in GEDI.

This script computes two position aggregations:
1. aligned_pos: position within the aligned portion (ignores soft-clips)
2. full_read_pos: position within the full read (includes soft-clips)
"""

import sys
import pysam
from collections import defaultdict

def get_clip_lengths(cigar):
    """Get 5' and 3' soft-clip lengths from CIGAR."""
    if not cigar:
        return 0, 0
    
    clip_5p = cigar[0][1] if cigar[0][0] == 4 else 0  # S at start
    clip_3p = cigar[-1][1] if cigar[-1][0] == 4 else 0  # S at end
    
    return clip_5p, clip_3p

def analyze_positions(bam_path, max_reads=50000):
    """
    Analyze T→C rates using two position definitions:
    1. aligned_pos: position in aligned portion
    2. full_pos: position in full read (aligned_pos + 5' clip offset)
    """
    
    bam = pysam.AlignmentFile(bam_path, 'rb')
    
    # Two aggregations
    aligned_by_pos = defaultdict(lambda: {'t': 0, 'tc': 0})
    full_by_pos = defaultdict(lambda: {'t': 0, 'tc': 0})
    
    read_count = 0
    
    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        
        read_count += 1
        if read_count > max_reads:
            break
        
        cigar = read.cigartuples
        if not cigar:
            continue
        
        clip_5p, clip_3p = get_clip_lengths(cigar)
        
        seq = read.query_sequence
        if not seq:
            continue
        
        is_reverse = read.is_reverse
        
        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
            if read_pos is None or ref_pos is None or ref_base is None:
                continue
            
            ref_base_upper = ref_base.upper()
            read_base = seq[read_pos].upper()
            
            # Position in aligned portion (what STAR likely uses)
            # This is the position WITHIN the M/I operations, not counting soft-clips
            aligned_pos = read_pos - clip_5p  # Adjust for 5' clip to get aligned-only position
            
            # Position in full read (what GEDI might use with mapToRead)
            full_pos = read_pos  # Direct position in the full read sequence
            
            # Only count positive aligned positions (after soft-clip)
            if aligned_pos < 0:
                continue
            
            # Track T positions
            if is_reverse:
                if ref_base_upper == 'A':
                    aligned_by_pos[aligned_pos]['t'] += 1
                    full_by_pos[full_pos]['t'] += 1
                    if read_base == 'G':
                        aligned_by_pos[aligned_pos]['tc'] += 1
                        full_by_pos[full_pos]['tc'] += 1
            else:
                if ref_base_upper == 'T':
                    aligned_by_pos[aligned_pos]['t'] += 1
                    full_by_pos[full_pos]['t'] += 1
                    if read_base == 'C':
                        aligned_by_pos[aligned_pos]['tc'] += 1
                        full_by_pos[full_pos]['tc'] += 1
    
    bam.close()
    
    return aligned_by_pos, full_by_pos, read_count

def main():
    bam_path = '/storage/SLAM-Seq-prod-compare-20260109/star/WDHD1_0h3_Aligned.sortedByCoord.out.bam'
    
    print("="*100)
    print("Soft-Clip Position Shift Analysis")
    print("="*100)
    
    aligned_by_pos, full_by_pos, read_count = analyze_positions(bam_path, max_reads=50000)
    
    print(f"\nAnalyzed {read_count} reads")
    
    print("\n" + "="*80)
    print("Comparison: Aligned Position vs Full Read Position")
    print("="*80)
    print(f"{'Pos':>4} | {'ALIGNED_T':>10} | {'ALIGNED_TC':>10} | {'ALIGN_Rate':>10} | {'FULL_T':>10} | {'FULL_TC':>10} | {'FULL_Rate':>10}")
    print("-"*95)
    
    all_positions = sorted(set(aligned_by_pos.keys()) | set(full_by_pos.keys()))
    
    for pos in all_positions:
        if pos > 55:
            continue
        if pos not in [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 49, 50]:
            continue
        
        at = aligned_by_pos[pos]['t']
        atc = aligned_by_pos[pos]['tc']
        ar = atc / at * 100 if at > 0 else 0
        
        ft = full_by_pos[pos]['t']
        ftc = full_by_pos[pos]['tc']
        fr = ftc / ft * 100 if ft > 0 else 0
        
        print(f"{pos:>4} | {at:>10} | {atc:>10} | {ar:>9.4f}% | {ft:>10} | {ftc:>10} | {fr:>9.4f}%")
    
    print("\n" + "="*80)
    print("Summary by Position Range")
    print("="*80)
    
    # Early positions (0-10) in aligned coordinates
    a_early_t = sum(aligned_by_pos[p]['t'] for p in range(11))
    a_early_tc = sum(aligned_by_pos[p]['tc'] for p in range(11))
    a_early_rate = a_early_tc / a_early_t * 100 if a_early_t > 0 else 0
    
    # Late positions (40-49) in aligned coordinates
    a_late_t = sum(aligned_by_pos[p]['t'] for p in range(40, 50))
    a_late_tc = sum(aligned_by_pos[p]['tc'] for p in range(40, 50))
    a_late_rate = a_late_tc / a_late_t * 100 if a_late_t > 0 else 0
    
    # Early in full coordinates
    f_early_t = sum(full_by_pos[p]['t'] for p in range(11))
    f_early_tc = sum(full_by_pos[p]['tc'] for p in range(11))
    f_early_rate = f_early_tc / f_early_t * 100 if f_early_t > 0 else 0
    
    # Late in full coordinates
    f_late_t = sum(full_by_pos[p]['t'] for p in range(40, 50))
    f_late_tc = sum(full_by_pos[p]['tc'] for p in range(40, 50))
    f_late_rate = f_late_tc / f_late_t * 100 if f_late_t > 0 else 0
    
    print(f"\nAligned Position Coordinates:")
    print(f"  Early (0-10): T={a_early_t}, T→C={a_early_tc}, Rate={a_early_rate:.4f}%")
    print(f"  Late (40-49): T={a_late_t}, T→C={a_late_tc}, Rate={a_late_rate:.4f}%")
    print(f"  Late/Early ratio: {a_late_rate/a_early_rate if a_early_rate > 0 else 0:.2f}")
    
    print(f"\nFull Read Position Coordinates:")
    print(f"  Early (0-10): T={f_early_t}, T→C={f_early_tc}, Rate={f_early_rate:.4f}%")
    print(f"  Late (40-49): T={f_late_t}, T→C={f_late_tc}, Rate={f_late_rate:.4f}%")
    print(f"  Late/Early ratio: {f_late_rate/f_early_rate if f_early_rate > 0 else 0:.2f}")
    
    print("\n" + "="*80)
    print("INTERPRETATION")
    print("="*80)
    print("""
If "Aligned Position" rates match STAR's mismatchdetails:
  → STAR uses positions within the aligned portion, ignoring soft-clips.

If "Full Read Position" rates match GEDI's mismatchdetails:
  → GEDI uses positions within the full read, including soft-clip offset.

The key insight: With 60% of reads having 3' soft-clips, position 49 in
"aligned coordinates" may represent positions 49+clip_len in "full read
coordinates", which could shift into positions beyond 49.

This could explain why STAR and GEDI disagree at late positions.
""")

if __name__ == '__main__':
    main()
