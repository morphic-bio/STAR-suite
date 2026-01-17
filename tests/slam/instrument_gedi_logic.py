#!/usr/bin/env python3
"""
Instrument GEDI's mismatch filtering logic to understand why positions are skipped.

This script replicates GEDI's SlamCollector logic to log:
- trim (mpos<trim5p or mpos>=readLen-trim5p)
- overlap (isPositionInOverlap)
- no4sU (no4su mask)
- masked (SNP mask)
- counted

For SE reads, overlap filtering doesn't apply.
Since trim5p=0 and trim3p=0 were used, trim filtering doesn't apply either.

This leaves: masked SNPs and the no4sU condition.
"""

import sys
import gzip
import pysam
from collections import defaultdict

def parse_snp_mask(snp_bed_path):
    """Parse SNP mask BED file."""
    masked = set()
    if not snp_bed_path:
        return masked
    
    opener = gzip.open if snp_bed_path.endswith('.gz') else open
    with opener(snp_bed_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                for pos in range(start, end):
                    masked.add((chrom, pos))
    return masked

def analyze_read_positions(bam_path, reference_path=None, snp_mask=None, max_reads=10000, trim5p=0, trim3p=0):
    """
    Analyze read positions and log why each T position is counted or skipped.
    
    Replicates GEDI's logic from SlamCollector.java
    """
    
    results = defaultdict(lambda: {'counted': 0, 'skipped_trim': 0, 'skipped_mask': 0, 'skipped_quality': 0, 'skipped_other': 0})
    position_details = []  # (read_name, genomic_pos, read_pos, status, reason)
    
    bam = pysam.AlignmentFile(bam_path, 'rb')
    
    read_count = 0
    for read in bam.fetch():
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        
        read_count += 1
        if read_count > max_reads:
            break
        
        # Get aligned pairs (read_pos, ref_pos, ref_base)
        # read_pos is 0-based position in the read
        # ref_pos is 0-based position in the reference
        
        read_len = read.query_length
        chrom = read.reference_name
        
        # GEDI's trim logic (disabled if trim5p=0 and trim3p=0)
        # if mpos<trim5p || mpos>=readLen-trim5p: skip
        # if (readLen-mpos<trim3p && mpos<readLen): skip  (for 3' trim within mate)
        
        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
            if read_pos is None or ref_pos is None:
                continue  # Insertion or deletion
            
            if ref_base is None:
                continue
            
            ref_base_upper = ref_base.upper()
            
            # Only look at T positions (sense strand)
            # For antisense, we'd look at A positions
            is_reverse = read.is_reverse
            
            if is_reverse:
                # On minus strand, T in reference means A in the sense transcript
                # Conversions would be A→G in the read (complement of T→C)
                target_base = 'A'
            else:
                target_base = 'T'
            
            if ref_base_upper != target_base:
                continue
            
            # Now check read base
            read_base = read.query_sequence[read_pos].upper() if read.query_sequence else None
            if read_base is None:
                continue
            
            # Check for T→C conversion (or A→G for reverse strand)
            if is_reverse:
                is_conversion = (ref_base_upper == 'A' and read_base == 'G')
            else:
                is_conversion = (ref_base_upper == 'T' and read_base == 'C')
            
            # GEDI's position mapping: mapToRead()
            # For SE reads without indels, this is just the read position
            # For reads with insertions: positions after insertion are shifted
            # For reads with deletions: positions are not shifted (deletion doesn't add read bases)
            
            # In GEDI, mapToRead1 adjusts for insertions before this position
            # Let's compute the "GEDI-style" position
            gedi_read_pos = read_pos
            
            # Check insertions before this position
            cigar_tuples = read.cigartuples
            if cigar_tuples:
                consumed_read = 0
                for op, length in cigar_tuples:
                    if consumed_read >= read_pos:
                        break
                    if op in [0, 1, 4, 7, 8]:  # M, I, S, =, X consume read
                        if op == 1:  # Insertion
                            # GEDI adds insertion length to position mapping
                            gedi_read_pos += length
                        consumed_read += length
            
            # Apply trim filtering (even though disabled in this run)
            status = 'counted'
            reason = None
            
            if trim5p > 0 or trim3p > 0:
                if gedi_read_pos < trim5p:
                    status = 'skipped'
                    reason = 'trim5p'
                elif gedi_read_pos >= read_len - trim5p:
                    status = 'skipped'
                    reason = 'trim_end'
            
            # Check SNP mask
            if status == 'counted' and snp_mask and (chrom, ref_pos) in snp_mask:
                status = 'skipped'
                reason = 'masked_snp'
            
            # Check base quality (GEDI doesn't seem to have explicit quality filter in SlamCollector)
            # But the aligner may have filtered low-quality bases already
            
            # Record result
            results[read_pos][status] += 1
            if status == 'skipped':
                results[read_pos][f'skipped_{reason}'] = results[read_pos].get(f'skipped_{reason}', 0) + 1
            
            # For detailed logging (first 100 reads only)
            if read_count <= 100:
                position_details.append({
                    'read_name': read.query_name,
                    'chrom': chrom,
                    'genomic_pos': ref_pos,
                    'read_pos': read_pos,
                    'gedi_read_pos': gedi_read_pos,
                    'ref_base': ref_base_upper,
                    'read_base': read_base,
                    'is_conversion': is_conversion,
                    'status': status,
                    'reason': reason,
                    'is_reverse': is_reverse
                })
    
    bam.close()
    return results, position_details

def main():
    bam_path = '/storage/SLAM-Seq-prod-compare-20260109/star/WDHD1_0h3_Aligned.sortedByCoord.out.bam'
    
    print("="*100)
    print("GEDI Logic Instrumentation - Why Mismatches are Skipped")
    print("="*100)
    print(f"\nAnalyzing: {bam_path}")
    print(f"trim5p=0, trim3p=0 (disabled)")
    print(f"Max reads: 10000")
    
    results, details = analyze_read_positions(bam_path, trim5p=0, trim3p=0, max_reads=10000)
    
    print("\n" + "="*80)
    print("Position Summary (T positions only)")
    print("="*80)
    print(f"{'Pos':>4} | {'Counted':>10} | {'Skip_Trim':>10} | {'Skip_Mask':>10} | {'Skip_Other':>10}")
    print("-"*60)
    
    for pos in sorted(results.keys()):
        r = results[pos]
        print(f"{pos:>4} | {r['counted']:>10} | {r.get('skipped_trim5p', 0) + r.get('skipped_trim_end', 0):>10} | {r.get('skipped_masked_snp', 0):>10} | {r.get('skipped_other', 0):>10}")
    
    print("\n" + "="*80)
    print("Sample Position Details (first 100 reads)")
    print("="*80)
    
    # Show T→C conversions
    conversions = [d for d in details if d['is_conversion']]
    print(f"\nT→C Conversions found: {len(conversions)}")
    print(f"{'Read':>20} | {'GenPos':>12} | {'ReadPos':>8} | {'GEDIPos':>8} | {'Status':>10} | {'Reason':>12}")
    print("-"*90)
    
    for d in conversions[:20]:
        print(f"{d['read_name'][:20]:>20} | {d['chrom']}:{d['genomic_pos']:>8} | {d['read_pos']:>8} | {d['gedi_read_pos']:>8} | {d['status']:>10} | {d['reason'] or 'n/a':>12}")
    
    # Check if read_pos != gedi_read_pos anywhere
    print("\n" + "="*80)
    print("Position Mapping Differences (read_pos != gedi_read_pos)")
    print("="*80)
    
    diffs = [d for d in details if d['read_pos'] != d['gedi_read_pos']]
    print(f"Reads with position mapping differences: {len(diffs)}")
    for d in diffs[:10]:
        print(f"  {d['read_name']}: read_pos={d['read_pos']}, gedi_read_pos={d['gedi_read_pos']}, diff={d['gedi_read_pos']-d['read_pos']}")

if __name__ == '__main__':
    main()
