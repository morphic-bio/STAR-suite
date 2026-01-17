#!/usr/bin/env python3
"""
Compare T→C mismatches per genomic position for high-divergence genes.

Uses mismatchdetails.tsv files to compare STAR vs GEDI at each position.
"""

import sys
import gzip
from collections import defaultdict

def load_gene_list():
    """Load high-divergence genes from overlap gap report."""
    # Top divergent genes from analysis
    return [
        'ENSG00000172053',  # QARS1
        'ENSG00000233954',  # UQCRHL
        'ENSG00000197061',  # H4C3
        'ENSG00000140525',  # FANCI
        'ENSG00000147403',  # RPL10
        'ENSG00000265681',  # RPL17
        'ENSG00000196365',  # LONP1
    ]

def parse_mismatch_details_by_gene(path: str, target_genes: set):
    """
    Parse mismatch details and group by gene.
    
    Returns: dict[gene_id] -> dict[position] -> {'coverage': float, 'mismatches': float}
    """
    gene_data = defaultdict(lambda: defaultdict(lambda: {'coverage': 0, 'mismatches': 0}))
    
    # First, we need to map genomic positions to genes
    # Since mismatchdetails doesn't have gene info, we'll need to infer from other sources
    # For now, we'll parse all T→C mismatches and note that we need gene mapping
    
    print(f"Note: mismatchdetails.tsv doesn't contain gene IDs.")
    print(f"Will analyze all T→C mismatches and report position distribution.")
    
    position_data = defaultdict(lambda: {'coverage': 0, 'mismatches': 0})
    
    with open(path) as f:
        header = f.readline().strip().split('\t')
        
        try:
            genomic_idx = header.index('Genomic')
            read_idx = header.index('Read')
            pos_idx = header.index('Position')
            cov_idx = header.index('Coverage')
            mismatch_idx = header.index('Mismatches')
        except ValueError as e:
            print(f"Error: Missing column: {e}")
            sys.exit(1)
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) <= max(genomic_idx, read_idx, pos_idx, cov_idx, mismatch_idx):
                continue
            
            genomic = parts[genomic_idx]
            read = parts[read_idx]
            position = int(parts[pos_idx])
            coverage = float(parts[cov_idx])
            mismatches = float(parts[mismatch_idx])
            
            if genomic == 'T' and read == 'C':
                position_data[position]['coverage'] += coverage
                position_data[position]['mismatches'] += mismatches
    
    return position_data

def load_gene_positions_from_gtf(gtf_path: str, target_genes: set):
    """Load gene positions from GTF."""
    gene_positions = {}
    
    opener = gzip.open if gtf_path.endswith('.gz') else open
    with opener(gtf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature = parts[2]
            if feature != 'gene':
                continue
            
            attrs = parts[8]
            gene_id = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1].split('.')[0]
                    break
            
            if gene_id in target_genes:
                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                gene_positions[gene_id] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                }
    
    return gene_positions

def compare_position_distributions(star_data: dict, gedi_data: dict):
    """Compare position distributions between STAR and GEDI."""
    print("\n" + "=" * 80)
    print("Position Distribution Comparison")
    print("=" * 80)
    
    all_positions = sorted(set(star_data.keys()) | set(gedi_data.keys()))
    
    print(f"\n{'Position':<10} {'STAR Cov':>12} {'STAR Mism':>12} {'STAR Rate':>12} {'GEDI Cov':>12} {'GEDI Mism':>12} {'GEDI Rate':>12} {'Rate Ratio':>12}")
    print("-" * 100)
    
    total_star_cov = 0
    total_star_mism = 0
    total_gedi_cov = 0
    total_gedi_mism = 0
    
    for pos in all_positions[:30]:  # Show first 30 positions
        star_cov = star_data.get(pos, {}).get('coverage', 0)
        star_mism = star_data.get(pos, {}).get('mismatches', 0)
        star_rate = star_mism / star_cov if star_cov > 0 else 0
        
        gedi_cov = gedi_data.get(pos, {}).get('coverage', 0)
        gedi_mism = gedi_data.get(pos, {}).get('mismatches', 0)
        gedi_rate = gedi_mism / gedi_cov if gedi_cov > 0 else 0
        
        ratio = star_rate / gedi_rate if gedi_rate > 0 else float('inf') if star_rate > 0 else 0
        
        total_star_cov += star_cov
        total_star_mism += star_mism
        total_gedi_cov += gedi_cov
        total_gedi_mism += gedi_mism
        
        print(f"{pos:<10} {star_cov:>12.0f} {star_mism:>12.1f} {star_rate:>12.6f} {gedi_cov:>12.0f} {gedi_mism:>12.1f} {gedi_rate:>12.6f} {ratio:>12.3f}")
    
    print("-" * 100)
    print(f"{'TOTAL':<10} {total_star_cov:>12.0f} {total_star_mism:>12.1f} {total_star_mism/total_star_cov if total_star_cov > 0 else 0:>12.6f} {total_gedi_cov:>12.0f} {total_gedi_mism:>12.1f} {total_gedi_mism/total_gedi_cov if total_gedi_cov > 0 else 0:>12.6f}")
    
    # Analyze differences
    print("\n" + "=" * 80)
    print("Difference Analysis")
    print("=" * 80)
    
    # Positions where STAR has significantly more
    star_higher = []
    gedi_higher = []
    
    for pos in all_positions:
        star_rate = star_data.get(pos, {}).get('mismatches', 0) / star_data.get(pos, {}).get('coverage', 1) if star_data.get(pos, {}).get('coverage', 0) > 0 else 0
        gedi_rate = gedi_data.get(pos, {}).get('mismatches', 0) / gedi_data.get(pos, {}).get('coverage', 1) if gedi_data.get(pos, {}).get('coverage', 0) > 0 else 0
        
        if star_rate > gedi_rate * 1.5 and star_data.get(pos, {}).get('coverage', 0) > 100:
            star_higher.append((pos, star_rate, gedi_rate))
        elif gedi_rate > star_rate * 1.5 and gedi_data.get(pos, {}).get('coverage', 0) > 100:
            gedi_higher.append((pos, star_rate, gedi_rate))
    
    print(f"\nPositions where STAR rate > 1.5x GEDI rate (with cov > 100): {len(star_higher)}")
    for pos, s_rate, g_rate in star_higher[:10]:
        print(f"  Position {pos}: STAR={s_rate:.6f}, GEDI={g_rate:.6f}, ratio={s_rate/g_rate:.2f}")
    
    print(f"\nPositions where GEDI rate > 1.5x STAR rate (with cov > 100): {len(gedi_higher)}")
    for pos, s_rate, g_rate in gedi_higher[:10]:
        print(f"  Position {pos}: STAR={s_rate:.6f}, GEDI={g_rate:.6f}, ratio={g_rate/s_rate:.2f}")

def main():
    star_path = 'test/tmp_prod_compat/default_SlamQuant.out.mismatchdetails.tsv'
    gedi_path = '/storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3.mismatchdetails.tsv'
    
    print("=" * 80)
    print("Per-Position Mismatch Comparison")
    print("=" * 80)
    
    print("\nNote: mismatchdetails.tsv files contain position-level data but not gene-level.")
    print("Comparing overall position distributions between STAR and GEDI.")
    
    print("\nParsing STAR mismatch details...")
    star_data = parse_mismatch_details_by_gene(star_path, set())
    
    print("Parsing GEDI mismatch details...")
    gedi_data = parse_mismatch_details_by_gene(gedi_path, set())
    
    compare_position_distributions(star_data, gedi_data)
    
    print("\n" + "=" * 80)
    print("CONCLUSION")
    print("=" * 80)
    print("""
This analysis compares T→C mismatch rates by read position across all genes.
For gene-specific analysis, we would need:
1. Genomic position mapping (from alignment files)
2. Gene assignment per read (from debug output)
3. Per-read mismatch calls (from alignment SAM/BAM)

The position-level comparison shows whether differences are:
- Concentrated at specific read positions (e.g., first/last bases)
- Distributed across all positions
- Systematic (consistent ratio) vs random
""")

if __name__ == '__main__':
    main()
