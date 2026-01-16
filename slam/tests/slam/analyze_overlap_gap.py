#!/usr/bin/env python3
"""
Analyze whether STAR vs GEDI divergence is driven by multi-gene assignment on overlapping loci.

Test A: Exclude overlapping genes and recompute correlations
Test B: Aggregate overlapping genes into loci and recompute correlations

Usage:
    python analyze_overlap_gap.py \
        --star-output test/tmp_prod_compat/default_SlamQuant.out \
        --gedi-output /storage/SLAM-Seq-prod-compare-20260109/gedi/WDHD1_0h3_Sense_rerun.tsv.gz \
        --gtf test/fixtures/slam/ref/genes.gtf \
        --output STAR_SLAM_OverlapGap_Report.md
"""

import argparse
import gzip
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple
import sys

# Try to import scipy for Spearman, fall back to manual implementation
try:
    from scipy.stats import pearsonr, spearmanr
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    import math

@dataclass
class GeneData:
    """Normalized gene data for comparison."""
    gene_id: str
    symbol: str
    read_count: float
    ntr: float
    conversions: float
    coverage: float

@dataclass
class GeneAnnotation:
    """Gene annotation from GTF."""
    gene_id: str
    symbol: str
    chrom: str
    strand: str
    start: int
    end: int
    exon_intervals: List[Tuple[int, int]]

def manual_pearson(x: List[float], y: List[float]) -> float:
    """Manual Pearson correlation implementation."""
    n = len(x)
    if n == 0:
        return 0.0
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    
    num = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    denom_x = math.sqrt(sum((xi - mean_x) ** 2 for xi in x))
    denom_y = math.sqrt(sum((yi - mean_y) ** 2 for yi in y))
    
    if denom_x == 0 or denom_y == 0:
        return 0.0
    return num / (denom_x * denom_y)

def manual_spearman(x: List[float], y: List[float]) -> float:
    """Manual Spearman correlation implementation."""
    def rank(data):
        sorted_indices = sorted(range(len(data)), key=lambda i: data[i])
        ranks = [0.0] * len(data)
        for rank_val, idx in enumerate(sorted_indices):
            ranks[idx] = rank_val + 1
        return ranks
    
    rank_x = rank(x)
    rank_y = rank(y)
    return manual_pearson(rank_x, rank_y)

def calc_pearson(x: List[float], y: List[float]) -> float:
    if HAS_SCIPY:
        r, _ = pearsonr(x, y)
        return r
    return manual_pearson(x, y)

def calc_spearman(x: List[float], y: List[float]) -> float:
    if HAS_SCIPY:
        r, _ = spearmanr(x, y)
        return r
    return manual_spearman(x, y)

def format_float(val, fmt=".6f") -> str:
    """Format a float, handling None and NaN gracefully."""
    import math
    if val is None:
        return "N/A"
    try:
        if math.isnan(val):
            return "NaN"
        return f"{val:{fmt}}"
    except (TypeError, ValueError):
        return "N/A"

def parse_gtf(gtf_path: str) -> Dict[str, GeneAnnotation]:
    """Parse GTF file to extract gene and exon information."""
    genes = {}
    exons_by_gene = defaultdict(list)
    
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            start, end = int(start), int(end)
            
            # Parse attributes
            attr_dict = {}
            for attr in attrs.split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, val = attr.split(' ', 1)
                    attr_dict[key] = val.strip('"')
            
            gene_id = attr_dict.get('gene_id', '').split('.')[0]
            gene_name = attr_dict.get('gene_name', gene_id)
            
            if not gene_id:
                continue
            
            if feature == 'gene':
                genes[gene_id] = GeneAnnotation(
                    gene_id=gene_id,
                    symbol=gene_name,
                    chrom=chrom,
                    strand=strand,
                    start=start,
                    end=end,
                    exon_intervals=[]
                )
            elif feature == 'exon':
                exons_by_gene[gene_id].append((start, end))
    
    # Add exons to genes
    for gene_id, exons in exons_by_gene.items():
        if gene_id in genes:
            genes[gene_id].exon_intervals = sorted(exons)
    
    return genes

def find_overlapping_genes(genes: Dict[str, GeneAnnotation], strand_aware: bool = True) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """
    Find genes with overlapping exons.
    
    Returns:
        - Set of all genes that overlap with at least one other gene
        - Dict mapping each gene to its overlapping partners (for clustering)
    """
    overlapping = set()
    overlap_graph = defaultdict(set)
    
    # Group genes by chromosome (and strand if strand-aware)
    groups = defaultdict(list)
    for gene_id, gene in genes.items():
        key = (gene.chrom, gene.strand) if strand_aware else (gene.chrom,)
        groups[key].append(gene)
    
    # Check for overlaps within each group
    for group_genes in groups.values():
        n = len(group_genes)
        for i in range(n):
            for j in range(i + 1, n):
                g1, g2 = group_genes[i], group_genes[j]
                
                # Check if gene ranges overlap first (quick check)
                if g1.end < g2.start or g2.end < g1.start:
                    continue
                
                # Check if any exons overlap
                for e1_start, e1_end in g1.exon_intervals:
                    for e2_start, e2_end in g2.exon_intervals:
                        if not (e1_end < e2_start or e2_end < e1_start):
                            overlapping.add(g1.gene_id)
                            overlapping.add(g2.gene_id)
                            overlap_graph[g1.gene_id].add(g2.gene_id)
                            overlap_graph[g2.gene_id].add(g1.gene_id)
                            break
                    else:
                        continue
                    break
    
    return overlapping, overlap_graph

def build_overlap_clusters(overlap_graph: Dict[str, Set[str]]) -> List[Set[str]]:
    """Build connected components (clusters) from overlap graph."""
    visited = set()
    clusters = []
    
    def dfs(gene_id: str, cluster: Set[str]):
        if gene_id in visited:
            return
        visited.add(gene_id)
        cluster.add(gene_id)
        for neighbor in overlap_graph.get(gene_id, []):
            dfs(neighbor, cluster)
    
    for gene_id in overlap_graph:
        if gene_id not in visited:
            cluster = set()
            dfs(gene_id, cluster)
            if cluster:
                clusters.append(cluster)
    
    return clusters

def parse_star_output(path: str) -> Dict[str, GeneData]:
    """Parse STAR SlamQuant.out file."""
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        gene_idx = header.index('Gene')
        symbol_idx = header.index('Symbol')
        rc_idx = header.index('ReadCount')
        ntr_idx = header.index('NTR')
        conv_idx = header.index('Conversions')
        cov_idx = header.index('Coverage')
        
        for line in f:
            parts = line.strip().split('\t')
            gene_id = parts[gene_idx].split('.')[0]
            data[gene_id] = GeneData(
                gene_id=gene_id,
                symbol=parts[symbol_idx],
                read_count=float(parts[rc_idx]),
                ntr=float(parts[ntr_idx]),
                conversions=float(parts[conv_idx]),
                coverage=float(parts[cov_idx])
            )
    return data

def parse_gedi_output(path: str, use_mean_for_ntr: bool = False) -> Dict[str, GeneData]:
    """Parse GEDI output file.
    
    Args:
        path: Path to GEDI output file
        use_mean_for_ntr: If True, use 'Mean' column for NTR comparison; if False, use 'MAP'
    """
    data = {}
    opener = gzip.open if path.endswith('.gz') else open
    
    with opener(path, 'rt') as f:
        header = f.readline().strip().split('\t')
        gene_idx = 0  # Gene
        symbol_idx = 1  # Symbol
        
        # Find column indices (GEDI has dynamic column names with path prefixes)
        # Be specific to avoid matching wrong columns
        # Take only the FIRST match for each column type
        rc_idx = mean_idx = map_idx = conv_idx = cov_idx = None
        for i, h in enumerate(header):
            h_suffix = h.split()[-1] if h.split() else h  # Get last word (the actual column name)
            if h_suffix == 'Readcount' and rc_idx is None:
                rc_idx = i
            elif h_suffix == 'Mean' and mean_idx is None:
                mean_idx = i
            elif h_suffix == 'MAP' and map_idx is None:
                map_idx = i
            elif h_suffix == 'Conversions' and conv_idx is None:
                conv_idx = i
            elif h_suffix == 'Coverage' and cov_idx is None:  # First Coverage only, not "Double-Hit Coverage"
                cov_idx = i
        
        # Choose NTR column based on flag
        ntr_idx = mean_idx if use_mean_for_ntr else map_idx
        
        for line in f:
            parts = line.strip().split('\t')
            gene_id = parts[gene_idx].split('.')[0]
            data[gene_id] = GeneData(
                gene_id=gene_id,
                symbol=parts[symbol_idx],
                read_count=float(parts[rc_idx]) if rc_idx is not None else 0,
                ntr=float(parts[ntr_idx]) if ntr_idx is not None else 0,
                conversions=float(parts[conv_idx]) if conv_idx is not None else 0,
                coverage=float(parts[cov_idx]) if cov_idx is not None else 0
            )
    return data

def compute_correlations(star_data: Dict[str, GeneData], 
                         gedi_data: Dict[str, GeneData],
                         include_genes: Set[str] = None,
                         thresholds: List[int] = [20, 50, 100]) -> Dict:
    """Compute correlation metrics at various thresholds."""
    results = {}
    
    # Get common genes
    common_genes = set(star_data.keys()) & set(gedi_data.keys())
    if include_genes is not None:
        common_genes &= include_genes
    
    for threshold in thresholds:
        # Filter by GEDI readcount threshold
        filtered_genes = [g for g in common_genes 
                          if gedi_data[g].read_count >= threshold]
        
        if len(filtered_genes) < 3:
            results[threshold] = {
                'n_genes': len(filtered_genes),
                'ntr_pearson': None,
                'ntr_spearman': None,
                'k_nt_pearson': None,
                'k_nt_spearman': None
            }
            continue
        
        # Extract values
        star_ntr = [star_data[g].ntr for g in filtered_genes]
        gedi_ntr = [gedi_data[g].ntr for g in filtered_genes]
        
        # k/nT = conversions / coverage (handle division by zero)
        def safe_k_nt(d: GeneData) -> float:
            return d.conversions / d.coverage if d.coverage > 0 else 0
        
        star_k_nt = [safe_k_nt(star_data[g]) for g in filtered_genes]
        gedi_k_nt = [safe_k_nt(gedi_data[g]) for g in filtered_genes]
        
        results[threshold] = {
            'n_genes': len(filtered_genes),
            'ntr_pearson': calc_pearson(star_ntr, gedi_ntr),
            'ntr_spearman': calc_spearman(star_ntr, gedi_ntr),
            'k_nt_pearson': calc_pearson(star_k_nt, gedi_k_nt),
            'k_nt_spearman': calc_spearman(star_k_nt, gedi_k_nt)
        }
    
    return results

def aggregate_to_loci(data: Dict[str, GeneData], 
                      clusters: List[Set[str]]) -> Dict[str, GeneData]:
    """Aggregate gene data to locus level for overlapping clusters."""
    result = {}
    clustered_genes = set()
    
    for i, cluster in enumerate(clusters):
        locus_id = f"LOCUS_{i:04d}"
        members = [data[g] for g in cluster if g in data]
        
        if not members:
            continue
        
        clustered_genes.update(cluster)
        
        # Aggregate: sum counts, weighted average NTR
        total_rc = sum(m.read_count for m in members)
        total_conv = sum(m.conversions for m in members)
        total_cov = sum(m.coverage for m in members)
        
        # Weighted NTR by read count
        if total_rc > 0:
            weighted_ntr = sum(m.ntr * m.read_count for m in members) / total_rc
        else:
            weighted_ntr = 0
        
        result[locus_id] = GeneData(
            gene_id=locus_id,
            symbol=f"cluster({len(members)})",
            read_count=total_rc,
            ntr=weighted_ntr,
            conversions=total_conv,
            coverage=total_cov
        )
    
    # Add non-clustered genes as-is
    for gene_id, gene_data in data.items():
        if gene_id not in clustered_genes:
            result[gene_id] = gene_data
    
    return result

def format_correlation_table(results: Dict, label: str) -> str:
    """Format correlation results as markdown table."""
    lines = [
        f"**{label}**\n",
        "| Threshold | N Genes | NTR Pearson | NTR Spearman | k/nT Pearson | k/nT Spearman |",
        "|-----------|---------|-------------|--------------|--------------|---------------|"
    ]
    
    for threshold in sorted(results.keys()):
        r = results[threshold]
        if r['ntr_pearson'] is None:
            lines.append(f"| >={threshold} | {r['n_genes']} | N/A | N/A | N/A | N/A |")
        else:
            lines.append(f"| >={threshold} | {r['n_genes']} | {r['ntr_pearson']:.6f} | {r['ntr_spearman']:.6f} | {r['k_nt_pearson']:.6f} | {r['k_nt_spearman']:.6f} |")
    
    return '\n'.join(lines)

def format_delta_table(baseline: Dict, test: Dict, label: str) -> str:
    """Format delta (improvement) table comparing test to baseline."""
    lines = [
        f"**{label} (Delta vs Baseline)**\n",
        "| Threshold | N Genes | ΔNTR Pearson | ΔNTR Spearman | Δk/nT Pearson | Δk/nT Spearman |",
        "|-----------|---------|--------------|---------------|---------------|----------------|"
    ]
    
    for threshold in sorted(baseline.keys()):
        b = baseline[threshold]
        t = test[threshold]
        
        if b['ntr_pearson'] is None or t['ntr_pearson'] is None:
            lines.append(f"| >={threshold} | {t['n_genes']} | N/A | N/A | N/A | N/A |")
        else:
            d_ntr_p = t['ntr_pearson'] - b['ntr_pearson']
            d_ntr_s = t['ntr_spearman'] - b['ntr_spearman']
            d_knt_p = t['k_nt_pearson'] - b['k_nt_pearson']
            d_knt_s = t['k_nt_spearman'] - b['k_nt_spearman']
            lines.append(f"| >={threshold} | {t['n_genes']} | {d_ntr_p:+.6f} | {d_ntr_s:+.6f} | {d_knt_p:+.6f} | {d_knt_s:+.6f} |")
    
    return '\n'.join(lines)

def main():
    parser = argparse.ArgumentParser(description='Analyze STAR vs GEDI overlap gap')
    parser.add_argument('--star-output', required=True, help='STAR SlamQuant.out file')
    parser.add_argument('--gedi-output', required=True, help='GEDI output file (can be .gz)')
    parser.add_argument('--gtf', required=True, help='GTF annotation file')
    parser.add_argument('--output', default='STAR_SLAM_OverlapGap_Report.md', help='Output report file')
    args = parser.parse_args()
    
    print("Loading data...")
    
    # Parse inputs
    print(f"  Parsing GTF: {args.gtf}")
    genes = parse_gtf(args.gtf)
    print(f"    Found {len(genes)} genes")
    
    print(f"  Parsing STAR output: {args.star_output}")
    star_data = parse_star_output(args.star_output)
    print(f"    Found {len(star_data)} genes")
    
    print(f"  Parsing GEDI output: {args.gedi_output}")
    gedi_data_map = parse_gedi_output(args.gedi_output, use_mean_for_ntr=False)
    gedi_data_mean = parse_gedi_output(args.gedi_output, use_mean_for_ntr=True)
    print(f"    Found {len(gedi_data_map)} genes")
    
    # Check coverage parsing
    sample_gene = next(iter(gedi_data_map.values()))
    print(f"    Sample coverage value: {sample_gene.coverage} (should be non-zero for expressed genes)")
    
    # Default to MAP for backward compatibility
    gedi_data = gedi_data_map
    
    # Find overlapping genes
    print("\nFinding overlapping genes...")
    overlapping_genes, overlap_graph = find_overlapping_genes(genes, strand_aware=True)
    print(f"  Found {len(overlapping_genes)} genes with exon overlaps")
    
    # Build overlap clusters
    clusters = build_overlap_clusters(overlap_graph)
    print(f"  Found {len(clusters)} overlap clusters")
    
    # Identify MT genes
    mt_genes = {g for g in genes if genes[g].symbol.startswith('MT-') or genes[g].chrom in ('chrM', 'MT')}
    print(f"  Found {len(mt_genes)} MT genes")
    
    # Identify MT genes that overlap
    mt_overlapping = mt_genes & overlapping_genes
    print(f"  MT genes with overlaps: {len(mt_overlapping)}")
    
    # All genes in both datasets
    all_genes = set(star_data.keys()) & set(gedi_data.keys())
    non_overlap_genes = all_genes - overlapping_genes
    non_mt_genes = all_genes - mt_genes
    non_mt_non_overlap = all_genes - mt_genes - overlapping_genes
    
    print(f"\nGene counts:")
    print(f"  Common genes: {len(all_genes)}")
    print(f"  Non-overlapping genes: {len(non_overlap_genes)}")
    print(f"  Non-MT genes: {len(non_mt_genes)}")
    print(f"  Non-MT, non-overlapping: {len(non_mt_non_overlap)}")
    
    # Compute correlations for each condition
    print("\nComputing correlations...")
    
    results = {}
    
    # First, compare MAP vs Mean to determine best NTR comparator
    print("  Comparing GEDI MAP vs Mean as NTR comparator...")
    results['gedi_map'] = compute_correlations(star_data, gedi_data_map, all_genes)
    results['gedi_mean'] = compute_correlations(star_data, gedi_data_mean, all_genes)
    
    map_ntr_corr = results['gedi_map'][20]['ntr_pearson']
    mean_ntr_corr = results['gedi_mean'][20]['ntr_pearson']
    print(f"    STAR NTR vs GEDI MAP:  {map_ntr_corr:.6f}")
    print(f"    STAR NTR vs GEDI Mean: {mean_ntr_corr:.6f}")
    
    # Use the better comparator
    if mean_ntr_corr and map_ntr_corr and mean_ntr_corr > map_ntr_corr:
        print("    -> Using GEDI Mean as NTR comparator (higher correlation)")
        gedi_data = gedi_data_mean
        ntr_source = "Mean"
    else:
        print("    -> Using GEDI MAP as NTR comparator")
        gedi_data = gedi_data_map
        ntr_source = "MAP"
    
    # Baseline: all genes
    print("  Baseline (all genes)...")
    results['baseline'] = compute_correlations(star_data, gedi_data, all_genes)
    
    # Test A1: Exclude all overlapping genes
    print("  Test A1: Exclude all overlapping genes...")
    results['exclude_overlap'] = compute_correlations(star_data, gedi_data, non_overlap_genes)
    
    # Test A2: Exclude only MT genes
    print("  Test A2: Exclude MT genes only...")
    results['exclude_mt'] = compute_correlations(star_data, gedi_data, non_mt_genes)
    
    # Test A3: Exclude both MT and overlapping genes
    print("  Test A3: Exclude MT and overlapping genes...")
    results['exclude_mt_overlap'] = compute_correlations(star_data, gedi_data, non_mt_non_overlap)
    
    # Test B: Locus-collapsed aggregation
    print("  Test B: Locus-collapsed aggregation...")
    star_loci = aggregate_to_loci(star_data, clusters)
    gedi_loci = aggregate_to_loci(gedi_data, clusters)
    results['locus_collapsed'] = compute_correlations(star_loci, gedi_loci)
    
    # Generate report
    print(f"\nGenerating report: {args.output}")
    
    report = f"""# STAR vs GEDI Overlap Gap Analysis

**Date:** Generated automatically  
**Purpose:** Test whether gene overlaps (especially MT genes) explain the ~0.92 correlation ceiling

---

## Inputs

- **STAR output:** `{args.star_output}`
- **GEDI output:** `{args.gedi_output}`
- **GTF annotation:** `{args.gtf}`

---

## NTR Metric Comparison

GEDI provides both MAP (maximum a posteriori) and Mean (posterior mean) estimates. Comparing to STAR's NTR:

| GEDI Metric | NTR Pearson (>=20) | k/nT Pearson (>=20) |
|-------------|-------------------|---------------------|
| MAP | {results['gedi_map'][20]['ntr_pearson']:.6f} | {format_float(results['gedi_map'][20]['k_nt_pearson'])} |
| Mean | {results['gedi_mean'][20]['ntr_pearson']:.6f} | {format_float(results['gedi_mean'][20]['k_nt_pearson'])} |

**Selected:** GEDI {ntr_source} for NTR comparison (based on higher correlation)

---

## Gene Overlap Statistics

| Category | Count |
|----------|-------|
| Total genes in GTF | {len(genes)} |
| Genes with exon overlaps | {len(overlapping_genes)} |
| Overlap clusters | {len(clusters)} |
| MT genes (chrM) | {len(mt_genes)} |
| MT genes with overlaps | {len(mt_overlapping)} |

### Overlap Clusters (Top 10 by size)

| Cluster | Size | Genes |
|---------|------|-------|
"""
    
    # Add top clusters
    sorted_clusters = sorted(clusters, key=len, reverse=True)[:10]
    for i, cluster in enumerate(sorted_clusters):
        gene_symbols = [genes[g].symbol for g in cluster if g in genes][:5]
        symbols_str = ', '.join(gene_symbols)
        if len(cluster) > 5:
            symbols_str += f", ... (+{len(cluster)-5} more)"
        report += f"| {i+1} | {len(cluster)} | {symbols_str} |\n"
    
    report += f"""
---

## Test A: Overlap-Gene Exclusion

### Baseline (All Genes)

{format_correlation_table(results['baseline'], f"All common genes (n={len(all_genes)})")}

### A1: Exclude All Overlapping Genes

{format_correlation_table(results['exclude_overlap'], f"Non-overlapping genes (n={len(non_overlap_genes)})")}

{format_delta_table(results['baseline'], results['exclude_overlap'], "A1")}

### A2: Exclude MT Genes Only

{format_correlation_table(results['exclude_mt'], f"Non-MT genes (n={len(non_mt_genes)})")}

{format_delta_table(results['baseline'], results['exclude_mt'], "A2")}

### A3: Exclude Both MT and Overlapping Genes

{format_correlation_table(results['exclude_mt_overlap'], f"Non-MT, non-overlapping genes (n={len(non_mt_non_overlap)})")}

{format_delta_table(results['baseline'], results['exclude_mt_overlap'], "A3")}

---

## Test B: Locus-Collapsed Aggregation

Overlapping genes are aggregated into single "locus" units by summing read counts and taking weighted-average NTR.

{format_correlation_table(results['locus_collapsed'], f"Locus-collapsed (clusters merged)")}

{format_delta_table(results['baseline'], results['locus_collapsed'], "B")}

---

## Summary and Interpretation

"""
    
    # Compute summary stats
    baseline_ntr = results['baseline'][20]['ntr_pearson']
    
    for test_name, test_label in [
        ('exclude_overlap', 'Excluding overlapping genes'),
        ('exclude_mt', 'Excluding MT genes'),
        ('exclude_mt_overlap', 'Excluding MT + overlapping'),
        ('locus_collapsed', 'Locus-collapsed aggregation')
    ]:
        test_ntr = results[test_name][20]['ntr_pearson']
        if test_ntr and baseline_ntr:
            delta = test_ntr - baseline_ntr
            report += f"- **{test_label}**: NTR Pearson = {test_ntr:.6f} (Δ = {delta:+.6f})\n"
    
    # Interpretation
    excl_overlap_delta = results['exclude_overlap'][20]['ntr_pearson'] - baseline_ntr if results['exclude_overlap'][20]['ntr_pearson'] else 0
    excl_mt_delta = results['exclude_mt'][20]['ntr_pearson'] - baseline_ntr if results['exclude_mt'][20]['ntr_pearson'] else 0
    locus_delta = results['locus_collapsed'][20]['ntr_pearson'] - baseline_ntr if results['locus_collapsed'][20]['ntr_pearson'] else 0
    
    report += f"""
### Key Findings

"""
    
    if excl_overlap_delta > 0.01:
        report += f"1. **Overlapping genes ARE a significant driver**: Excluding them improves NTR Pearson by {excl_overlap_delta:+.4f}\n"
    else:
        report += f"1. **Overlapping genes are NOT the main driver**: Excluding them changes NTR Pearson by only {excl_overlap_delta:+.4f}\n"
    
    if excl_mt_delta > 0.005:
        report += f"2. **MT genes contribute to divergence**: Excluding them improves NTR Pearson by {excl_mt_delta:+.4f}\n"
    else:
        report += f"2. **MT genes have minimal impact on NTR correlation**: Excluding them changes NTR Pearson by only {excl_mt_delta:+.4f}\n"
    
    if locus_delta > 0.01:
        report += f"3. **Locus aggregation helps**: Collapsing overlaps improves NTR Pearson by {locus_delta:+.4f}\n"
    else:
        report += f"3. **Locus aggregation has limited impact**: NTR Pearson changes by only {locus_delta:+.4f}\n"
    
    report += f"""
### Conclusion

"""
    
    if max(excl_overlap_delta, excl_mt_delta, locus_delta) < 0.02:
        report += """The ~0.92 correlation ceiling is **NOT primarily driven by gene overlap handling**. 
The remaining divergence is likely due to:
- Per-position conversion counting differences
- Read weighting/assignment differences within non-overlapping genes
- Solver/model parameter differences

The multi-gene assignment hypothesis is **not supported** by this analysis.
"""
    else:
        max_improvement = max(excl_overlap_delta, excl_mt_delta, locus_delta)
        report += f"""Gene overlap handling **does contribute** to the divergence (up to {max_improvement:+.4f} improvement).
However, this does not fully close the gap to perfect correlation, indicating other factors also contribute.
"""
    
    # Write report
    with open(args.output, 'w') as f:
        f.write(report)
    
    print(f"\nReport saved to: {args.output}")
    print("\n" + "=" * 60)
    print("QUICK SUMMARY")
    print("=" * 60)
    print(f"Baseline NTR Pearson (>=20): {baseline_ntr:.6f}")
    print(f"Exclude overlapping:         {results['exclude_overlap'][20]['ntr_pearson']:.6f} (Δ = {excl_overlap_delta:+.6f})")
    print(f"Exclude MT:                  {results['exclude_mt'][20]['ntr_pearson']:.6f} (Δ = {excl_mt_delta:+.6f})")
    print(f"Locus-collapsed:             {results['locus_collapsed'][20]['ntr_pearson']:.6f} (Δ = {locus_delta:+.6f})")

if __name__ == '__main__':
    main()
