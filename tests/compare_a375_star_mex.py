#!/usr/bin/env python3
"""
Compare A375 Cell Ranger MEX output with STAR MEX output.

Harmonizes cell barcodes, compares gene counts, and calculates correlations.
For Spearman correlation, only includes genes with >= 10 counts in both datasets.
"""

import argparse
import gzip
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from collections import defaultdict


def open_file(path):
    """Open a file, handling gzip compression."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def read_lines(path):
    """Read lines from a file (handles gzip)."""
    with open_file(path) as f:
        return [line.rstrip("\n") for line in f]


def read_mtx(path):
    """Read Matrix Market format file. Returns (nrows, ncols, nnz, entries)."""
    entries = []
    with open_file(path) as handle:
        header = handle.readline().strip()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"Invalid MatrixMarket header: {header}")
        line = handle.readline().strip()
        while line.startswith("%"):
            line = handle.readline().strip()
        if not line:
            raise ValueError("Missing matrix dimensions")
        nrows, ncols, nnz = [int(x) for x in line.split()]
        for raw in handle:
            raw = raw.strip()
            if not raw:
                continue
            parts = raw.split()
            if len(parts) < 3:
                continue
            row = int(parts[0])  # 1-based
            col = int(parts[1])  # 1-based
            val = int(float(parts[2]))
            entries.append((row, col, val))
    return nrows, ncols, nnz, entries


def resolve_mex_file(mex_dir, basename):
    """Resolve MEX file path, checking for .gz extension."""
    plain = os.path.join(mex_dir, basename)
    gz = plain + ".gz"
    if os.path.exists(plain):
        return plain
    if os.path.exists(gz):
        return gz
    raise ValueError(f"Missing {basename}(.gz) in {mex_dir}")


def harmonize_barcode(bc):
    """Harmonize barcode by removing common suffixes like -1, -2, etc."""
    # Remove trailing -N pattern (e.g., -1, -2)
    if '-' in bc:
        parts = bc.rsplit('-', 1)
        if parts[1].isdigit():
            return parts[0]
    return bc


def load_mex(mex_dir):
    """Load MEX directory. Returns (features_df, barcodes_list, matrix_dict)."""
    # Read features.tsv
    features_path = resolve_mex_file(mex_dir, "features.tsv")
    features_lines = read_lines(features_path)
    features_data = []
    for line in features_lines:
        parts = line.split("\t")
        if len(parts) >= 1:
            feature_id = parts[0]
            feature_name = parts[1] if len(parts) > 1 else feature_id
            feature_type = parts[2] if len(parts) > 2 else "Gene Expression"
            features_data.append({
                'id': feature_id,
                'name': feature_name,
                'type': feature_type
            })
    features_df = pd.DataFrame(features_data)
    
    # Read barcodes.tsv
    barcodes_path = resolve_mex_file(mex_dir, "barcodes.tsv")
    barcodes_list = read_lines(barcodes_path)
    
    # Read matrix.mtx
    matrix_path = resolve_mex_file(mex_dir, "matrix.mtx")
    nrows, ncols, nnz, entries = read_mtx(matrix_path)
    
    # Build sparse matrix as dict: (gene_idx, barcode_idx) -> count
    matrix_dict = {}
    for row, col, val in entries:
        gene_idx = row - 1  # Convert to 0-based
        barcode_idx = col - 1  # Convert to 0-based
        if gene_idx < len(features_df) and barcode_idx < len(barcodes_list):
            matrix_dict[(gene_idx, barcode_idx)] = val
    
    return features_df, barcodes_list, matrix_dict


def filter_mex_features(features_df, matrix_dict, keep_types):
    """Filter features/matrix to selected feature types. Returns (features_df, matrix_dict)."""
    if not keep_types:
        return features_df, matrix_dict
    keep_mask = features_df['type'].isin(keep_types)
    if keep_mask.all():
        return features_df, matrix_dict
    keep_indices = [idx for idx, keep in enumerate(keep_mask.values) if keep]
    old_to_new = {old_idx: new_idx for new_idx, old_idx in enumerate(keep_indices)}
    filtered_features = features_df.loc[keep_mask].reset_index(drop=True)
    filtered_matrix = {}
    for (gene_idx, barcode_idx), val in matrix_dict.items():
        new_idx = old_to_new.get(gene_idx)
        if new_idx is not None:
            filtered_matrix[(new_idx, barcode_idx)] = val
    return filtered_features, filtered_matrix


def build_gene_barcode_matrix(features_df, barcodes_list, matrix_dict, harmonized_bc_map):
    """
    Build a dense matrix: genes (rows) x harmonized barcodes (cols).
    Only includes barcodes that actually have counts in the matrix.
    Returns (gene_counts_df, harmonized_bc_list).
    """
    # First, find which harmonized barcodes actually have counts
    harmonized_bcs_with_counts = set()
    gene_bc_counts = defaultdict(lambda: defaultdict(int))
    
    for (gene_idx, orig_bc_idx), count in matrix_dict.items():
        if orig_bc_idx < len(barcodes_list):
            orig_bc = barcodes_list[orig_bc_idx]
            harmonized_bc = harmonized_bc_map[orig_bc]
            harmonized_bcs_with_counts.add(harmonized_bc)
            gene_bc_counts[gene_idx][harmonized_bc] += count
    
    # Only use barcodes that have counts
    harmonized_bc_list = sorted(harmonized_bcs_with_counts)
    harmonized_idx_map = {bc: idx for idx, bc in enumerate(harmonized_bc_list)}
    
    # Convert to DataFrame: rows = genes, cols = harmonized barcodes
    n_genes = len(features_df)
    n_harmonized_bcs = len(harmonized_bc_list)
    
    # Initialize matrix with zeros
    matrix_data = np.zeros((n_genes, n_harmonized_bcs), dtype=np.int32)
    
    for gene_idx in range(n_genes):
        for harmonized_bc, count in gene_bc_counts[gene_idx].items():
            harmonized_idx = harmonized_idx_map[harmonized_bc]
            matrix_data[gene_idx, harmonized_idx] = count
    
    # Create DataFrame with gene IDs/names as index
    gene_counts_df = pd.DataFrame(
        matrix_data,
        index=features_df['id'].values,
        columns=harmonized_bc_list
    )
    
    return gene_counts_df, harmonized_bc_list


def compare_mex(a375_dir, star_dir, min_cells_pct, min_counts, feature_types):
    """Compare A375 Cell Ranger MEX with STAR MEX."""
    print("=" * 70)
    print("A375 Cell Ranger vs STAR MEX Comparison")
    print("=" * 70)
    print()
    
    # Load both MEX directories
    print("Loading A375 Cell Ranger MEX...")
    a375_features, a375_barcodes, a375_matrix = load_mex(a375_dir)
    print(f"  Features: {len(a375_features)}")
    print(f"  Barcodes: {len(a375_barcodes)}")
    print(f"  Matrix entries: {len(a375_matrix)}")

    if feature_types:
        a375_features, a375_matrix = filter_mex_features(
            a375_features, a375_matrix, feature_types
        )
        print(f"  Features after filter: {len(a375_features)}")
        print(f"  Matrix entries after filter: {len(a375_matrix)}")
    
    print("\nLoading STAR MEX...")
    star_features, star_barcodes, star_matrix = load_mex(star_dir)
    print(f"  Features: {len(star_features)}")
    print(f"  Barcodes: {len(star_barcodes)}")
    print(f"  Matrix entries: {len(star_matrix)}")

    if feature_types:
        star_features, star_matrix = filter_mex_features(
            star_features, star_matrix, feature_types
        )
        print(f"  Features after filter: {len(star_features)}")
        print(f"  Matrix entries after filter: {len(star_matrix)}")
    print()
    
    # Harmonize barcodes
    print("Harmonizing barcodes...")
    a375_harmonized = {bc: harmonize_barcode(bc) for bc in a375_barcodes}
    star_harmonized = {bc: harmonize_barcode(bc) for bc in star_barcodes}
    
    a375_harmonized_set = set(a375_harmonized.values())
    star_harmonized_set = set(star_harmonized.values())
    common_barcodes = a375_harmonized_set & star_harmonized_set
    
    print(f"  A375 unique harmonized barcodes: {len(a375_harmonized_set)}")
    print(f"  STAR unique harmonized barcodes: {len(star_harmonized_set)}")
    print(f"  Common harmonized barcodes: {len(common_barcodes)}")
    print()
    
    # Build gene-barcode matrices (only for common barcodes)
    print("Building gene-barcode matrices...")
    a375_gene_df, a375_harm_bcs = build_gene_barcode_matrix(
        a375_features, a375_barcodes, a375_matrix, a375_harmonized
    )
    star_gene_df, star_harm_bcs = build_gene_barcode_matrix(
        star_features, star_barcodes, star_matrix, star_harmonized
    )
    
    # Filter to common barcodes only (and drop barcodes absent from either matrix)
    common_bc_list = sorted(common_barcodes)
    a375_cols = [bc for bc in common_bc_list if bc in a375_gene_df.columns]
    star_cols = [bc for bc in common_bc_list if bc in star_gene_df.columns]
    common_bc_list = sorted(set(a375_cols) & set(star_cols))
    a375_gene_df = a375_gene_df[common_bc_list]
    star_gene_df = star_gene_df[common_bc_list]
    
    print(f"  A375 matrix shape: {a375_gene_df.shape}")
    print(f"  STAR matrix shape: {star_gene_df.shape}")
    print()
    
    # Find common genes
    a375_genes = set(a375_gene_df.index)
    star_genes = set(star_gene_df.index)
    common_genes = sorted(a375_genes & star_genes)
    
    print("Gene Statistics:")
    print(f"  A375 unique genes: {len(a375_genes)}")
    print(f"  STAR unique genes: {len(star_genes)}")
    print(f"  Common genes: {len(common_genes)}")
    print()
    
    # Align matrices to common genes and barcodes
    a375_aligned = a375_gene_df.loc[common_genes, common_bc_list].fillna(0).astype(int)
    star_aligned = star_gene_df.loc[common_genes, common_bc_list].fillna(0).astype(int)
    
    # Ensure same column order
    a375_aligned = a375_aligned[common_bc_list]
    star_aligned = star_aligned[common_bc_list]
    
    # Calculate per-gene sums across all barcodes
    a375_gene_sums = a375_aligned.sum(axis=1)
    star_gene_sums = star_aligned.sum(axis=1)

    # Filter genes by minimum cells and counts in both datasets
    min_cells = max(1, int(np.ceil(len(common_bc_list) * min_cells_pct)))
    a375_gene_cells = (a375_aligned > 0).sum(axis=1)
    star_gene_cells = (star_aligned > 0).sum(axis=1)
    filter_mask = (
        (a375_gene_sums >= min_counts)
        & (star_gene_sums >= min_counts)
        & (a375_gene_cells >= min_cells)
        & (star_gene_cells >= min_cells)
    )
    filtered_genes = [gene for gene in common_genes if filter_mask.loc[gene]]

    print("Filtering for correlations:")
    print(f"  min_cells_per_gene: {min_cells} ({min_cells_pct:.2%} of {len(common_bc_list)} cells)")
    print(f"  min_counts_per_gene: {min_counts}")
    print(f"  Genes passing filter: {len(filtered_genes)}")
    print()
    
    # Calculate correlations
    print("Correlation Analysis:")
    print()
    
    # Pearson: filtered genes
    if len(filtered_genes) > 1:
        pearson_r, pearson_p = pearsonr(a375_gene_sums[filtered_genes], star_gene_sums[filtered_genes])
        print(f"Pearson (filtered {len(filtered_genes)} genes):")
        print(f"  r = {pearson_r:.6f}")
        print(f"  p = {pearson_p:.2e}")
    else:
        print("Pearson: Insufficient data (need >1 filtered gene)")
        pearson_r = np.nan
    
    print()
    
    # Spearman: filtered genes
    if len(filtered_genes) > 1:
        spearman_r, spearman_p = spearmanr(
            a375_gene_sums[filtered_genes],
            star_gene_sums[filtered_genes]
        )
        print(f"Spearman (filtered {len(filtered_genes)} genes):")
        print(f"  ρ = {spearman_r:.6f}")
        print(f"  p = {spearman_p:.2e}")
    else:
        print("Spearman: Insufficient data (need >1 filtered gene)")
        spearman_r = np.nan
    
    print()
    
    # Summary
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"Common harmonized barcodes: {len(common_barcodes)}")
    print(f"Common genes: {len(common_genes)}")
    print(f"Filtered genes: {len(filtered_genes)}")
    if not np.isnan(pearson_r):
        print(f"Pearson correlation: {pearson_r:.6f}")
    if not np.isnan(spearman_r):
        print(f"Spearman correlation: {spearman_r:.6f}")
    print()
    
    return {
        'common_barcodes': len(common_barcodes),
        'common_genes': len(common_genes),
        'filtered_genes': len(filtered_genes),
        'pearson_r': pearson_r,
        'spearman_r': spearman_r
    }


def main():
    parser = argparse.ArgumentParser(
        description='Compare A375 Cell Ranger MEX with STAR MEX output'
    )
    parser.add_argument('a375_dir', help='A375 Cell Ranger MEX directory')
    parser.add_argument('star_dir', help='STAR MEX directory')
    parser.add_argument(
        '--min-cells-pct',
        type=float,
        default=0.01,
        help='Minimum fraction of common cells expressing a gene (default: 0.01)',
    )
    parser.add_argument(
        '--min-counts',
        type=int,
        default=20,
        help='Minimum total counts per gene in both datasets (default: 20)',
    )
    parser.add_argument(
        '--feature-types',
        type=str,
        default='Gene Expression',
        help='Comma-separated feature types to include (default: Gene Expression). '
             'Use "all" to disable filtering.',
    )
    args = parser.parse_args()
    
    if not os.path.isdir(args.a375_dir):
        print(f"Error: A375 directory not found: {args.a375_dir}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.isdir(args.star_dir):
        print(f"Error: STAR directory not found: {args.star_dir}", file=sys.stderr)
        sys.exit(1)
    
    feature_types = None
    if args.feature_types and args.feature_types.lower() not in ('all', '*', 'any', 'none'):
        feature_types = [ft.strip() for ft in args.feature_types.split(',') if ft.strip()]

    try:
        results = compare_mex(
            args.a375_dir,
            args.star_dir,
            args.min_cells_pct,
            args.min_counts,
            feature_types,
        )
        sys.exit(0)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
