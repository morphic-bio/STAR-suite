#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csr_matrix, issparse


def read_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle if line.strip()]


def read_feature_ids(path: Path) -> list[str]:
    feature_rows = [row.split("\t") for row in read_lines(path)]
    return [row[0] for row in feature_rows]


def read_sparse_matrix(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        matrix = mmread(handle)
    if not issparse(matrix):
        raise TypeError(f"Expected sparse matrix input: {path}")
    return matrix.tocsr()


def strip_barcode_suffix(barcode: str) -> str:
    if barcode.endswith("-1"):
        return barcode[:-2]
    return barcode


def ensure_csr(matrix):
    if issparse(matrix):
        return matrix.tocsr()
    return csr_matrix(matrix)


def build_indexer(keys: list[str]) -> dict[str, int]:
    return {key: idx for idx, key in enumerate(keys)}


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build counts.h5ad from GeneFull Cell Ranger-style MEX plus packaged raw velocyto layers."
        )
    )
    parser.add_argument("--run-dir", required=True, help="STAR/CR-compat run directory containing outs/")
    parser.add_argument("--output-h5ad", required=True, help="Output AnnData path")
    parser.add_argument(
        "--feature-raw-dir",
        help="Override raw GeneFull MEX directory (default: <run-dir>/outs/raw_feature_bc_matrix)",
    )
    parser.add_argument(
        "--feature-filtered-dir",
        help="Override filtered GeneFull MEX directory (default: <run-dir>/outs/filtered_feature_bc_matrix)",
    )
    parser.add_argument(
        "--velocyto-raw-dir",
        help="Override packaged raw velocyto MEX directory (default: <run-dir>/outs/raw_velocyto_feature_bc_matrix)",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    feature_raw_dir = Path(args.feature_raw_dir).resolve() if args.feature_raw_dir else run_dir / "outs" / "raw_feature_bc_matrix"
    feature_filtered_dir = (
        Path(args.feature_filtered_dir).resolve()
        if args.feature_filtered_dir
        else run_dir / "outs" / "filtered_feature_bc_matrix"
    )
    velocyto_raw_dir = (
        Path(args.velocyto_raw_dir).resolve()
        if args.velocyto_raw_dir
        else run_dir / "outs" / "raw_velocyto_feature_bc_matrix"
    )
    output_h5ad = Path(args.output_h5ad).resolve()
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)

    for path in [feature_raw_dir, feature_filtered_dir, velocyto_raw_dir]:
        if not path.exists():
            raise FileNotFoundError(f"Missing required directory: {path}")

    adata_raw = sc.read_10x_mtx(feature_raw_dir, var_names="gene_ids", cache=False)
    adata_filtered = sc.read_10x_mtx(feature_filtered_dir, var_names="gene_ids", cache=False)

    adata_raw.X = ensure_csr(adata_raw.X)
    adata_raw.obs["is_cell"] = adata_raw.obs_names.isin(adata_filtered.obs_names)
    adata_raw.obs["filter"] = adata_raw.obs["is_cell"]

    velocyto_barcodes = read_lines(velocyto_raw_dir / "barcodes.tsv.gz")
    velocyto_features = read_feature_ids(velocyto_raw_dir / "features.tsv.gz")

    raw_feature_ids = adata_raw.var_names.tolist()
    if raw_feature_ids != velocyto_features:
        velocyto_feature_index = build_indexer(velocyto_features)
        missing_genes = [gene_id for gene_id in raw_feature_ids if gene_id not in velocyto_feature_index]
        if missing_genes:
            preview = ", ".join(missing_genes[:5])
            raise ValueError(f"Velocyto features missing GeneFull IDs: {preview}")
        row_indices = np.array([velocyto_feature_index[gene_id] for gene_id in raw_feature_ids], dtype=np.int64)
    else:
        row_indices = np.arange(len(raw_feature_ids), dtype=np.int64)

    velocyto_barcode_keys = [strip_barcode_suffix(barcode) for barcode in velocyto_barcodes]
    if len(set(velocyto_barcode_keys)) != len(velocyto_barcode_keys):
        raise ValueError("Velocyto barcodes are not unique after -1 suffix normalization")
    velocyto_barcode_index = build_indexer(velocyto_barcode_keys)
    raw_barcodes = [strip_barcode_suffix(barcode) for barcode in adata_raw.obs_names.tolist()]
    missing_barcodes = [barcode for barcode in raw_barcodes if barcode not in velocyto_barcode_index]
    if missing_barcodes:
        preview = ", ".join(missing_barcodes[:5])
        raise ValueError(f"Velocyto barcodes missing GeneFull raw barcodes: {preview}")
    col_indices = np.array([velocyto_barcode_index[barcode] for barcode in raw_barcodes], dtype=np.int64)

    for layer_name in ["spliced", "unspliced", "ambiguous"]:
        matrix = read_sparse_matrix(velocyto_raw_dir / f"{layer_name}.mtx.gz")
        subset = matrix[row_indices, :][:, col_indices].T.tocsr()
        if subset.shape != adata_raw.shape:
            raise ValueError(
                f"{layer_name} shape mismatch after subsetting: {subset.shape} vs {adata_raw.shape}"
            )
        adata_raw.layers[layer_name] = subset

    adata_raw.uns["gene_expression_source"] = str(feature_raw_dir)
    adata_raw.uns["velocyto_source"] = str(velocyto_raw_dir)
    adata_raw.uns["gene_expression_feature_kind"] = "GeneFull"

    adata_raw.write_h5ad(output_h5ad)
    print(f"Wrote {output_h5ad}")
    print(adata_raw)


if __name__ == "__main__":
    main()
