#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp


def read_barcode_set(path: Path) -> set[str]:
    with open(path, "r", encoding="utf-8") as handle:
        return {line.strip() for line in handle if line.strip()}


def compute_n_genes(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return np.asarray((matrix > 0).sum(axis=1)).ravel()
    return np.asarray((matrix > 0).sum(axis=1)).ravel()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute adaptive max_genes threshold from STAR-called singlets."
    )
    parser.add_argument("--counts-h5ad", required=True)
    parser.add_argument("--non-empty-barcodes", required=True)
    parser.add_argument("--doublet-barcodes", required=True)
    parser.add_argument("--min-genes", required=True, type=int)
    parser.add_argument("--n-mad", required=True, type=float)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()

    counts_h5ad = Path(args.counts_h5ad).resolve()
    non_empty_barcodes = read_barcode_set(Path(args.non_empty_barcodes).resolve())
    doublet_barcodes = read_barcode_set(Path(args.doublet_barcodes).resolve())

    adata = ad.read_h5ad(counts_h5ad)
    obs_names = np.asarray(adata.obs_names.astype(str))

    non_empty_mask = np.isin(obs_names, list(non_empty_barcodes))
    singlet_mask = non_empty_mask & ~np.isin(obs_names, list(doublet_barcodes))
    singlet_count = int(singlet_mask.sum())
    if singlet_count == 0:
        raise ValueError("No singlet STAR-called cells available for adaptive thresholding")

    n_genes = compute_n_genes(adata.X[singlet_mask])
    median = float(np.median(n_genes))
    mad = float(np.median(np.abs(n_genes - median)))
    raw_adaptive_max = int(median + args.n_mad * mad)
    effective_max = max(raw_adaptive_max, int(args.min_genes))

    stats = {
        "counts_h5ad": str(counts_h5ad),
        "star_called_cells": int(non_empty_mask.sum()),
        "singlet_cells": singlet_count,
        "n_mad": args.n_mad,
        "min_genes": int(args.min_genes),
        "median_n_genes": median,
        "mad_n_genes": mad,
        "raw_adaptive_max_genes": raw_adaptive_max,
        "effective_max_genes": effective_max,
        "max_genes_was_clamped": bool(effective_max != raw_adaptive_max),
    }

    output_json = Path(args.output_json).resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w", encoding="utf-8") as handle:
        json.dump(stats, handle, indent=2, sort_keys=True)

    print(json.dumps(stats, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
