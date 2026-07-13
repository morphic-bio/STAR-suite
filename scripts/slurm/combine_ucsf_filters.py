#!/usr/bin/env python3
"""Combine STAR cell calls, doublets, and baseline QC in one AnnData file."""

from __future__ import annotations

import argparse
from pathlib import Path


DEFAULT_MITO_GENES = (
    "ENSG00000198888",
    "ENSG00000198763",
    "ENSG00000198804",
    "ENSG00000198712",
    "ENSG00000228253",
    "ENSG00000198899",
    "ENSG00000198938",
    "ENSG00000198840",
    "ENSG00000212907",
    "ENSG00000198886",
    "ENSG00000198786",
    "ENSG00000198695",
    "ENSG00000198727",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add STAR cell, doublet, and baseline QC masks to an AnnData file."
    )
    parser.add_argument("--input-file", required=True, type=Path)
    parser.add_argument("--non-empty-barcodes", required=True, type=Path)
    parser.add_argument("--doublet-barcodes", required=True, type=Path)
    parser.add_argument("--doublet-scores", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mito-genes", type=Path)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-genes", type=int, default=2500)
    parser.add_argument("--mt-pct-cutoff", type=float, default=5.0)
    return parser.parse_args()


def _read_lines(path: Path) -> set[str]:
    if not path.is_file():
        raise SystemExit(f"missing barcode file: {path}")
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def main() -> int:
    args = parse_args()
    if args.min_genes < 0 or args.max_genes < args.min_genes:
        raise SystemExit("invalid min/max gene thresholds")

    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    args.output_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.input_file)
    non_empty = _read_lines(args.non_empty_barcodes)
    doublets = _read_lines(args.doublet_barcodes)
    if not args.doublet_scores.is_file():
        raise SystemExit(f"missing doublet score file: {args.doublet_scores}")
    doublet_scores = pd.read_csv(
        args.doublet_scores,
        sep="\t",
        dtype={"Barcode": "string", "Classification": "string"},
    )
    required_score_columns = {"Barcode", "Classification", "Score"}
    if not required_score_columns.issubset(doublet_scores.columns):
        missing = sorted(required_score_columns - set(doublet_scores.columns))
        raise SystemExit(f"doublet score file is missing columns: {', '.join(missing)}")
    if doublet_scores["Barcode"].duplicated().any():
        raise SystemExit("doublet score file contains duplicate barcodes")
    doublet_scores["Score"] = pd.to_numeric(doublet_scores["Score"], errors="coerce")
    score_barcodes = set(doublet_scores["Barcode"].astype(str))
    obs_barcode_set = set(adata.obs_names.astype(str))
    unknown_non_empty = non_empty - obs_barcode_set
    if unknown_non_empty:
        raise SystemExit("non-empty barcode file contains barcodes absent from the AnnData")
    if not doublets.issubset(non_empty):
        raise SystemExit("doublet barcode file is not a subset of the non-empty barcode file")
    if score_barcodes != non_empty:
        raise SystemExit("doublet score file must cover exactly the non-empty barcode set")
    score_doublets = set(
        doublet_scores.loc[doublet_scores["Classification"] == "doublet", "Barcode"].astype(str)
    )
    if score_doublets != doublets:
        raise SystemExit("doublet score classifications do not match the doublet barcode file")
    if not np.isfinite(doublet_scores["Score"].to_numpy(dtype=float)).all():
        raise SystemExit("doublet score file contains non-finite scores")
    score_by_barcode = doublet_scores.set_index("Barcode")["Score"]
    obs_names = adata.obs_names.astype(str)

    adata.obs["non_empty"] = np.isin(obs_names, list(non_empty))
    adata.obs["doublet"] = np.isin(obs_names, list(doublets))
    adata.obs["doublet_scores"] = obs_names.map(score_by_barcode)
    if sparse.issparse(adata.X):
        n_genes = np.asarray(adata.X.getnnz(axis=1)).ravel()
        total_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    else:
        dense_counts = np.asarray(adata.X)
        n_genes = np.count_nonzero(dense_counts, axis=1)
        total_counts = dense_counts.sum(axis=1)
    adata.obs["n_genes"] = n_genes
    adata.obs["total_counts"] = total_counts

    mito_genes = (
        _read_lines(args.mito_genes) if args.mito_genes is not None else set(DEFAULT_MITO_GENES)
    )
    present_mito = [name for name in adata.var_names.astype(str) if name in mito_genes]
    if present_mito:
        mt_matrix = adata[:, present_mito].X
        mt_counts = np.asarray(mt_matrix.sum(axis=1)).ravel()
    else:
        mt_counts = np.zeros(adata.n_obs, dtype=float)
    adata.obs["mt_counts"] = mt_counts
    denominator = np.maximum(total_counts, 1)
    adata.obs["mt_pct"] = mt_counts / denominator * 100.0

    baseline_filter = (
        (n_genes >= args.min_genes)
        & (n_genes <= args.max_genes)
        & (adata.obs["mt_pct"].to_numpy() <= args.mt_pct_cutoff)
    )
    singlet = adata.obs["non_empty"].to_numpy() & ~adata.obs["doublet"].to_numpy()
    adata.obs["filter"] = baseline_filter
    adata.obs["singlet"] = singlet
    adata.obs["singlet_filtered"] = singlet & baseline_filter

    unfiltered_path = args.output_dir / "unfiltered_counts.h5ad"
    filtered_path = args.output_dir / "filtered_counts.h5ad"
    adata.write_h5ad(unfiltered_path)
    adata[singlet].copy().write_h5ad(filtered_path)
    print(f"Non-empty cells: {int(adata.obs['non_empty'].sum())}")
    print(f"Doublets: {int(adata.obs['doublet'].sum())}")
    print(f"Baseline QC pass: {int(adata.obs['filter'].sum())}")
    print(f"Wrote {unfiltered_path}")
    print(f"Wrote {filtered_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
