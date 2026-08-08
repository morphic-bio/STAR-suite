#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp


def ensure_layer(source: ad.AnnData, layer_name: str):
    if layer_name not in source.layers:
        raise KeyError(f"Missing layer '{layer_name}' in source AnnData")
    return source.layers[layer_name]


def subset_layer_by_obs(source: ad.AnnData, target: ad.AnnData, layer_name: str):
    if source.n_vars != target.n_vars:
        raise ValueError(
            f"Source/target var dimensions differ: {source.n_vars} vs {target.n_vars}"
        )
    if not np.array_equal(np.asarray(source.var_names), np.asarray(target.var_names)):
        raise ValueError("Source and target var_names differ; refusing to propagate layer")

    obs_name_to_idx = {name: idx for idx, name in enumerate(source.obs_names.astype(str))}
    target_obs = target.obs_names.astype(str)
    missing = [name for name in target_obs if name not in obs_name_to_idx]
    if missing:
        preview = ", ".join(missing[:5])
        raise KeyError(
            f"{len(missing)} target barcodes missing from source layer surface; "
            f"examples: {preview}"
        )

    row_idx = np.fromiter((obs_name_to_idx[name] for name in target_obs), dtype=np.int64)
    layer = ensure_layer(source, layer_name)
    subset = layer[row_idx, :]
    if not sp.issparse(subset):
        subset = sp.csr_matrix(subset)
    else:
        subset = subset.tocsr()
    return subset


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Propagate an AnnData layer from a full source object into a target subset "
            "by matching obs_names and preserving target order."
        )
    )
    parser.add_argument("--source-h5ad", required=True)
    parser.add_argument("--target-h5ad", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--layer-name", required=True)
    args = parser.parse_args()

    source_h5ad = Path(args.source_h5ad).resolve()
    target_h5ad = Path(args.target_h5ad).resolve()
    output_h5ad = Path(args.output_h5ad).resolve()

    source = ad.read_h5ad(source_h5ad)
    target = ad.read_h5ad(target_h5ad)
    target.layers[args.layer_name] = subset_layer_by_obs(source, target, args.layer_name)
    target.write_h5ad(output_h5ad)

    print(f"Source: {source_h5ad} ({source.n_obs} x {source.n_vars})")
    print(f"Target: {target_h5ad} ({target.n_obs} x {target.n_vars})")
    print(f"Layer: {args.layer_name}")
    print(f"Wrote {output_h5ad}")


if __name__ == "__main__":
    main()
