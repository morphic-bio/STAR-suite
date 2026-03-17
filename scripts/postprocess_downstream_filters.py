#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np


def require_obs_column(adata: ad.AnnData, column: str) -> np.ndarray:
    if column not in adata.obs.columns:
        raise KeyError(f"Missing required obs column: {column}")
    return adata.obs[column].astype(bool).to_numpy()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Rewrite downstream filtered outputs so QC-only and default singlet-filtered "
            "views are both preserved."
        )
    )
    parser.add_argument("--unfiltered-h5ad", required=True)
    parser.add_argument("--qc-output-h5ad", required=True)
    parser.add_argument("--default-singlet-output-h5ad", required=True)
    args = parser.parse_args()

    unfiltered_h5ad = Path(args.unfiltered_h5ad).resolve()
    qc_output_h5ad = Path(args.qc_output_h5ad).resolve()
    default_singlet_output_h5ad = Path(args.default_singlet_output_h5ad).resolve()

    adata = ad.read_h5ad(unfiltered_h5ad)
    qc_mask = require_obs_column(adata, "filter")
    singlet_mask = require_obs_column(adata, "singlet_filtered")
    non_empty_mask = require_obs_column(adata, "non_empty")

    # Restrict exported QC views to STAR-called cells only. The upstream
    # combineFilters.py computes QC metrics on the full raw-backed object.
    qc_filtered = adata[qc_mask & non_empty_mask].copy()
    qc_filtered.write_h5ad(qc_output_h5ad)

    singlet_filtered = adata[singlet_mask].copy()
    singlet_filtered.write_h5ad(default_singlet_output_h5ad)

    print(f"Unfiltered cells: {adata.n_obs}")
    print(f"QC-only filtered cells: {qc_filtered.n_obs}")
    print(f"Default singlet-filtered cells: {singlet_filtered.n_obs}")
    print(f"Wrote {qc_output_h5ad}")
    print(f"Wrote {default_singlet_output_h5ad}")


if __name__ == "__main__":
    main()
