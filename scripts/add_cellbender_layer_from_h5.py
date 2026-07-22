#!/usr/bin/env python3

import argparse

import anndata as ad
import numpy as np
import scipy.sparse as sp
from cellbender.remove_background.downstream import anndata_from_h5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add a CellBender-denoised sparse layer to an AnnData object."
    )
    parser.add_argument("--cellbender-h5", required=True, help="CellBender .h5 output")
    parser.add_argument("--input-h5ad", required=True, help="Input AnnData file")
    parser.add_argument("--output-h5ad", required=True, help="Output AnnData file")
    parser.add_argument("--layer-name", default="denoised", help="Layer name to write")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    original = ad.read_h5ad(args.input_h5ad)
    cb = anndata_from_h5(args.cellbender_h5, analyzed_barcodes_only=True)

    original_barcodes = {bc: idx for idx, bc in enumerate(original.obs_names)}
    cb_rows = []
    target_rows = []
    for cb_idx, bc in enumerate(cb.obs_names):
      idx = original_barcodes.get(bc)
      if idx is not None:
        cb_rows.append(cb_idx)
        target_rows.append(idx)

    if len(cb_rows) == 0:
        raise SystemExit("No overlapping barcodes between CellBender output and input h5ad")

    cb_matrix = cb.X.tocoo()
    target_rows = np.asarray(target_rows, dtype=np.int64)
    mapped_rows = target_rows[cb_matrix.row]
    layer = sp.csr_matrix(
        (cb_matrix.data, (mapped_rows, cb_matrix.col)),
        shape=original.shape,
        dtype=cb_matrix.data.dtype,
    )
    original.layers[args.layer_name] = layer
    original.write_h5ad(args.output_h5ad)


if __name__ == "__main__":
    main()
