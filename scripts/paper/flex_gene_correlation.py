#!/usr/bin/env python3
"""Gene-level concordance between per-sample Flex MEX outputs (STAR Suite vs Cell Ranger).

For each sample pair, per-gene count totals are summed over the barcodes called by both tools and compared
over every gene present in both feature lists, zero-count genes included:
  - Spearman correlation (tied ranks averaged)
  - Pearson correlation on the raw (untransformed) totals
The script prints one line per sample and the mean over samples, which is how the manuscript reports Flex
(Supplementary Table S3). The definitions match `spearman_all_genes` / `pearson_all_genes` in
scripts/report_additional_parity_metrics.py.

Barcodes: a trailing "-1" (Cell Ranger GEM-well suffix) is stripped. Within a per-sample directory the barcode
identifies the cell; for pooled-tag samples Cell Ranger and STAR Suite both write CB16 + tag barcodes.
Only "Gene Expression" features are used when features.tsv has a type column.

Usage:
  flex_gene_correlation.py LABEL CR_MEX STAR_MEX [CR_MEX STAR_MEX ...]

  CR_MEX   Cell Ranger per-sample matrix, e.g. outs/per_sample_outs/<sample>/count/sample_filtered_feature_bc_matrix
  STAR_MEX STAR Suite per-sample matrix, e.g. <prefix>/per_sample/<tag>/Gene/filtered
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


def mex_file(d, name):
    for n in (name, name + ".gz"):
        path = os.path.join(d, n)
        if os.path.exists(path):
            return path
    raise FileNotFoundError(os.path.join(d, name))


def load_mex(d):
    feats = pd.read_csv(mex_file(d, "features.tsv"), sep="\t", header=None)
    keep = feats[2] == "Gene Expression" if feats.shape[1] > 2 else pd.Series(True, index=feats.index)
    rows = {i + 1: gene for i, gene in enumerate(feats[0]) if keep[i]}
    path = mex_file(d, "barcodes.tsv")
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        barcodes = [line.strip()[:-2] if line.strip().endswith("-1") else line.strip() for line in fh]
    entries = pd.read_csv(mex_file(d, "matrix.mtx"), sep=" ", comment="%", header=None)
    entries = entries.iloc[1:].astype({0: np.int64, 1: np.int64, 2: np.float64})  # drop the dimension line
    return rows, barcodes, entries


def gene_totals(rows, barcodes, entries, allowed):
    in_set = np.array([b in allowed for b in barcodes])[entries[1].values - 1]
    sums = np.bincount(entries[0].values[in_set] - 1, weights=entries[2].values[in_set], minlength=max(rows))
    return {gene: sums[row - 1] for row, gene in rows.items()}


def main(argv):
    if len(argv) < 4 or len(argv) % 2:
        sys.exit(__doc__)
    label, pairs = argv[1], argv[2:]
    spearmans, pearsons = [], []
    for cr_dir, star_dir in zip(pairs[0::2], pairs[1::2]):
        cr, star = load_mex(cr_dir), load_mex(star_dir)
        shared = set(cr[1]) & set(star[1])
        cr_tot, star_tot = gene_totals(*cr, shared), gene_totals(*star, shared)
        genes = sorted(set(cr_tot) & set(star_tot))
        a = np.array([cr_tot[g] for g in genes])
        b = np.array([star_tot[g] for g in genes])
        rho, r = spearmanr(a, b).statistic, pearsonr(a, b).statistic
        spearmans.append(rho)
        pearsons.append(r)
        print(f"{label}\t{star_dir}\tcr_cells={len(cr[1])}\tstar_cells={len(star[1])}\tshared={len(shared)}"
              f"\tgenes={len(genes)}\tspearman={rho:.6f}\tpearson={r:.6f}")
    print(f"{label}\tMEAN\tsamples={len(spearmans)}\tspearman={np.mean(spearmans):.6f}\tpearson={np.mean(pearsons):.6f}")


if __name__ == "__main__":
    main(sys.argv)
