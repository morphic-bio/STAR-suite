#!/usr/bin/env python3
"""
Project STAR-Flex full-align and no-align outputs into a CR9 reference embedding.

This reproduces the SC2300771 parity figures used in the README:
1. build a CR9 reference embedding and Leiden labels
2. project STAR full-align and no-align into that fixed CR9 space with ingest
3. write projected AnnData, PNGs, and projected-label agreement tables
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

matplotlib.use("Agg")

STAR_TAG_TO_SAMPLE = {
    "BC004": "WT-Day-7",
    "BC006": "PAX6-PTC-D9-Day7",
    "BC007": "WT-Day-8",
    "BC008": "PAX6-PTC-D9-Day8",
}
CR9_SAMPLES = ["WT-Day-7", "PAX6-PTC-D9-Day7", "WT-Day-8", "PAX6-PTC-D9-Day8"]


def read_star_sample(sample_tag: str, root: Path, pipeline: str) -> ad.AnnData:
    sample_dir = root / sample_tag / "Gene" / "filtered"
    matrix = mmread(sample_dir / "matrix.mtx").tocsr().transpose().tocsr()
    barcodes = pd.read_csv(sample_dir / "barcodes.tsv", sep="\t", header=None)[0].astype(str)
    features = pd.read_csv(
        sample_dir / "features.tsv",
        sep="\t",
        header=None,
        names=["gene_id", "gene_name", "feature_type"],
    )

    obs = pd.DataFrame(index=[f"{pipeline}:{sample_tag}:{bc}" for bc in barcodes])
    obs["sample_id"] = sample_tag
    obs["sample_name"] = STAR_TAG_TO_SAMPLE[sample_tag]
    obs["pipeline"] = pipeline

    var = pd.DataFrame(index=features["gene_id"].astype(str).tolist())
    var["gene_name"] = features["gene_name"].astype(str).values
    var["feature_type"] = features["feature_type"].astype(str).values
    return ad.AnnData(X=matrix, obs=obs, var=var)


def read_cr9_sample(sample_name: str, root: Path) -> ad.AnnData:
    h5_path = root / sample_name / "count" / "sample_filtered_feature_bc_matrix.h5"
    adata = sc.read_10x_h5(h5_path)
    adata.var_names = adata.var["gene_ids"].astype(str)
    adata.var_names_make_unique()
    adata.obs_names = [f"cr9:{sample_name}:{bc}" for bc in adata.obs_names.astype(str)]
    adata.obs["sample_id"] = sample_name
    adata.obs["sample_name"] = sample_name
    adata.obs["pipeline"] = "cr9"
    if "feature_types" not in adata.var.columns:
        adata.var["feature_types"] = "Gene Expression"
    if "gene_name" not in adata.var.columns:
        if "gene_symbols" in adata.var.columns:
            adata.var["gene_name"] = adata.var["gene_symbols"].astype(str)
        else:
            adata.var["gene_name"] = adata.var_names.astype(str)
    return adata


def build_cr9_reference(adata: ad.AnnData) -> ad.AnnData:
    adata = adata.copy()
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata, random_state=0)
    sc.tl.leiden(adata)
    return adata


def prepare_query(adata: ad.AnnData, ref_hvgs: pd.Index) -> ad.AnnData:
    adata = adata.copy()
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    common = ref_hvgs.intersection(adata.var_names)
    adata = adata[:, common].copy()
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    return adata[:, ref_hvgs].copy()


def write_projection_summary(outdir: Path, fullalign: ad.AnnData, noalign: ad.AnnData) -> None:
    fa_obs = fullalign.obs[["sample_name", "leiden"]].copy()
    na_obs = noalign.obs[["sample_name", "leiden"]].copy()
    fa_obs["cell_key"] = [":".join(x.split(":")[1:]) for x in fa_obs.index]
    na_obs["cell_key"] = [":".join(x.split(":")[1:]) for x in na_obs.index]
    fa_obs = fa_obs.set_index("cell_key")
    na_obs = na_obs.set_index("cell_key")
    shared = fa_obs.index.intersection(na_obs.index)

    ari = adjusted_rand_score(fa_obs.loc[shared, "leiden"], na_obs.loc[shared, "leiden"])
    nmi = normalized_mutual_info_score(fa_obs.loc[shared, "leiden"], na_obs.loc[shared, "leiden"])

    pd.DataFrame(
        [
            {
                "shared_cells": int(len(shared)),
                "ari_projected_leiden": float(ari),
                "nmi_projected_leiden": float(nmi),
            }
        ]
    ).to_csv(outdir / "fullalign_vs_noalign_projected_leiden_agreement.tsv", sep="\t", index=False)
    pd.crosstab(fa_obs.loc[shared, "leiden"], na_obs.loc[shared, "leiden"]).to_csv(
        outdir / "fullalign_vs_noalign_projected_leiden_crosstab.tsv",
        sep="\t",
    )


def plot_umap(adata: ad.AnnData, outdir: Path, stem: str) -> None:
    sc.settings.figdir = str(outdir)
    sc.settings.set_figure_params(dpi_save=200)
    sc.pl.umap(adata, color=["leiden", "sample_name"], show=False, save=f"_{stem}.png")
    plt.close("all")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full-root",
        default="/storage/flex_fullalign_2024ref_20260325_020044/per_sample_filtered",
        help="STAR-Flex full-align filtered per-sample root",
    )
    parser.add_argument(
        "--noalign-root",
        default="/storage/flex_noalign_full_20260325_022459/per_sample_filtered",
        help="STAR-Flex no-align filtered per-sample root",
    )
    parser.add_argument(
        "--cr9-root",
        default="/mnt/pikachu/benchmark_cr9_flex_full/outs/per_sample_outs",
        help="CellRanger 9 per-sample output root",
    )
    parser.add_argument(
        "--outdir",
        default="comparisons/flex_leiden_cr9_projection_20260330",
        help="Output directory",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    fullalign = ad.concat(
        [read_star_sample(tag, Path(args.full_root), "fullalign") for tag in STAR_TAG_TO_SAMPLE],
        join="outer",
        merge="same",
    )
    noalign = ad.concat(
        [read_star_sample(tag, Path(args.noalign_root), "noalign") for tag in STAR_TAG_TO_SAMPLE],
        join="outer",
        merge="same",
    )
    cr9 = ad.concat(
        [read_cr9_sample(sample, Path(args.cr9_root)) for sample in CR9_SAMPLES],
        join="outer",
        merge="same",
    )

    ref = build_cr9_reference(cr9)
    ref.write_h5ad(outdir / "sc2300771_cr9_reference.h5ad")
    plot_umap(ref, outdir, "sc2300771_cr9_reference")

    results = []
    for label, query_raw in [("fullalign", fullalign), ("noalign", noalign)]:
        query = prepare_query(query_raw, ref.var_names)
        sc.tl.ingest(query, ref, obs="leiden")
        query.write_h5ad(outdir / f"sc2300771_{label}_projected_to_cr9.h5ad")
        plot_umap(query, outdir, f"sc2300771_{label}_projected_to_cr9")
        query.obs["projected_leiden"] = query.obs["leiden"].astype(str)
        query.obs[["sample_name", "projected_leiden"]].value_counts().rename("n_cells").reset_index().to_csv(
            outdir / f"sc2300771_{label}_projected_leiden_by_sample.tsv",
            sep="\t",
            index=False,
        )
        results.append(
            {
                "pipeline": label,
                "cells_projected": int(query.n_obs),
                "genes_projected": int(query.n_vars),
                "projected_clusters_used": int(query.obs["projected_leiden"].nunique()),
            }
        )

    pd.DataFrame(results).to_csv(outdir / "projection_summary.tsv", sep="\t", index=False)

    write_projection_summary(
        outdir,
        ad.read_h5ad(outdir / "sc2300771_fullalign_projected_to_cr9.h5ad", backed="r"),
        ad.read_h5ad(outdir / "sc2300771_noalign_projected_to_cr9.h5ad", backed="r"),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
