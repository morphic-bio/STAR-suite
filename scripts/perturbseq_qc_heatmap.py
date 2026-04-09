#!/usr/bin/env python3
"""
Perturb-seq QC: cross-sample knockdown heatmap and single-sample waterfall.

For each sample, computes log2FC of each CRISPRi target gene in perturbed
cells vs NTC (Non-Targeting) control cells using CP10K + log1p normalization
and Wilcoxon rank-sum test with BH correction.

Outputs:
  - waterfall_EBs1_1.png   Single-sample bar chart sorted by log2FC
  - knockdown_heatmap.png  Cross-sample heatmap (rows=genes, cols=samples)
  - perturbseq_qc.tsv      Full results table
"""

import csv
import sys
import warnings
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from scipy import sparse
from scipy.stats import rankdata, mannwhitneyu

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

PROD_ROOT = Path(
    "/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_prod_globus_fixtruwl_20260330_225009"
)
DS_SUFFIX = "downstream_genefull_velocyto_cellbender_guarded_rerun"
FEAT_REF = Path(
    "/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv"
)
OUT_DIR = Path("/mnt/pikachu/ucsf_emptydrops_guarded_release_20260403/qc_plots")

SAMPLES = [
    "EBs1_1", "EBs1_2", "EBs1_3", "EBs1_4", "EBs1_5",
    "EBs2_1", "EBs2_2", "EBs2_3", "EBs2_4", "EBs2_5",
    "iPSC1_1", "iPSC1_2", "iPSC1_3",
    "iPSC2_1", "iPSC2_2", "iPSC2_3",
]

MIN_CELLS = 5        # minimum perturbed cells to test
FDR_THRESH = 0.05
LOG2FC_THRESH = -0.5  # negative = knockdown

# ---------------------------------------------------------------------------
# Guide-to-gene mapping
# ---------------------------------------------------------------------------

def load_guide_to_gene(feat_ref_path):
    """Return dict mapping guide name -> target gene name."""
    g2g = {}
    with open(feat_ref_path) as f:
        for row in csv.DictReader(f):
            g2g[row["name"]] = row["target_gene_name"]
    return g2g


def get_ntc_guides(guide2gene):
    """Return set of guide names that target Non-Targeting."""
    return {g for g, t in guide2gene.items() if t == "Non-Targeting"}


# ---------------------------------------------------------------------------
# BH correction
# ---------------------------------------------------------------------------

def bh_adjust(pvals):
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    adjusted = pvals * n / ranked
    # enforce monotonicity (working backwards from largest rank)
    adjusted = np.minimum.accumulate(adjusted[np.argsort(-ranked)])[np.argsort(np.argsort(-ranked))]
    return np.minimum(adjusted, 1.0)


# ---------------------------------------------------------------------------
# Per-sample analysis
# ---------------------------------------------------------------------------

def analyze_sample(sample_id, guide2gene, ntc_guides):
    """
    Compute log2FC and p-value for each target gene in one sample.
    Returns list of dicts with keys: gene, log2fc, pval, n_perturbed, n_ntc, mean_perturbed, mean_ntc
    """
    h5ad_path = PROD_ROOT / "samples" / sample_id / DS_SUFFIX / "filtered_counts.h5ad"

    with h5py.File(h5ad_path, "r") as f:
        # --- Gene symbols and index (categorical encoding) ---
        gs_group = f["var"]["gene_symbols"]
        gs_cats = np.array([
            x.decode() if isinstance(x, bytes) else x
            for x in gs_group["categories"][:]
        ])
        gs_codes = gs_group["codes"][:]
        gene_symbols = gs_cats[gs_codes]
        symbol_to_idx = {s: i for i, s in enumerate(gene_symbols)}

        # --- Feature call (guide assignment per cell, categorical) ---
        fc_codes = f["obs"]["feature_call"]["codes"][:]
        fc_cats = np.array([
            x.decode() if isinstance(x, bytes) else x
            for x in f["obs"]["feature_call"]["categories"][:]
        ])
        feature_calls = fc_cats[fc_codes]

        # --- Count matrix (CSR) ---
        X_group = f["X"]
        data = X_group["data"][:]
        indices = X_group["indices"][:]
        indptr = X_group["indptr"][:]
        shape = (len(indptr) - 1, len(gene_symbols))
        X = sparse.csr_matrix((data, indices, indptr), shape=shape)

    n_cells = X.shape[0]

    # --- CP10K + log1p normalization ---
    lib_sizes = np.array(X.sum(axis=1)).ravel()
    # Avoid division by zero
    lib_sizes[lib_sizes == 0] = 1.0
    scale_factors = 1e4 / lib_sizes

    # --- Identify NTC cells ---
    ntc_mask = np.array([fc in ntc_guides for fc in feature_calls])
    n_ntc = ntc_mask.sum()
    if n_ntc < MIN_CELLS:
        print(f"  WARNING: {sample_id} has only {n_ntc} NTC cells, skipping")
        return []

    # --- Map guides to target genes, aggregate ---
    # Build gene -> cell indices mapping
    gene_to_cells = {}
    for cell_idx, fc in enumerate(feature_calls):
        tg = guide2gene.get(fc)
        if tg and tg != "Non-Targeting":
            gene_to_cells.setdefault(tg, []).append(cell_idx)

    # --- Test each target gene ---
    results = []
    for gene, cell_idxs in sorted(gene_to_cells.items()):
        if len(cell_idxs) < MIN_CELLS:
            continue
        if gene not in symbol_to_idx:
            continue

        gene_idx = symbol_to_idx[gene]
        cell_idxs = np.array(cell_idxs)

        # Extract expression for this gene, normalize
        perturbed_raw = np.array(X[cell_idxs, gene_idx].todense()).ravel()
        perturbed_norm = np.log1p(perturbed_raw * scale_factors[cell_idxs])

        ntc_raw = np.array(X[ntc_mask, gene_idx].todense()).ravel()
        ntc_norm = np.log1p(ntc_raw * scale_factors[ntc_mask])

        mean_p = perturbed_norm.mean()
        mean_n = ntc_norm.mean()

        # log2 fold change (in log1p-CP10K space)
        # Use pseudocount to avoid log2(0)
        pseudo = 0.01
        log2fc = np.log2(mean_p + pseudo) - np.log2(mean_n + pseudo)

        # Wilcoxon rank-sum (Mann-Whitney U)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                stat, pval = mannwhitneyu(perturbed_norm, ntc_norm, alternative="two-sided")
            except ValueError:
                pval = 1.0

        results.append({
            "gene": gene,
            "log2fc": log2fc,
            "pval": pval,
            "n_perturbed": len(cell_idxs),
            "n_ntc": n_ntc,
            "mean_perturbed": mean_p,
            "mean_ntc": mean_n,
        })

    # BH correction
    if results:
        pvals = np.array([r["pval"] for r in results])
        padj = bh_adjust(pvals)
        for r, p in zip(results, padj):
            r["padj"] = p
            r["significant"] = (p < FDR_THRESH) and (r["log2fc"] < LOG2FC_THRESH)

    return results


# ---------------------------------------------------------------------------
# Waterfall plot (single sample)
# ---------------------------------------------------------------------------

def plot_waterfall(results, sample_id, out_path):
    """Bar chart of log2FC per target gene, sorted by effect size."""
    if not results:
        return

    # Sort by log2FC
    results_sorted = sorted(results, key=lambda r: r["log2fc"])
    genes = [r["gene"] for r in results_sorted]
    log2fcs = [r["log2fc"] for r in results_sorted]
    sig = [r.get("significant", False) for r in results_sorted]

    colors = ["#c0392b" if s else "#bdc3c7" for s in sig]

    fig, ax = plt.subplots(figsize=(14, max(8, len(genes) * 0.18)))
    y_pos = np.arange(len(genes))
    ax.barh(y_pos, log2fcs, color=colors, edgecolor="none", height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_xlabel("log2 Fold Change (vs NTC)", fontsize=11)
    ax.set_title(f"CRISPRi knockdown — {sample_id}", fontsize=13)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.axvline(LOG2FC_THRESH, color="#c0392b", linewidth=0.8, linestyle="--", alpha=0.5)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#c0392b", label=f"Significant (BH p<{FDR_THRESH}, log2FC<{LOG2FC_THRESH})"),
        Patch(facecolor="#bdc3c7", label="Not significant"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9)

    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


# ---------------------------------------------------------------------------
# Cross-sample heatmap
# ---------------------------------------------------------------------------

def plot_heatmap(all_results, out_path):
    """
    Heatmap: rows = target genes (sorted by median log2FC), columns = samples.
    Color = log2FC, significance markers overlaid.
    """
    # Collect all genes tested in at least 1 sample
    all_genes = set()
    for sample_id, results in all_results.items():
        for r in results:
            all_genes.add(r["gene"])

    # Build matrix: genes x samples
    genes_list = sorted(all_genes)
    gene_to_row = {g: i for i, g in enumerate(genes_list)}
    n_genes = len(genes_list)
    n_samples = len(SAMPLES)

    log2fc_mat = np.full((n_genes, n_samples), np.nan)
    sig_mat = np.zeros((n_genes, n_samples), dtype=bool)

    for j, sample_id in enumerate(SAMPLES):
        for r in all_results.get(sample_id, []):
            i = gene_to_row[r["gene"]]
            log2fc_mat[i, j] = r["log2fc"]
            sig_mat[i, j] = r.get("significant", False)

    # Sort genes by median log2FC (most knocked down at top)
    median_fc = np.nanmedian(log2fc_mat, axis=1)
    sort_idx = np.argsort(median_fc)
    log2fc_mat = log2fc_mat[sort_idx]
    sig_mat = sig_mat[sort_idx]
    genes_sorted = [genes_list[i] for i in sort_idx]

    # Clamp for color scale
    vmin, vmax = -3.0, 1.0

    # Color map: blue (knockdown) -> white (no change) -> red (upregulation)
    cmap = plt.cm.RdBu_r  # red = up, blue = down... we want blue=down
    # Actually RdBu has red=high blue=low. We want knockdown (negative) to be blue.
    # RdBu_r: blue=high, red=low. So negative values would be red.
    # Let's use RdBu (not reversed): red=high (positive), blue=low (negative) — knockdown is blue. Good.
    cmap = plt.cm.RdBu

    fig, ax = plt.subplots(figsize=(max(10, n_samples * 0.7 + 3), max(10, n_genes * 0.16 + 2)))

    im = ax.imshow(
        log2fc_mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
        interpolation="nearest"
    )

    # Mark significant cells with a dot
    for i in range(n_genes):
        for j in range(n_samples):
            if sig_mat[i, j]:
                ax.plot(j, i, marker="o", color="black", markersize=2.5, markeredgewidth=0)

    # Labels
    ax.set_xticks(np.arange(n_samples))
    ax.set_xticklabels(SAMPLES, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(n_genes))
    ax.set_yticklabels(genes_sorted, fontsize=5.5)

    # Separator between EBs and iPSC
    ax.axvline(9.5, color="black", linewidth=1.5, linestyle="-")

    # Group labels
    ax.text(4.5, -2, "EBs", ha="center", fontsize=10, fontweight="bold")
    ax.text(12.5, -2, "iPSC", ha="center", fontsize=10, fontweight="bold")

    ax.set_title("CRISPRi target gene knockdown across samples\n(log2FC vs NTC, dot = BH p<0.05 & log2FC<-0.5)", fontsize=12)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("log2 Fold Change", fontsize=10)

    # Mark NaN cells (not tested)
    nan_mask = np.isnan(log2fc_mat)
    for i in range(n_genes):
        for j in range(n_samples):
            if nan_mask[i, j]:
                ax.plot(j, i, marker="x", color="gray", markersize=3, markeredgewidth=0.5)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    guide2gene = load_guide_to_gene(FEAT_REF)
    ntc_guides = get_ntc_guides(guide2gene)
    target_genes = sorted(set(v for v in guide2gene.values() if v != "Non-Targeting"))
    print(f"Loaded {len(guide2gene)} guides -> {len(target_genes)} target genes, {len(ntc_guides)} NTC guides")

    all_results = {}
    for sample_id in SAMPLES:
        print(f"Analyzing {sample_id}...")
        results = analyze_sample(sample_id, guide2gene, ntc_guides)
        all_results[sample_id] = results
        n_sig = sum(1 for r in results if r.get("significant"))
        print(f"  {len(results)} genes tested, {n_sig} significant knockdowns")

    # --- Save full results table ---
    tsv_path = OUT_DIR / "perturbseq_qc.tsv"
    with open(tsv_path, "w") as f:
        header = ["sample", "gene", "log2fc", "pval", "padj", "significant",
                  "n_perturbed", "n_ntc", "mean_perturbed", "mean_ntc"]
        f.write("\t".join(header) + "\n")
        for sample_id in SAMPLES:
            for r in sorted(all_results[sample_id], key=lambda x: x["log2fc"]):
                f.write("\t".join([
                    sample_id, r["gene"],
                    f"{r['log2fc']:.4f}", f"{r['pval']:.2e}", f"{r['padj']:.2e}",
                    str(r["significant"]),
                    str(r["n_perturbed"]), str(r["n_ntc"]),
                    f"{r['mean_perturbed']:.4f}", f"{r['mean_ntc']:.4f}",
                ]) + "\n")
    print(f"\nSaved {tsv_path}")

    # --- Waterfall for EBs1_1 ---
    plot_waterfall(all_results["EBs1_1"], "EBs1_1", OUT_DIR / "waterfall_EBs1_1.png")

    # --- Cross-sample heatmap ---
    plot_heatmap(all_results, OUT_DIR / "knockdown_heatmap.png")

    # --- Summary ---
    print("\n=== Summary ===")
    for sample_id in SAMPLES:
        n_sig = sum(1 for r in all_results[sample_id] if r.get("significant"))
        n_tested = len(all_results[sample_id])
        print(f"  {sample_id:10s}  {n_sig:3d}/{n_tested:3d} significant knockdowns")


if __name__ == "__main__":
    main()
