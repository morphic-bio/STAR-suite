#!/usr/bin/env python3
"""
Perturb-seq QC dotplot: log2FC heatmap with dot size = cell count.

Shows only genes with sufficient cells (median >= MIN_MEDIAN_CELLS across
samples), giving a power-aware view of knockdown effects.

Also produces a cell-count distribution histogram and a power-ranked
summary table.
"""

import csv
import warnings
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from scipy import sparse
from scipy.stats import mannwhitneyu

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

MIN_CELLS = 5
FDR_THRESH = 0.05
LOG2FC_THRESH = -0.5
MIN_MEDIAN_CELLS = 20  # for the focused dotplot
EXPR_THRESH = 1.0      # CP10K threshold for "expressed"

EBS_SAMPLES = [s for s in SAMPLES if s.startswith("EBs")]
IPSC_SAMPLES = [s for s in SAMPLES if s.startswith("iPSC")]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_guide_to_gene(feat_ref_path):
    g2g = {}
    with open(feat_ref_path) as f:
        for row in csv.DictReader(f):
            g2g[row["name"]] = row["target_gene_name"]
    return g2g


def get_ntc_guides(guide2gene):
    return {g for g, t in guide2gene.items() if t == "Non-Targeting"}


def bh_adjust(pvals):
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    adjusted = pvals * n / ranked
    adjusted = np.minimum.accumulate(adjusted[np.argsort(-ranked)])[np.argsort(np.argsort(-ranked))]
    return np.minimum(adjusted, 1.0)


def analyze_sample(sample_id, guide2gene, ntc_guides):
    """Returns list of dicts: gene, log2fc, pval, padj, significant, n_perturbed, n_ntc."""
    h5ad_path = PROD_ROOT / "samples" / sample_id / DS_SUFFIX / "filtered_counts.h5ad"

    with h5py.File(h5ad_path, "r") as f:
        gs_cats = np.array([x.decode() if isinstance(x, bytes) else x for x in f["var"]["gene_symbols"]["categories"][:]])
        gs_codes = f["var"]["gene_symbols"]["codes"][:]
        gene_symbols = gs_cats[gs_codes]
        symbol_to_idx = {s: i for i, s in enumerate(gene_symbols)}

        fc_codes = f["obs"]["feature_call"]["codes"][:]
        fc_cats = np.array([x.decode() if isinstance(x, bytes) else x for x in f["obs"]["feature_call"]["categories"][:]])
        feature_calls = fc_cats[fc_codes]

        X_group = f["X"]
        data = X_group["data"][:]
        indices = X_group["indices"][:]
        indptr = X_group["indptr"][:]
        shape = (len(indptr) - 1, len(gene_symbols))
        X = sparse.csr_matrix((data, indices, indptr), shape=shape)

    lib_sizes = np.array(X.sum(axis=1)).ravel()
    lib_sizes[lib_sizes == 0] = 1.0
    scale_factors = 1e4 / lib_sizes

    ntc_mask = np.array([fc in ntc_guides for fc in feature_calls])
    n_ntc = ntc_mask.sum()
    if n_ntc < MIN_CELLS:
        return []

    gene_to_cells = {}
    for cell_idx, fc in enumerate(feature_calls):
        tg = guide2gene.get(fc)
        if tg and tg != "Non-Targeting":
            gene_to_cells.setdefault(tg, []).append(cell_idx)

    results = []
    for gene, cell_idxs in sorted(gene_to_cells.items()):
        n_pert = len(cell_idxs)
        if gene not in symbol_to_idx:
            # Gene not in expression matrix — still record cell count
            results.append({
                "gene": gene, "log2fc": np.nan, "pval": np.nan,
                "n_perturbed": n_pert, "n_ntc": n_ntc,
            })
            continue

        gene_idx = symbol_to_idx[gene]
        cell_idxs = np.array(cell_idxs)

        perturbed_raw = np.array(X[cell_idxs, gene_idx].todense()).ravel()
        perturbed_norm = np.log1p(perturbed_raw * scale_factors[cell_idxs])
        ntc_raw = np.array(X[ntc_mask, gene_idx].todense()).ravel()
        ntc_norm = np.log1p(ntc_raw * scale_factors[ntc_mask])

        mean_p = perturbed_norm.mean()
        mean_n = ntc_norm.mean()
        ntc_mean_cp10k = (ntc_raw * scale_factors[ntc_mask]).mean()
        pseudo = 0.01
        log2fc = np.log2(mean_p + pseudo) - np.log2(mean_n + pseudo)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                _, pval = mannwhitneyu(perturbed_norm, ntc_norm, alternative="two-sided")
            except ValueError:
                pval = 1.0

        results.append({
            "gene": gene, "log2fc": log2fc, "pval": pval,
            "n_perturbed": n_pert, "n_ntc": n_ntc,
            "mean_perturbed": mean_p, "mean_ntc": mean_n,
            "ntc_mean_cp10k": ntc_mean_cp10k,
        })

    # BH correction
    testable = [r for r in results if not np.isnan(r.get("pval", np.nan))]
    if testable:
        pvals = np.array([r["pval"] for r in testable])
        padj = bh_adjust(pvals)
        for r, p in zip(testable, padj):
            r["padj"] = p
            r["significant"] = (p < FDR_THRESH) and (r["log2fc"] < LOG2FC_THRESH)

    for r in results:
        r.setdefault("padj", np.nan)
        r.setdefault("significant", False)

    return results


# ---------------------------------------------------------------------------
# Cell count histogram
# ---------------------------------------------------------------------------

def plot_cell_count_distribution(all_results, out_path):
    """Histogram of median cell count per target gene across samples."""
    # Collect per-gene median cell counts
    gene_counts = {}  # gene -> list of counts across samples
    for sample_id, results in all_results.items():
        for r in results:
            gene_counts.setdefault(r["gene"], []).append(r["n_perturbed"])

    medians = {g: np.median(cs) for g, cs in gene_counts.items()}
    vals = sorted(medians.values())

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram
    ax = axes[0]
    counts_hist, bin_edges, patches = ax.hist(vals, bins=40, color="#3498db", edgecolor="white", linewidth=0.5)
    ax.axvline(MIN_MEDIAN_CELLS, color="#c0392b", linewidth=1.5, linestyle="--",
               label=f"Threshold ({MIN_MEDIAN_CELLS} cells)")
    ax.set_xlabel("Median cells per target gene (across 16 samples)", fontsize=10)
    ax.set_ylabel("Number of target genes", fontsize=10)
    ax.set_title("Cell count distribution", fontsize=12)
    ax.legend(fontsize=9)
    n_above = sum(1 for v in vals if v >= MIN_MEDIAN_CELLS)
    ax.text(0.95, 0.95, f"{n_above}/{len(vals)} genes\nabove threshold",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8))

    # Label the top 4 singleton genes on their bars
    top4 = sorted(medians.items(), key=lambda x: -x[1])[:4]
    for gene, med in top4:
        # Find which bin this gene falls in
        bin_idx = np.searchsorted(bin_edges[1:], med, side="left")
        bin_center = (bin_edges[bin_idx] + bin_edges[min(bin_idx + 1, len(bin_edges) - 1)]) / 2
        bar_height = counts_hist[min(bin_idx, len(counts_hist) - 1)]
        ax.annotate(
            gene, xy=(bin_center, bar_height), xytext=(0, 5),
            textcoords="offset points", ha="center", va="bottom",
            fontsize=8, fontweight="bold", color="#2c3e50",
        )

    # Cumulative
    ax = axes[1]
    ax.plot(vals, np.arange(1, len(vals) + 1) / len(vals) * 100,
            color="#2c3e50", linewidth=2)
    ax.axvline(MIN_MEDIAN_CELLS, color="#c0392b", linewidth=1.5, linestyle="--")
    ax.set_xlabel("Median cells per target gene", fontsize=10)
    ax.set_ylabel("Cumulative % of genes", fontsize=10)
    ax.set_title("Cumulative distribution", fontsize=12)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


# ---------------------------------------------------------------------------
# Dotplot: focused on well-powered genes
# ---------------------------------------------------------------------------

def plot_dotplot(all_results, out_path):
    """
    Dotplot: rows = target genes (median cells >= threshold), cols = samples.
    Dot color = log2FC, dot size = cell count. Significance ring overlay.
    """
    # Collect per-gene data
    gene_data = {}  # gene -> {sample -> result}
    gene_median_cells = {}
    for sample_id, results in all_results.items():
        for r in results:
            gene_data.setdefault(r["gene"], {})[sample_id] = r

    for gene, sample_results in gene_data.items():
        counts = [sr["n_perturbed"] for sr in sample_results.values()]
        gene_median_cells[gene] = np.median(counts)

    # Filter to well-powered genes
    powered_genes = [g for g, mc in gene_median_cells.items() if mc >= MIN_MEDIAN_CELLS]

    # Sort by median log2FC (strongest knockdown first)
    def median_fc(g):
        fcs = [gene_data[g][s]["log2fc"] for s in SAMPLES if s in gene_data[g] and not np.isnan(gene_data[g][s].get("log2fc", np.nan))]
        return np.median(fcs) if fcs else 0
    powered_genes.sort(key=median_fc)

    n_genes = len(powered_genes)
    n_samples = len(SAMPLES)

    print(f"  {n_genes} genes with median >= {MIN_MEDIAN_CELLS} cells")

    # Build matrices
    fc_mat = np.full((n_genes, n_samples), np.nan)
    count_mat = np.zeros((n_genes, n_samples))
    sig_mat = np.zeros((n_genes, n_samples), dtype=bool)

    for i, gene in enumerate(powered_genes):
        for j, sample_id in enumerate(SAMPLES):
            if sample_id in gene_data[gene]:
                r = gene_data[gene][sample_id]
                fc_mat[i, j] = r.get("log2fc", np.nan)
                count_mat[i, j] = r["n_perturbed"]
                sig_mat[i, j] = r.get("significant", False)

    # Normalize dot sizes: map cell count to area
    max_count = count_mat.max()
    min_dot = 8
    max_dot = 220
    size_mat = min_dot + (count_mat / max(max_count, 1)) * (max_dot - min_dot)

    # Color: diverging blue-white-red
    vmin, vmax = -3.0, 1.0
    cmap = plt.cm.RdBu
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(max(12, n_samples * 0.8 + 4), max(12, n_genes * 0.22 + 2)))

    # Background grid
    for i in range(n_genes):
        ax.axhline(i, color="#f0f0f0", linewidth=0.5, zorder=0)
    for j in range(n_samples):
        ax.axvline(j, color="#f0f0f0", linewidth=0.5, zorder=0)

    # Plot dots
    for i in range(n_genes):
        for j in range(n_samples):
            if np.isnan(fc_mat[i, j]):
                ax.plot(j, i, marker="x", color="#cccccc", markersize=4, zorder=1)
                continue
            color = cmap(norm(np.clip(fc_mat[i, j], vmin, vmax)))
            ax.scatter(j, i, s=size_mat[i, j], c=[color], edgecolors="black" if sig_mat[i, j] else "none",
                       linewidths=1.2 if sig_mat[i, j] else 0, zorder=2)

    # Separator between EBs and iPSC
    ax.axvline(9.5, color="black", linewidth=2, linestyle="-", zorder=3)

    # Labels
    ax.set_xticks(np.arange(n_samples))
    ax.set_xticklabels(SAMPLES, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(n_genes))

    # Gene labels with median cell count annotation
    gene_labels = []
    for gene in powered_genes:
        mc = int(gene_median_cells[gene])
        gene_labels.append(f"{gene} ({mc})")
    ax.set_yticklabels(gene_labels, fontsize=7)

    ax.set_xlim(-0.5, n_samples - 0.5)
    ax.set_ylim(-0.5, n_genes - 0.5)
    ax.invert_yaxis()

    # Group labels
    ax.text(4.5, -1.5, "EBs", ha="center", fontsize=11, fontweight="bold")
    ax.text(12.5, -1.5, "iPSC", ha="center", fontsize=11, fontweight="bold")

    ax.set_title(
        f"CRISPRi knockdown dotplot — {n_genes} genes with median >= {MIN_MEDIAN_CELLS} cells\n"
        f"(dot size = cell count, color = log2FC vs NTC, black ring = significant)",
        fontsize=11
    )

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("log2 Fold Change", fontsize=10)

    # Size legend
    legend_counts = [10, 50, 100, 300, 600]
    legend_sizes = [min_dot + (c / max(max_count, 1)) * (max_dot - min_dot) for c in legend_counts]
    # Place legend below colorbar
    for k, (cnt, sz) in enumerate(zip(legend_counts, legend_sizes)):
        ax.scatter([], [], s=sz, c="gray", edgecolors="none", label=f"{cnt} cells")
    ax.legend(loc="lower right", title="Cell count", fontsize=7, title_fontsize=8,
              framealpha=0.9, scatterpoints=1)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


# ---------------------------------------------------------------------------
# Power-ranked summary table
# ---------------------------------------------------------------------------

def write_power_summary(all_results, out_path):
    """Write a summary table ranked by median cell count and knockdown consistency."""
    gene_data = {}
    for sample_id, results in all_results.items():
        for r in results:
            gene_data.setdefault(r["gene"], {})[sample_id] = r

    rows = []
    for gene, sample_results in gene_data.items():
        counts = [sr["n_perturbed"] for sr in sample_results.values()]
        fcs = [sr["log2fc"] for sr in sample_results.values() if not np.isnan(sr.get("log2fc", np.nan))]
        n_sig = sum(1 for sr in sample_results.values() if sr.get("significant", False))
        n_tested = sum(1 for sr in sample_results.values() if not np.isnan(sr.get("log2fc", np.nan)))

        rows.append({
            "gene": gene,
            "median_cells": int(np.median(counts)),
            "min_cells": min(counts),
            "max_cells": max(counts),
            "median_log2fc": round(np.median(fcs), 3) if fcs else None,
            "n_sig": n_sig,
            "n_tested": n_tested,
            "pct_sig": round(100 * n_sig / n_tested, 1) if n_tested > 0 else 0,
        })

    rows.sort(key=lambda r: (-r["median_cells"],))

    with open(out_path, "w") as f:
        header = ["gene", "median_cells", "min_cells", "max_cells",
                  "median_log2fc", "n_sig_of_16", "pct_sig"]
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join([
                r["gene"],
                str(r["median_cells"]),
                str(r["min_cells"]),
                str(r["max_cells"]),
                str(r["median_log2fc"]) if r["median_log2fc"] is not None else "NA",
                f"{r['n_sig']}/{r['n_tested']}",
                f"{r['pct_sig']}%",
            ]) + "\n")

    print(f"  Saved {out_path}")

    # Print top 30
    print("\n  Top 30 genes by cell count (with knockdown stats):")
    print(f"  {'Gene':20s} {'Med cells':>10s} {'Med log2FC':>10s} {'Sig/Tested':>10s}")
    print(f"  {'-'*20} {'-'*10} {'-'*10} {'-'*10}")
    for r in rows[:30]:
        fc_str = f"{r['median_log2fc']:.3f}" if r['median_log2fc'] is not None else "NA"
        print(f"  {r['gene']:20s} {r['median_cells']:10d} {fc_str:>10s} {r['n_sig']:>3d}/{r['n_tested']:<3d}")


# ---------------------------------------------------------------------------
# Expression-filtered dotplot
# ---------------------------------------------------------------------------

def plot_expressed_dotplot(all_results, out_path):
    """
    Dotplot restricted to genes with NTC baseline expression >= EXPR_THRESH
    in the relevant condition (EBs or iPSC). Condition-aware: a gene is shown
    for EBs columns if expressed in any EBs sample, and for iPSC if expressed
    in any iPSC sample.
    """
    # Collect per-gene per-sample data
    gene_data = {}
    for sample_id, results in all_results.items():
        for r in results:
            gene_data.setdefault(r["gene"], {})[sample_id] = r

    # Determine expression per gene per condition
    gene_ebs_expr = {}   # gene -> median NTC CP10K across EBs
    gene_ipsc_expr = {}  # gene -> median NTC CP10K across iPSC
    for gene, sample_results in gene_data.items():
        ebs_vals = [sr.get("ntc_mean_cp10k", 0) for s, sr in sample_results.items() if s in EBS_SAMPLES]
        ipsc_vals = [sr.get("ntc_mean_cp10k", 0) for s, sr in sample_results.items() if s in IPSC_SAMPLES]
        gene_ebs_expr[gene] = np.median(ebs_vals) if ebs_vals else 0
        gene_ipsc_expr[gene] = np.median(ipsc_vals) if ipsc_vals else 0

    # Genes expressed in at least one condition
    expressed_genes = sorted(set(
        g for g in gene_data
        if gene_ebs_expr.get(g, 0) >= EXPR_THRESH or gene_ipsc_expr.get(g, 0) >= EXPR_THRESH
    ))

    # Sort by median log2FC among expressed conditions
    def sort_key(g):
        fcs = []
        for s, sr in gene_data[g].items():
            if np.isnan(sr.get("log2fc", np.nan)):
                continue
            # Only include if expressed in this condition
            if s in EBS_SAMPLES and gene_ebs_expr.get(g, 0) >= EXPR_THRESH:
                fcs.append(sr["log2fc"])
            elif s in IPSC_SAMPLES and gene_ipsc_expr.get(g, 0) >= EXPR_THRESH:
                fcs.append(sr["log2fc"])
        return np.median(fcs) if fcs else 0
    expressed_genes.sort(key=sort_key)

    n_genes = len(expressed_genes)
    n_samples = len(SAMPLES)
    print(f"  {n_genes} genes expressed (NTC CP10K >= {EXPR_THRESH}) in at least one condition")

    # Build matrices
    fc_mat = np.full((n_genes, n_samples), np.nan)
    count_mat = np.zeros((n_genes, n_samples))
    sig_mat = np.zeros((n_genes, n_samples), dtype=bool)
    not_expressed_mat = np.zeros((n_genes, n_samples), dtype=bool)  # gene not expressed in this condition

    for i, gene in enumerate(expressed_genes):
        for j, sample_id in enumerate(SAMPLES):
            if sample_id not in gene_data.get(gene, {}):
                continue
            r = gene_data[gene][sample_id]
            # Check if gene is expressed in this condition
            if sample_id in EBS_SAMPLES and gene_ebs_expr.get(gene, 0) < EXPR_THRESH:
                not_expressed_mat[i, j] = True
                count_mat[i, j] = r["n_perturbed"]
                continue
            if sample_id in IPSC_SAMPLES and gene_ipsc_expr.get(gene, 0) < EXPR_THRESH:
                not_expressed_mat[i, j] = True
                count_mat[i, j] = r["n_perturbed"]
                continue
            fc_mat[i, j] = r.get("log2fc", np.nan)
            count_mat[i, j] = r["n_perturbed"]
            sig_mat[i, j] = r.get("significant", False)

    # Dot sizes
    max_count = count_mat.max()
    min_dot = 15
    max_dot = 350
    size_mat = min_dot + (count_mat / max(max_count, 1)) * (max_dot - min_dot)

    # Color
    vmin, vmax = -3.5, 1.5
    cmap = plt.cm.RdBu
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    fig_h = max(8, n_genes * 0.38 + 3)
    fig_w = max(12, n_samples * 0.8 + 5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Background grid
    for i in range(n_genes):
        ax.axhline(i, color="#f5f5f5", linewidth=0.5, zorder=0)
    for j in range(n_samples):
        ax.axvline(j, color="#f5f5f5", linewidth=0.5, zorder=0)

    # Plot dots
    for i in range(n_genes):
        for j in range(n_samples):
            if not_expressed_mat[i, j]:
                # Gene not expressed in this condition — grey dash
                ax.plot(j, i, marker="_", color="#cccccc", markersize=6, markeredgewidth=1.5, zorder=1)
                continue
            if np.isnan(fc_mat[i, j]):
                ax.plot(j, i, marker="x", color="#cccccc", markersize=5, zorder=1)
                continue
            color = cmap(norm(np.clip(fc_mat[i, j], vmin, vmax)))
            ax.scatter(j, i, s=size_mat[i, j], c=[color],
                       edgecolors="black" if sig_mat[i, j] else "#888888",
                       linewidths=1.5 if sig_mat[i, j] else 0.3, zorder=2)

    # Separator between EBs and iPSC
    ax.axvline(9.5, color="black", linewidth=2, linestyle="-", zorder=3)

    # Labels
    ax.set_xticks(np.arange(n_samples))
    ax.set_xticklabels(SAMPLES, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(np.arange(n_genes))

    # Gene labels: annotate expression category
    gene_labels = []
    for gene in expressed_genes:
        ebs_e = gene_ebs_expr.get(gene, 0) >= EXPR_THRESH
        ipsc_e = gene_ipsc_expr.get(gene, 0) >= EXPR_THRESH
        if ebs_e and ipsc_e:
            tag = ""
        elif ebs_e:
            tag = " [EBs]"
        else:
            tag = " [iPSC]"
        median_cells = int(np.median([
            gene_data[gene][s]["n_perturbed"]
            for s in SAMPLES if s in gene_data.get(gene, {})
        ]))
        gene_labels.append(f"{gene}{tag}  ({median_cells})")
    ax.set_yticklabels(gene_labels, fontsize=8.5, fontfamily="monospace")

    ax.set_xlim(-0.5, n_samples - 0.5)
    ax.set_ylim(-0.5, n_genes - 0.5)
    ax.invert_yaxis()

    # Group labels
    ax.text(4.5, -1.3, "EBs", ha="center", fontsize=12, fontweight="bold")
    ax.text(12.5, -1.3, "iPSC", ha="center", fontsize=12, fontweight="bold")

    ax.set_title(
        f"CRISPRi knockdown — {n_genes} expressed target genes\n"
        f"dot size = cell count, color = log2FC, black ring = significant (BH p<0.05 & log2FC<-0.5)\n"
        f"grey dash = gene not expressed in that condition (NTC CP10K < {EXPR_THRESH})",
        fontsize=10, pad=15
    )

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.4, pad=0.02, aspect=20)
    cbar.set_label("log2 Fold Change vs NTC", fontsize=10)

    # Size legend
    legend_counts = [10, 50, 100, 300, 600]
    legend_sizes = [min_dot + (c / max(max_count, 1)) * (max_dot - min_dot) for c in legend_counts]
    for cnt, sz in zip(legend_counts, legend_sizes):
        ax.scatter([], [], s=sz, c="gray", edgecolors="none", label=f"{cnt} cells")
    ax.legend(loc="lower right", title="Cell count", fontsize=8, title_fontsize=9,
              framealpha=0.9, scatterpoints=1)

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
        print(f"  {len(results)} genes, {n_sig} significant knockdowns")

    # --- Cell count distribution ---
    print("\n=== Cell count distribution ===")
    plot_cell_count_distribution(all_results, OUT_DIR / "cell_count_distribution.png")

    # --- Power-ranked summary ---
    print("\n=== Power-ranked summary ===")
    write_power_summary(all_results, OUT_DIR / "power_ranked_summary.tsv")

    # --- Dotplot (all genes, cell-count filtered) ---
    print("\n=== Dotplot (well-powered genes) ===")
    plot_dotplot(all_results, OUT_DIR / "knockdown_dotplot.png")

    # --- Dotplot (expression-filtered) ---
    print("\n=== Dotplot (expressed genes only) ===")
    plot_expressed_dotplot(all_results, OUT_DIR / "knockdown_expressed_only.png")


if __name__ == "__main__":
    main()
