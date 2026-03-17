#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate QC histogram from unfiltered downstream AnnData output."
    )
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--min-genes", required=True, type=int)
    parser.add_argument("--max-genes", required=True, type=int)
    parser.add_argument("--mt-pct-cutoff", required=True, type=float)
    parser.add_argument("--n-mad", type=float)
    parser.add_argument("--raw-adaptive-max", type=int)
    args = parser.parse_args()

    adata = ad.read_h5ad(Path(args.input_h5ad).resolve())
    if "singlet" not in adata.obs or "singlet_filtered" not in adata.obs:
        raise KeyError("input h5ad must contain singlet and singlet_filtered obs columns")

    singlet_data = adata[adata.obs["singlet"].astype(bool)].copy()
    filtered_data = adata[adata.obs["singlet_filtered"].astype(bool)].copy()

    if singlet_data.n_obs == 0:
        raise ValueError("Cannot generate QC plot with zero singlet cells")

    gene_counts = singlet_data.obs["n_genes"].to_numpy()
    mt_pcts = singlet_data.obs["mt_pct"].to_numpy()
    filtered_gene_counts = filtered_data.obs["n_genes"].to_numpy() if filtered_data.n_obs else np.array([])

    max_gene_count = max(float(gene_counts.max()), float(args.max_genes), float(args.min_genes), 1.0)
    bin_counts, bin_edges = np.histogram(gene_counts, bins=20, range=(0, max_gene_count))
    max_bin_height = max(int(bin_counts.max()), 1)

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(
        go.Histogram(
            x=gene_counts,
            xbins=dict(start=0, end=max_gene_count, size=max_gene_count / 20),
            marker=dict(color="blue", line=dict(color="white", width=1)),
            opacity=0.7,
            name="All singlets",
        ),
        secondary_y=False,
    )
    fig.add_trace(
        go.Histogram(
            x=filtered_gene_counts,
            xbins=dict(start=0, end=max_gene_count, size=max_gene_count / 20),
            marker=dict(color="green", line=dict(color="white", width=1)),
            opacity=0.7,
            name="Filtered singlets",
        ),
        secondary_y=False,
    )

    mt_per_bin = []
    bin_centers = []
    for i in range(len(bin_edges) - 1):
        mask = (gene_counts >= bin_edges[i]) & (gene_counts < bin_edges[i + 1])
        mt_per_bin.append(float(mt_pcts[mask].mean()) if np.any(mask) else 0.0)
        bin_centers.append((bin_edges[i] + bin_edges[i + 1]) / 2)

    fig.add_trace(
        go.Scatter(
            x=bin_centers,
            y=mt_per_bin,
            mode="lines+markers",
            marker=dict(color="red", size=8),
            line=dict(color="red", width=2),
            name="MT % (mean)",
        ),
        secondary_y=True,
    )
    fig.add_trace(
        go.Scatter(
            x=[args.min_genes, args.min_genes],
            y=[0, max_bin_height * 1.1],
            mode="lines",
            line=dict(color="red", width=2, dash="dash"),
            name=f"Min: {args.min_genes}",
        ),
        secondary_y=False,
    )
    fig.add_trace(
        go.Scatter(
            x=[args.max_genes, args.max_genes],
            y=[0, max_bin_height * 1.1],
            mode="lines",
            line=dict(color="orange", width=2, dash="dash"),
            name=f"Max: {args.max_genes}",
        ),
        secondary_y=False,
    )
    fig.add_trace(
        go.Scatter(
            x=[0, max_gene_count],
            y=[args.mt_pct_cutoff, args.mt_pct_cutoff],
            mode="lines",
            line=dict(color="purple", width=2, dash="dash"),
            name=f"MT% cutoff: {args.mt_pct_cutoff}%",
        ),
        secondary_y=True,
    )

    title = "Gene Distribution and MT Percentage by Quantile"
    if args.n_mad is not None:
        title += f" (Max genes: {args.n_mad} MADs)"
        if args.raw_adaptive_max is not None and args.raw_adaptive_max != args.max_genes:
            title += f", raw={args.raw_adaptive_max}, applied={args.max_genes}"

    fig.update_layout(
        height=600,
        width=900,
        title=dict(text=title, x=0.5, xanchor="center"),
        xaxis_title="Number of genes",
        margin=dict(l=50, r=50, t=80, b=50),
        showlegend=True,
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        barmode="overlay",
    )
    fig.update_yaxes(title_text="Number of cells", range=[0, max_bin_height * 1.1], secondary_y=False)
    fig.update_yaxes(
        title_text="Mitochondrial percentage (%)",
        range=[0, max(float(mt_pcts.max()) * 1.1, float(args.mt_pct_cutoff) * 1.2, 1.0)],
        secondary_y=True,
    )

    filtered_median = float(np.median(filtered_gene_counts)) if filtered_gene_counts.size else float("nan")
    filtered_median_text = "NA" if np.isnan(filtered_median) else str(int(filtered_median))
    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=0.01,
        y=0.99,
        text=f"Min: {args.min_genes} | Max: {args.max_genes} | Passed: {filtered_data.n_obs}/{singlet_data.n_obs}",
        showarrow=False,
        font=dict(size=10),
        align="left",
        bgcolor="white",
        bordercolor="black",
        borderwidth=1,
    )
    fig.add_annotation(
        xref="paper",
        yref="paper",
        x=0.01,
        y=0.93,
        text=f"All: median={int(np.median(gene_counts))} | Filtered: median={filtered_median_text}",
        showarrow=False,
        font=dict(size=10),
        align="left",
        bgcolor="white",
        bordercolor="black",
        borderwidth=1,
    )

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    html_path = output_dir / "gene_quantile_histogram.html"
    png_path = output_dir / "gene_quantile_histogram.png"
    fig.write_html(html_path)
    fig.write_image(png_path, scale=2)
    print(f"Wrote {html_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
