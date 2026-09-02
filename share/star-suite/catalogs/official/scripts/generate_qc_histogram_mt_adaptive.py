#!/usr/bin/env python3
"""Generate the QC histogram for the adaptive mt% guard (two mt rejection lines).

Companion to generate_qc_histogram.py. Used by BOTH the new-implementation and
the conversion paths to (re)draw gene_quantile_histogram.{html,png}.

Difference from the legacy graph: the single 5% MT cutoff line is replaced by
TWO horizontal lines on the mitochondrial-% axis —

  * MT floor                — the strict 5% guard (cells at/below are always kept)
  * MT adaptive (median+k*MAD) — the per-sample MAD soft guard

The effective cut is the higher of the two; it is annotated on the plot. All
threshold values are read from adaptive_qc_threshold.json (the single source of
truth written by compute_adaptive_qc_threshold.py + apply_adaptive_mt_filter.py
/ convert_h5ad_mt_adaptive.py).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def write_mt_histogram(
    mt_pcts,
    output_dir,
    *,
    mt_floor: float,
    mt_raw: float,
    mt_threshold: float,
    n_mad: float,
    high_pct: float,
    high_fraction: float,
    high_fraction_limit: float,
) -> None:
    """Emit mt_quantile_histogram.{html,png} — only for flagged samples.

    A single-panel mt% distribution over singlets with the strict floor and the
    adaptive median+n_mad*MAD cut drawn on it, plus the high-mt review band.
    Generated only when the sample-level mt% flag fires, so a reviewer has a
    graph to judge the adaptive cut against the actual distribution.
    """
    mt = np.asarray(mt_pcts, dtype=float)
    mt = mt[np.isfinite(mt)]
    upper = max(100.0, float(mt.max()) if mt.size else 100.0)
    counts, _ = np.histogram(mt, bins=max(int(upper), 1), range=(0, upper))
    top = max(int(counts.max()), 1)

    fig = go.Figure()
    fig.add_shape(
        type="rect", x0=high_pct, x1=upper, y0=0, y1=top * 1.1,
        fillcolor="#A43B4A", opacity=0.08, line_width=0, layer="below",
    )
    fig.add_trace(
        go.Histogram(
            x=mt,
            xbins=dict(start=0, end=upper, size=1.0),
            marker=dict(color="#2E86AB", line=dict(color="white", width=0.5)),
            opacity=0.85,
            name="Singlets",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[mt_floor, mt_floor], y=[0, top * 1.1], mode="lines",
            line=dict(color="purple", width=2, dash="dash"),
            name=f"MT floor: {mt_floor:g}%",
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[mt_raw, mt_raw], y=[0, top * 1.1], mode="lines",
            line=dict(color="darkred", width=2, dash="dot"),
            name=f"MT adaptive (median+{n_mad:g}·MAD): {mt_raw:.2f}%",
        )
    )
    fig.update_layout(
        height=600,
        width=900,
        title=dict(
            text="Mitochondrial % Distribution (singlets) — FLAGGED for review",
            x=0.5, xanchor="center",
        ),
        xaxis_title="mt% per cell",
        yaxis_title="Number of singlet cells",
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(l=50, r=50, t=90, b=50),
    )
    fig.update_xaxes(range=[0, upper])
    fig.update_yaxes(range=[0, top * 1.1])
    fig.add_annotation(
        xref="paper", yref="paper", x=0.99, y=0.99,
        text=(
            f"FLAGGED: {high_fraction:.1%} of singlets above {high_pct:g}% mt "
            f"(limit {high_fraction_limit:.0%}) | effective cut {mt_threshold:.2f}%"
        ),
        showarrow=False, font=dict(size=10, color="#A43B4A"), align="right",
        bgcolor="white", bordercolor="#A43B4A", borderwidth=1,
    )

    html_path = output_dir / "mt_quantile_histogram.html"
    png_path = output_dir / "mt_quantile_histogram.png"
    fig.write_html(html_path)
    fig.write_image(png_path, scale=2)
    print(f"Wrote {html_path}")
    print(f"Wrote {png_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--threshold-json", required=True,
        help="adaptive_qc_threshold.json with n_genes bounds + mt_pct_* keys.",
    )
    args = parser.parse_args()

    with open(Path(args.threshold_json).resolve(), "r", encoding="utf-8") as handle:
        thr = json.load(handle)
    required = [
        "min_genes", "effective_max_genes",
        "mt_pct_floor", "mt_pct_raw_threshold", "mt_pct_threshold", "mt_pct_n_mad",
    ]
    missing = [k for k in required if k not in thr]
    if missing:
        raise KeyError(
            f"threshold JSON missing keys {missing} — run compute_adaptive_qc_threshold.py "
            f"and apply_adaptive_mt_filter.py / convert_h5ad_mt_adaptive.py first"
        )
    min_genes = int(thr["min_genes"])
    max_genes = int(thr["effective_max_genes"])
    mt_floor = float(thr["mt_pct_floor"])
    mt_raw = float(thr["mt_pct_raw_threshold"])
    mt_threshold = float(thr["mt_pct_threshold"])
    n_mad = float(thr["mt_pct_n_mad"])

    adata = ad.read_h5ad(Path(args.input_h5ad).resolve())
    if "singlet" not in adata.obs or "singlet_filtered" not in adata.obs:
        raise KeyError("input h5ad must contain singlet and singlet_filtered obs columns")

    singlet_data = adata[adata.obs["singlet"].astype(bool)].copy()
    filtered_data = adata[adata.obs["singlet_filtered"].astype(bool)].copy()
    if singlet_data.n_obs == 0:
        raise ValueError("Cannot generate QC plot with zero singlet cells")

    gene_counts = singlet_data.obs["n_genes"].to_numpy()
    mt_pcts = singlet_data.obs["mt_pct"].to_numpy()
    filtered_gene_counts = (
        filtered_data.obs["n_genes"].to_numpy() if filtered_data.n_obs else np.array([])
    )

    max_gene_count = max(
        float(gene_counts.max()), float(max_genes), float(min_genes), 1.0
    )
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
    # n_genes bounds (primary axis).
    fig.add_trace(
        go.Scatter(
            x=[min_genes, min_genes],
            y=[0, max_bin_height * 1.1],
            mode="lines",
            line=dict(color="red", width=2, dash="dash"),
            name=f"Min genes: {min_genes}",
        ),
        secondary_y=False,
    )
    fig.add_trace(
        go.Scatter(
            x=[max_genes, max_genes],
            y=[0, max_bin_height * 1.1],
            mode="lines",
            line=dict(color="orange", width=2, dash="dash"),
            name=f"Max genes: {max_genes}",
        ),
        secondary_y=False,
    )
    # Two mt% rejection lines (secondary axis).
    effective_is_floor = mt_threshold <= mt_floor + 1e-9
    fig.add_trace(
        go.Scatter(
            x=[0, max_gene_count],
            y=[mt_floor, mt_floor],
            mode="lines",
            line=dict(color="purple", width=2, dash="dash"),
            name=(
                f"MT floor: {mt_floor:.1f}%"
                + (" (effective)" if effective_is_floor else "")
            ),
        ),
        secondary_y=True,
    )
    fig.add_trace(
        go.Scatter(
            x=[0, max_gene_count],
            y=[mt_raw, mt_raw],
            mode="lines",
            line=dict(color="darkred", width=2, dash="dot"),
            name=(
                f"MT adaptive (median+{n_mad:g}·MAD): {mt_raw:.2f}%"
                + ("" if effective_is_floor else " (effective)")
            ),
        ),
        secondary_y=True,
    )

    title = "Gene Distribution and MT Percentage by Quantile"
    fig.update_layout(
        height=600,
        width=900,
        title=dict(text=title, x=0.5, xanchor="center", y=0.97),
        xaxis_title="Number of genes",
        margin=dict(l=50, r=50, t=110, b=50),
        showlegend=True,
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        barmode="overlay",
    )
    fig.update_yaxes(
        title_text="Number of cells", range=[0, max_bin_height * 1.1], secondary_y=False
    )
    fig.update_yaxes(
        title_text="Mitochondrial percentage (%)",
        range=[0, max(float(mt_pcts.max()) * 1.1, mt_threshold * 1.2, 1.0)],
        secondary_y=True,
    )

    filtered_median = (
        float(np.median(filtered_gene_counts)) if filtered_gene_counts.size else float("nan")
    )
    filtered_median_text = "NA" if np.isnan(filtered_median) else str(int(filtered_median))
    fig.add_annotation(
        xref="paper", yref="paper", x=0.01, y=0.99,
        text=(
            f"Min genes: {min_genes} | Max genes: {max_genes} | "
            f"MT effective cut: {mt_threshold:.2f}% | "
            f"Passed: {filtered_data.n_obs}/{singlet_data.n_obs}"
        ),
        showarrow=False, font=dict(size=10), align="left",
        bgcolor="white", bordercolor="black", borderwidth=1,
    )
    fig.add_annotation(
        xref="paper", yref="paper", x=0.01, y=0.93,
        text=(
            f"All: median={int(np.median(gene_counts))} | "
            f"Filtered: median={filtered_median_text}"
        ),
        showarrow=False, font=dict(size=10), align="left",
        bgcolor="white", bordercolor="black", borderwidth=1,
    )

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    html_path = output_dir / "gene_quantile_histogram.html"
    png_path = output_dir / "gene_quantile_histogram.png"
    fig.write_html(html_path)
    fig.write_image(png_path, scale=2)
    print(f"Wrote {html_path}")
    print(f"Wrote {png_path}")

    # Standalone mt% distribution histogram — emitted ONLY when the sample-level
    # mt% flag fired, so flagged samples have a graph to support human review.
    if bool(thr.get("mt_pct_flag", False)):
        write_mt_histogram(
            mt_pcts,
            output_dir,
            mt_floor=mt_floor,
            mt_raw=mt_raw,
            mt_threshold=mt_threshold,
            n_mad=n_mad,
            high_pct=float(thr.get("mt_pct_flag_high_pct", 20.0)),
            high_fraction=float(thr.get("mt_pct_flag_high_fraction", 0.0)),
            high_fraction_limit=float(thr.get("mt_pct_flag_high_fraction_limit", 0.10)),
        )
    else:
        print("Sample not flagged (mt_pct_flag false/absent) — mt% histogram skipped")


if __name__ == "__main__":
    main()
