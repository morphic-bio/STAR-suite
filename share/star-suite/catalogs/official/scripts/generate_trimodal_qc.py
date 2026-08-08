#!/usr/bin/env python3
"""Render a single unified trimodal QC figure (Plotly HTML+JSON+PNG) over a multiome
MuData with RNA(GEX) + ATAC-peak + (optional) CRISPR-guide modalities.

This is the GO/NO-GO endpoint instrument for the single-pass multiome path: one figure
that *requires* all three modalities and answers "did GEX, ATAC, and guide capture all
work, on the same cells, well enough to proceed?". It is a view over an assembled MuData
(build it with build_multiome_mudata.py, plus a guide modality), so it cannot be produced
without the trimodal object -- which is the point.

The cell basis must be the *finalized* cell call, not a permissive candidate set. For
STAR/STARsolo EmptyDrops_CR, the matrix written to outs/filtered_feature_bc_matrix is the
permissive candidate set (knee + a rescued low-UMI tail); use the EmptyDrops knee
(is_simple_cell==1) or pass it to build_multiome_mudata.py via --filtered-barcodes so the
QC reports real cells. See docs/trimodal-qc.md.

Panels (2x3): GEX UMIs vs genes (colour %mito); ATAC in-peak depth vs #peaks (colour
FRiP); guide assignment status; cells-per-modality overlap; joint GEX UMAP coloured by
guide (perturbation); joint GEX UMAP coloured by ATAC depth/FRiP.

Outputs: <prefix>.html, <prefix>.json, <prefix>.png (PNG via kaleido).
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import mudata as md
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render the unified trimodal (GEX+ATAC+guide) QC figure from a MuData."
    )
    parser.add_argument("--mudata", required=True, help="Input .h5mu with RNA + ATAC (+ guide) modalities")
    parser.add_argument("--output-prefix", required=True, help="Output path prefix; writes <prefix>.{html,json,png}")
    parser.add_argument("--rna-mod", default="rna", help="RNA/GEX modality name (default: rna)")
    parser.add_argument("--atac-mod", default="atac", help="ATAC peak modality name (default: atac)")
    parser.add_argument("--guide-mod", default="guide", help="CRISPR-guide modality name (default: guide; skipped if absent)")
    parser.add_argument("--guide-gene-sep", default="-", help="Separator splitting a guide name into its target gene (default '-', e.g. GENE-1)")
    parser.add_argument("--frip-column", default="atac_frip", help="Per-cell FRiP column in atac.obs, if present")
    parser.add_argument("--atac-metrics", help="Optional per-barcode ATAC metrics TSV (provides FRiP if not in atac.obs)")
    parser.add_argument("--atac-metrics-barcode-column", default="barcode")
    parser.add_argument("--atac-metrics-frip-column", default="atac_peak_fraction")
    parser.add_argument("--strip-barcode-suffix", action="store_true", help="Strip a trailing -1 from metrics barcodes before joining")
    parser.add_argument("--umap-obsm", default="X_umap", help="Reuse this rna.obsm embedding if present, else compute a UMAP")
    parser.add_argument("--title", default="Single-pass trimodal QC (GEX + ATAC + CRISPR guide)")
    parser.add_argument("--no-png", action="store_true", help="Skip the static PNG export (kaleido)")
    return parser.parse_args()


def _strip1(barcodes) -> list[str]:
    return [b[:-2] if str(b).endswith("-1") else str(b) for b in barcodes]


def _colsum(X) -> np.ndarray:
    return np.asarray(X.sum(1)).ravel()


def _nnz_per_row(X) -> np.ndarray:
    return np.asarray((X > 0).sum(1)).ravel()


def compute_qc(mdata: md.MuData, args: argparse.Namespace) -> tuple[pd.DataFrame, Optional[str]]:
    """Per-cell, per-modality QC on the RNA cell index. Returns (qc, frip_column_or_None)."""
    rna = mdata.mod[args.rna_mod]
    qc = pd.DataFrame(index=rna.obs_names)
    qc["gex_umis"] = _colsum(rna.X)
    qc["gex_genes"] = _nnz_per_row(rna.X)
    mito = rna.var_names.str.upper().str.startswith(("MT-", "MT."))
    if mito.any():
        qc["gex_pct_mito"] = 100 * _colsum(rna[:, mito].X) / np.maximum(qc["gex_umis"], 1)
    else:
        qc["gex_pct_mito"] = 0.0

    frip: Optional[str] = None
    if args.atac_mod in mdata.mod:
        atac = mdata.mod[args.atac_mod]
        aidx = pd.Index(atac.obs_names)
        qc["atac_inpeak"] = pd.Series(_colsum(atac.X), index=aidx).reindex(qc.index).values
        qc["atac_npeaks"] = pd.Series(_nnz_per_row(atac.X), index=aidx).reindex(qc.index).values
        if args.frip_column in atac.obs.columns:
            qc["atac_frip"] = pd.Series(atac.obs[args.frip_column].values, index=aidx).reindex(qc.index).values
            frip = "atac_frip"
        elif args.atac_metrics:
            m = pd.read_csv(args.atac_metrics, sep="\t")
            bc = _strip1(m[args.atac_metrics_barcode_column]) if args.strip_barcode_suffix else m[args.atac_metrics_barcode_column].astype(str)
            m = m.assign(_bc=list(bc)).set_index("_bc")
            if args.atac_metrics_frip_column in m.columns and len(set(m.index) & set(qc.index)) > 0.5 * len(qc):
                qc["atac_frip"] = m[args.atac_metrics_frip_column].reindex(qc.index).values
                frip = "atac_frip"

    if args.guide_mod in mdata.mod:
        guide = mdata.mod[args.guide_mod]
        gx = guide.X.toarray() if sp.issparse(guide.X) else np.asarray(guide.X)
        gdf = pd.DataFrame(gx > 0, index=guide.obs_names, columns=guide.var_names).reindex(qc.index)
        considered = gdf.notna().any(axis=1)
        gbool = gdf.fillna(False).astype(bool)
        sep = args.guide_gene_sep
        genes_per_cell = gbool.apply(lambda r: sorted({c.split(sep)[0] for c in gbool.columns[r.values]}), axis=1)
        qc["n_guides"] = gbool.sum(1)
        n_genes = genes_per_cell.apply(len)
        qc["guide_status"] = np.where(~considered, "not_in_guide_set",
                              np.where(n_genes == 0, "unassigned",
                              np.where(n_genes == 1, "singlet", "multiplet")))
        qc["assigned_gene"] = genes_per_cell.apply(lambda s: s[0] if len(s) == 1 else ("multiplet" if len(s) > 1 else "none"))
    else:
        qc["guide_status"] = "not_in_guide_set"
        qc["assigned_gene"] = "none"
    return qc, frip


def ensure_umap(mdata: md.MuData, qc: pd.DataFrame, args: argparse.Namespace) -> None:
    """Reuse rna.obsm[umap_obsm] if present, else compute a quick GEX UMAP."""
    rna = mdata.mod[args.rna_mod]
    if args.umap_obsm in rna.obsm:
        um = pd.DataFrame(rna.obsm[args.umap_obsm][:, :2], index=rna.obs_names, columns=["umap1", "umap2"])
    else:
        gx = rna.copy()
        sc.pp.normalize_total(gx, target_sum=1e4)
        sc.pp.log1p(gx)
        sc.pp.highly_variable_genes(gx, n_top_genes=2000)
        gx = gx[:, gx.var.highly_variable].copy()
        sc.pp.scale(gx, max_value=10)
        sc.tl.pca(gx, n_comps=30)
        sc.pp.neighbors(gx, n_neighbors=15)
        sc.tl.umap(gx)
        um = pd.DataFrame(gx.obsm["X_umap"], index=gx.obs_names, columns=["umap1", "umap2"])
    um = um.reindex(qc.index)
    qc["umap1"], qc["umap2"] = um["umap1"], um["umap2"]


def render_qc(mdata: md.MuData, output_prefix: str, args: argparse.Namespace) -> dict:
    """Compute QC, render the 2x3 dashboard, write {prefix}.{html,json,png}. Returns count summary."""
    qc, frip = compute_qc(mdata, args)
    ensure_umap(mdata, qc, args)

    n_gex = int(mdata.mod[args.rna_mod].n_obs)
    n_atac = int(qc["atac_inpeak"].notna().sum()) if "atac_inpeak" in qc else 0
    assigned_mask = qc["guide_status"].isin(["singlet", "multiplet"])
    n_assigned = int(assigned_mask.sum())
    triple = int((qc.get("atac_inpeak", pd.Series(index=qc.index)).notna() & assigned_mask).sum())
    rate = 100 * n_assigned / max(n_gex, 1)

    # Portrait 3x2 layout: colorbar panels (GEX, ATAC, UMAP-depth) fill the right
    # column so colorbars stack down the right edge; bars + guide UMAP fill the
    # left column; guide legend goes below.
    ROWS, COLS, VS = 3, 2, 0.10
    _rh = (1.0 - (ROWS - 1) * VS) / ROWS
    row_center = {r: 1.0 - (r - 1) * (_rh + VS) - _rh / 2 for r in range(1, ROWS + 1)}

    def cbar(row, title):
        return dict(title=dict(text=title, side="right", font=dict(size=12)),
                    tickfont=dict(size=11), x=1.005, xanchor="left",
                    y=row_center[row], yanchor="middle", len=_rh * 0.9, thickness=12)

    fig = make_subplots(rows=ROWS, cols=COLS, vertical_spacing=VS, horizontal_spacing=0.15,
        subplot_titles=(
            "Cells per modality / overlap", "GEX: UMIs vs genes (color=%mito)",
            "Guide: assignment status", "ATAC: in-peak depth vs #peaks",
            "Joint UMAP — guide (perturbation)", "Joint UMAP — ATAC depth"))
    MS = 4

    # (1,1) cells per modality / overlap
    fig.add_trace(go.Bar(x=["GEX", "ATAC", "guide-assigned", "triple"],
        y=[n_gex, n_atac, n_assigned, triple], showlegend=False,
        marker_color=["#1f77b4", "#9467bd", "#2ca02c", "#000000"]), 1, 1)
    fig.update_yaxes(title="cells", row=1, col=1)

    # (1,2) GEX
    fig.add_trace(go.Scattergl(x=qc.gex_umis, y=qc.gex_genes, mode="markers",
        marker=dict(size=MS, color=qc.gex_pct_mito, colorscale="Viridis", cmax=20, cmin=0,
                    colorbar=cbar(1, "%mito")), showlegend=False), 1, 2)
    fig.update_xaxes(type="log", row=1, col=2, title="UMIs"); fig.update_yaxes(type="log", row=1, col=2, title="genes")

    # (2,1) guide assignment status
    vc = qc["guide_status"].value_counts()
    fig.add_trace(go.Bar(x=vc.index.tolist(), y=vc.values.tolist(), showlegend=False,
        marker_color=["#2ca02c", "#1f77b4", "#ff7f0e", "#999999"][:len(vc)]), 2, 1)

    # (2,2) ATAC
    if "atac_inpeak" in qc:
        fig.add_trace(go.Scattergl(x=qc.atac_inpeak, y=qc.atac_npeaks, mode="markers",
            marker=dict(size=MS, color=(qc[frip] if frip else qc.atac_inpeak), colorscale="Cividis",
                        colorbar=cbar(2, "FRiP" if frip else "depth")),
            showlegend=False), 2, 2)
        fig.update_xaxes(type="log", row=2, col=2, title="in-peak fragments"); fig.update_yaxes(type="log", row=2, col=2, title="# peaks")

    # (3,1) joint UMAP — guide
    sub = qc.dropna(subset=["umap1"])
    topg = sub.loc[sub.guide_status == "singlet", "assigned_gene"].value_counts().head(8).index.tolist()
    for gene in topg:
        s = sub[sub.assigned_gene == gene]
        fig.add_trace(go.Scattergl(x=s.umap1, y=s.umap2, mode="markers", name=gene,
            marker=dict(size=MS), legendgroup="g"), 3, 1)
    bg = sub[~sub.assigned_gene.isin(topg)]
    fig.add_trace(go.Scattergl(x=bg.umap1, y=bg.umap2, mode="markers", name="other/none",
        marker=dict(size=3, color="#dddddd"), legendgroup="g"), 3, 1)

    # (3,2) joint UMAP — ATAC depth
    cvar = qc[frip] if frip else qc.get("atac_inpeak")
    if cvar is not None:
        fig.add_trace(go.Scattergl(x=sub.umap1, y=sub.umap2, mode="markers", showlegend=False,
            marker=dict(size=MS, color=cvar.reindex(sub.index), colorscale="Cividis",
                        colorbar=cbar(3, "FRiP" if frip else "ATAC depth"))), 3, 2)

    fig.update_layout(height=1000, width=820, template="plotly_white",
        font=dict(size=15),
        legend=dict(itemsizing="constant", font=dict(size=12), orientation="h",
                    yanchor="top", y=-0.05, xanchor="center", x=0.5),
        margin=dict(l=78, r=112, t=92, b=80),
        title=dict(text=(f"{args.title}<br>GEX={n_gex}, ATAC={n_atac}, "
                    f"guide-assigned={n_assigned} ({rate:.1f}%), triple={triple}"),
                   font=dict(size=15), x=0.5, xanchor="center", y=0.99, yanchor="top"))
    fig.update_annotations(font_size=13)

    out = Path(output_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(f"{out}.html", include_plotlyjs="cdn")
    Path(f"{out}.json").write_text(fig.to_json())
    if not args.no_png:
        fig.write_image(f"{out}.png", scale=2)
    summary = {"gex": n_gex, "atac": n_atac, "guide_assigned": n_assigned, "triple": triple, "assignment_rate_pct": round(rate, 1)}
    print(f">> QC summary: {summary}")
    print(f">> wrote {out}.html, {out}.json" + ("" if args.no_png else f", {out}.png"))
    return summary


def main() -> None:
    args = parse_args()
    print(f">> reading MuData: {args.mudata}")
    mdata = md.read(args.mudata)
    render_qc(mdata, args.output_prefix, args)


if __name__ == "__main__":
    main()
