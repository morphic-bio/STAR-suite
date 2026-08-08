#!/usr/bin/env python3
"""Render four-factor QC for RNA + ATAC + protein/ADT + identity MuData."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import mudata as md
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import scipy.sparse as sp
from plotly.subplots import make_subplots


DEFAULT_IDENTITY_MODS = ("guide", "state", "hash")
DEFAULT_MARKER_PAIRS = (
    ("ITGB1", "CD29"),
    ("CD46", "CD46"),
    ("MS4A1", "CD20"),
    ("CD3D", "CD3"),
    ("CD4", "CD4"),
    ("CD8A", "CD8"),
    ("LYZ", "CD14"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render protein-aware four-factor QC from a MuData with RNA, ATAC, "
            "protein/ADT, and optional guide/hash/state modalities."
        )
    )
    parser.add_argument("--mudata", required=True, help="Input .h5mu")
    parser.add_argument("--output-prefix", required=True, help="Output path prefix")
    parser.add_argument("--rna-mod", default="rna")
    parser.add_argument("--atac-mod", default="atac")
    parser.add_argument("--protein-mod", default="protein")
    parser.add_argument(
        "--identity-mod",
        action="append",
        default=None,
        help=(
            "Feature-barcode identity modality to summarize. May be repeated. "
            "Defaults to guide, state, and hash when present."
        ),
    )
    parser.add_argument(
        "--min-identity-umis",
        type=int,
        default=1,
        help="Minimum count for a feature to count as an identity assignment",
    )
    parser.add_argument(
        "--guide-gene-sep",
        default="-",
        help="Optional separator used to collapse guide feature names to target labels",
    )
    parser.add_argument("--frip-column", default="atac_frip")
    parser.add_argument("--atac-metrics", help="Optional per-barcode ATAC metrics TSV")
    parser.add_argument("--atac-metrics-barcode-column", default="barcode")
    parser.add_argument("--atac-metrics-frip-column", default="atac_peak_fraction")
    parser.add_argument("--strip-barcode-suffix", action="store_true")
    parser.add_argument("--umap-obsm", default="X_umap")
    parser.add_argument("--protein-marker", help="Protein/ADT marker to show on UMAP")
    parser.add_argument("--rna-marker", help="RNA gene to compare with the selected protein marker")
    parser.add_argument("--title", default="Four-factor QC (RNA + ATAC + protein + identity)")
    parser.add_argument("--no-png", action="store_true", help="Skip static PNG export")
    return parser.parse_args()


def read_table_auto(path: str | Path) -> pd.DataFrame:
    name = str(path).lower()
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def strip_suffix(barcodes: Any) -> list[str]:
    return [value[:-2] if str(value).endswith("-1") else str(value) for value in barcodes]


def row_sum(matrix: Any) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel()


def row_nnz(matrix: Any) -> np.ndarray:
    return np.asarray((matrix > 0).sum(axis=1)).ravel()


def dense(matrix: Any) -> np.ndarray:
    return matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)


def counts_matrix(adata: Any) -> Any:
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def feature_labels(adata: Any) -> list[str]:
    if "feature_names" in adata.var:
        return adata.var["feature_names"].astype(str).tolist()
    if "gene_symbols" in adata.var:
        return adata.var["gene_symbols"].astype(str).tolist()
    return adata.var_names.astype(str).tolist()


def reindex_values(values: Any, source_index: Any, target_index: Any, fill_value: Any = np.nan) -> np.ndarray:
    return pd.Series(values, index=pd.Index(source_index)).reindex(pd.Index(target_index)).fillna(fill_value).values


def obs_bool_or_counts(adata: Any, column: str, target_index: Any) -> np.ndarray:
    if column in adata.obs:
        return reindex_values(adata.obs[column].astype(bool).values, adata.obs_names, target_index, False).astype(bool)
    return reindex_values(row_sum(counts_matrix(adata)) > 0, adata.obs_names, target_index, False).astype(bool)


def aligned_dense_counts(adata: Any, target_index: Any) -> tuple[np.ndarray, np.ndarray]:
    indexer = adata.obs_names.get_indexer(pd.Index(target_index))
    present = indexer >= 0
    out = np.zeros((len(target_index), adata.n_vars), dtype=float)
    if present.any():
        out[present, :] = dense(counts_matrix(adata))[indexer[present], :]
    return out, present


def find_feature_index(adata: Any, marker: str | None) -> int | None:
    if not marker:
        return None
    marker_lower = marker.lower()
    candidates = [adata.var_names.astype(str).tolist()]
    if "feature_names" in adata.var:
        candidates.append(adata.var["feature_names"].astype(str).tolist())
    if "gene_symbols" in adata.var:
        candidates.append(adata.var["gene_symbols"].astype(str).tolist())
    for labels in candidates:
        for idx, label in enumerate(labels):
            if label == marker or label.lower() == marker_lower:
                return idx
    return None


def select_marker_pair(rna: Any, protein: Any, args: argparse.Namespace) -> tuple[str | None, str | None]:
    if args.protein_marker:
        if args.rna_marker:
            return args.rna_marker, args.protein_marker
        mapping = {protein_marker.lower(): rna_marker for rna_marker, protein_marker in DEFAULT_MARKER_PAIRS}
        mapped_rna_marker = mapping.get(args.protein_marker.lower())
        if mapped_rna_marker and find_feature_index(rna, mapped_rna_marker) is not None:
            return mapped_rna_marker, args.protein_marker
        return None, args.protein_marker
    for rna_marker, protein_marker in DEFAULT_MARKER_PAIRS:
        if find_feature_index(rna, rna_marker) is not None and find_feature_index(protein, protein_marker) is not None:
            return rna_marker, protein_marker
    labels = feature_labels(protein)
    totals = row_sum(counts_matrix(protein).T)
    if len(labels) and len(totals):
        return None, labels[int(np.argmax(totals))]
    return None, None


def add_identity_qc(mdata: md.MuData, qc: pd.DataFrame, args: argparse.Namespace) -> list[str]:
    mods = tuple(args.identity_mod) if args.identity_mod else DEFAULT_IDENTITY_MODS
    labels_per_cell: list[list[str]] = [[] for _ in range(qc.shape[0])]
    considered = np.zeros(qc.shape[0], dtype=bool)
    n_called = np.zeros(qc.shape[0], dtype=int)
    used_mods: list[str] = []

    for mod in mods:
        if mod not in mdata.mod:
            continue
        used_mods.append(mod)
        adata = mdata.mod[mod]
        matrix, _ = aligned_dense_counts(adata, qc.index)
        present = obs_bool_or_counts(adata, f"{mod}_barcode_present", qc.index)
        considered |= present
        call_mask = matrix >= args.min_identity_umis
        labels = feature_labels(adata)
        if mod == "guide" and args.guide_gene_sep:
            labels = [label.split(args.guide_gene_sep)[0] for label in labels]
        for row_idx in np.flatnonzero(call_mask.any(axis=1)):
            called = sorted({labels[col_idx] for col_idx in np.flatnonzero(call_mask[row_idx, :])})
            labels_per_cell[row_idx].extend(f"{mod}:{label}" for label in called)
            n_called[row_idx] += len(called)

    status = np.full(qc.shape[0], "not_in_identity_set", dtype=object)
    status[considered] = "unassigned"
    status[n_called == 1] = "singlet"
    status[n_called > 1] = "multiplet"
    display_labels = []
    for labels in labels_per_cell:
        if not labels:
            display_labels.append("none")
        elif len(labels) == 1:
            display_labels.append(labels[0])
        else:
            display_labels.append("multiplet")
    qc["identity_considered"] = considered
    qc["identity_assigned"] = n_called > 0
    qc["identity_n_features"] = n_called
    qc["identity_status"] = status
    qc["identity_label"] = display_labels
    return used_mods


def compute_qc(mdata: md.MuData, args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, Any]]:
    rna = mdata.mod[args.rna_mod]
    qc = pd.DataFrame(index=rna.obs_names)
    qc["gex_umis"] = row_sum(counts_matrix(rna))
    qc["gex_genes"] = row_nnz(counts_matrix(rna))
    mito = pd.Index(rna.var_names.astype(str)).str.upper().str.startswith(("MT-", "MT."))
    qc["gex_pct_mito"] = 100 * row_sum(counts_matrix(rna[:, mito])) / np.maximum(qc["gex_umis"], 1) if mito.any() else 0.0

    frip_column: str | None = None
    qc["atac_present"] = False
    if args.atac_mod in mdata.mod:
        atac = mdata.mod[args.atac_mod]
        qc["atac_inpeak"] = reindex_values(row_sum(counts_matrix(atac)), atac.obs_names, qc.index)
        qc["atac_npeaks"] = reindex_values(row_nnz(counts_matrix(atac)), atac.obs_names, qc.index)
        qc["atac_present"] = qc["atac_inpeak"].notna().values
        if args.frip_column in atac.obs:
            qc["atac_frip"] = reindex_values(atac.obs[args.frip_column].values, atac.obs_names, qc.index)
            frip_column = "atac_frip"
        elif args.atac_metrics:
            metrics = read_table_auto(args.atac_metrics)
            barcodes = metrics[args.atac_metrics_barcode_column].astype(str)
            if args.strip_barcode_suffix:
                barcodes = strip_suffix(barcodes)
            metrics = metrics.assign(_barcode=barcodes).set_index("_barcode")
            if args.atac_metrics_frip_column in metrics:
                qc["atac_frip"] = metrics[args.atac_metrics_frip_column].reindex(qc.index).values
                frip_column = "atac_frip"

    protein_marker: str | None = None
    rna_marker: str | None = None
    n_isotype_features = 0
    protein_specific_mean = 0.0
    protein_isotype_mean = 0.0
    qc["protein_present"] = False
    if args.protein_mod in mdata.mod:
        protein = mdata.mod[args.protein_mod]
        counts = counts_matrix(protein)
        qc["protein_umis"] = reindex_values(row_sum(counts), protein.obs_names, qc.index, 0)
        qc["protein_features_detected"] = reindex_values(row_nnz(counts), protein.obs_names, qc.index, 0)
        qc["protein_present"] = obs_bool_or_counts(protein, "protein_barcode_present", qc.index)
        qc["protein_module_call"] = obs_bool_or_counts(protein, "protein_module_call", qc.index)
        top_fraction = protein.obs["protein_top_feature_fraction"].values if "protein_top_feature_fraction" in protein.obs else np.zeros(protein.n_obs)
        qc["protein_top_feature_fraction"] = reindex_values(top_fraction, protein.obs_names, qc.index, 0)
        # Protein staining specificity: per-cell specific-antibody signal vs isotype-control
        # (background) signal. Isotype controls come from protein.var['is_isotype'] when the
        # builder marks them, otherwise from an IgG/isotype/control name match.
        labels = pd.Index(feature_labels(protein))
        if "is_isotype" in protein.var:
            iso_mask = np.asarray(protein.var["is_isotype"].astype(bool))
        else:
            iso_mask = np.asarray(labels.str.contains(r"igg|isotype|\bctrl\b|control", case=False, regex=True))
        dense_counts = dense(counts)
        iso_total = dense_counts[:, iso_mask].sum(axis=1) if iso_mask.any() else np.zeros(dense_counts.shape[0])
        spec_total = dense_counts[:, ~iso_mask].sum(axis=1)
        qc["protein_isotype_total"] = reindex_values(np.asarray(iso_total).ravel(), protein.obs_names, qc.index, 0)
        qc["protein_specific_total"] = reindex_values(np.asarray(spec_total).ravel(), protein.obs_names, qc.index, 0)
        n_isotype_features = int(iso_mask.sum())
        protein_specific_mean = float(np.mean(spec_total)) if spec_total.size else 0.0
        protein_isotype_mean = float(np.mean(iso_total)) if iso_mask.any() else 0.0
        rna_marker, protein_marker = select_marker_pair(rna, protein, args)
        protein_idx = find_feature_index(protein, protein_marker)
        if protein_idx is not None:
            qc[f"protein_marker_{protein_marker}"] = reindex_values(
                dense(counts)[:, protein_idx],
                protein.obs_names,
                qc.index,
                0,
            )
        rna_idx = find_feature_index(rna, rna_marker)
        if rna_idx is not None and rna_marker:
            qc[f"rna_marker_{rna_marker}"] = dense(counts_matrix(rna[:, [rna_idx]])).ravel()

    identity_mods = add_identity_qc(mdata, qc, args)
    four_factor_overlap = qc["atac_present"] & qc["protein_present"] & qc["identity_considered"]
    four_factor_assigned = qc["atac_present"] & qc["protein_present"] & qc["identity_assigned"]
    summary = {
        "n_rna_cells": int(qc.shape[0]),
        "n_atac_cells": int(qc["atac_present"].sum()),
        "n_protein_barcodes_present": int(qc["protein_present"].sum()),
        "n_identity_barcodes_present": int(qc["identity_considered"].sum()),
        "n_identity_assigned": int(qc["identity_assigned"].sum()),
        "n_four_factor_overlap": int(four_factor_overlap.sum()),
        "n_four_factor_assigned": int(four_factor_assigned.sum()),
        "identity_modalities": identity_mods,
        "protein_marker": protein_marker or "",
        "rna_marker": rna_marker or "",
        "frip_column": frip_column or "",
        "n_protein_isotype_features": n_isotype_features,
        "protein_specific_mean": round(protein_specific_mean, 4),
        "protein_isotype_mean": round(protein_isotype_mean, 4),
    }
    if isinstance(mdata.uns.get("multiome"), dict):
        summary["builder_multiome"] = {
            key: int(value) if isinstance(value, (np.integer, int)) else value
            for key, value in mdata.uns["multiome"].items()
            if key.startswith("n_") or key.endswith("_n_vars")
        }
    return qc, summary


def simple_embedding(adata: Any) -> pd.DataFrame:
    matrix = dense(counts_matrix(adata)).astype(float, copy=False)
    coords = np.zeros((adata.n_obs, 2), dtype=float)
    if matrix.size == 0:
        return pd.DataFrame(coords, index=adata.obs_names, columns=["umap1", "umap2"])

    totals = np.maximum(matrix.sum(axis=1, keepdims=True), 1.0)
    matrix = np.log1p(matrix / totals * 1e4)
    matrix = matrix - matrix.mean(axis=0, keepdims=True)
    try:
        u, singular_values, _ = np.linalg.svd(matrix, full_matrices=False)
        n_dims = min(2, u.shape[1], singular_values.shape[0])
        if n_dims:
            coords[:, :n_dims] = u[:, :n_dims] * singular_values[:n_dims]
    except np.linalg.LinAlgError:
        coords[:, 0] = np.arange(adata.n_obs, dtype=float)
    if not np.any(coords[:, 1]) and adata.n_obs > 1:
        coords[:, 1] = np.arange(adata.n_obs, dtype=float)
    return pd.DataFrame(coords, index=adata.obs_names, columns=["umap1", "umap2"])


def ensure_umap(mdata: md.MuData, qc: pd.DataFrame, args: argparse.Namespace) -> None:
    rna = mdata.mod[args.rna_mod]
    if args.umap_obsm in rna.obsm:
        umap = pd.DataFrame(rna.obsm[args.umap_obsm][:, :2], index=rna.obs_names, columns=["umap1", "umap2"])
    else:
        import scanpy as sc

        work = rna.copy()
        if work.n_obs < 3 or work.n_vars < 3:
            umap = simple_embedding(work)
        else:
            sc.pp.normalize_total(work, target_sum=1e4)
            sc.pp.log1p(work)
            sc.pp.highly_variable_genes(work, n_top_genes=min(2000, work.n_vars))
            work = work[:, work.var.highly_variable].copy()
            if work.n_vars < 3:
                umap = simple_embedding(rna)
            else:
                sc.pp.scale(work, max_value=10)
                n_comps = min(30, work.n_obs - 1, work.n_vars - 1)
                if n_comps < 2:
                    umap = simple_embedding(rna)
                else:
                    sc.tl.pca(work, n_comps=n_comps)
                    sc.pp.neighbors(work, n_neighbors=min(15, max(2, work.n_obs - 1)))
                    sc.tl.umap(work)
                    umap = pd.DataFrame(work.obsm["X_umap"][:, :2], index=work.obs_names, columns=["umap1", "umap2"])
    umap = umap.reindex(qc.index)
    qc["umap1"] = umap["umap1"]
    qc["umap2"] = umap["umap2"]


def add_empty_panel(fig: Any, row: int, col: int, text: str) -> None:
    fig.add_annotation(text=text, row=row, col=col, showarrow=False, x=0.5, y=0.5, xref="x domain", yref="y domain")


def add_identity_umap(fig: Any, qc: pd.DataFrame, row: int, col: int) -> None:
    sub = qc.dropna(subset=["umap1", "umap2"])
    top_labels = sub.loc[sub["identity_status"] == "singlet", "identity_label"].value_counts().head(8).index.tolist()
    for label in top_labels:
        cells = sub[sub["identity_label"] == label]
        fig.add_trace(
            go.Scattergl(x=cells["umap1"], y=cells["umap2"], mode="markers", name=label, marker=dict(size=3)),
            row,
            col,
        )
    background = sub[~sub["identity_label"].isin(top_labels)]
    fig.add_trace(
        go.Scattergl(
            x=background["umap1"],
            y=background["umap2"],
            mode="markers",
            name="other/none",
            marker=dict(size=2, color="#d0d0d0"),
        ),
        row,
        col,
    )


def render_qc(mdata: md.MuData, output_prefix: str, args: argparse.Namespace) -> dict[str, Any]:
    qc, summary = compute_qc(mdata, args)
    ensure_umap(mdata, qc, args)

    # Portrait 4x2 layout: all four colorbar panels in the right column so the
    # colorbars stack down the right edge (one per row) without overlapping; the
    # four bar/UMAP panels fill the left column. Identity legend goes below.
    ROWS, COLS, VS = 4, 2, 0.09
    _rh = (1.0 - (ROWS - 1) * VS) / ROWS
    row_center = {r: 1.0 - (r - 1) * (_rh + VS) - _rh / 2 for r in range(1, ROWS + 1)}

    def cbar(row: int, title: str) -> dict:
        return dict(title=dict(text=title, side="right", font=dict(size=12)),
                    tickfont=dict(size=11), x=1.005, xanchor="left",
                    y=row_center[row], yanchor="middle", len=_rh * 0.92, thickness=12)

    fig = make_subplots(
        rows=ROWS, cols=COLS, vertical_spacing=VS, horizontal_spacing=0.15,
        subplot_titles=(
            "Cells per modality / overlap", "GEX: UMIs vs genes",
            "Identity assignment status", "ATAC: in-peak depth vs peaks",
            "UMAP: identity", "Protein: UMIs vs detected ADTs",
            "Protein: signal vs isotype background", "UMAP: protein marker",
        ),
    )
    MS = 5

    # (1,1) cells per modality / overlap
    fig.add_trace(
        go.Bar(
            x=["RNA/ATAC", "protein", "identity present", "4-factor", "4-factor assigned"],
            y=[summary["n_atac_cells"], summary["n_protein_barcodes_present"],
               summary["n_identity_barcodes_present"], summary["n_four_factor_overlap"],
               summary["n_four_factor_assigned"]],
            showlegend=False,
            marker_color=["#4c78a8", "#f58518", "#54a24b", "#000000", "#b279a2"],
        ), 1, 1)
    fig.update_yaxes(title="cells", row=1, col=1)

    # (1,2) GEX: UMIs vs genes
    fig.add_trace(
        go.Scattergl(
            x=np.maximum(qc["gex_umis"], 1), y=np.maximum(qc["gex_genes"], 1), mode="markers",
            marker=dict(size=MS, color=qc["gex_pct_mito"], colorscale="Viridis", colorbar=cbar(1, "%mito")),
            showlegend=False,
        ), 1, 2)
    fig.update_xaxes(type="log", title="UMIs", row=1, col=2)
    fig.update_yaxes(type="log", title="genes", row=1, col=2)

    # (2,1) identity assignment status
    status_order = ["not_in_identity_set", "unassigned", "singlet", "multiplet"]
    status_disp = {"not_in_identity_set": "not in set", "unassigned": "unassigned",
                   "singlet": "singlet", "multiplet": "multiplet"}
    status_counts = qc["identity_status"].value_counts().reindex(status_order).fillna(0).astype(int)
    fig.add_trace(
        go.Bar(x=[status_disp.get(s, s) for s in status_counts.index], y=status_counts.values.tolist(),
               showlegend=False, marker_color="#4c78a8"),
        2, 1)

    # (2,2) ATAC: in-peak depth vs peaks
    if "atac_inpeak" in qc:
        atac_color = qc["atac_frip"] if "atac_frip" in qc else qc["atac_inpeak"]
        fig.add_trace(
            go.Scattergl(
                x=np.maximum(qc["atac_inpeak"].fillna(0), 1), y=np.maximum(qc["atac_npeaks"].fillna(0), 1), mode="markers",
                marker=dict(size=MS, color=atac_color, colorscale="Cividis",
                            colorbar=cbar(2, "FRiP" if "atac_frip" in qc else "depth")),
                showlegend=False,
            ), 2, 2)
        fig.update_xaxes(type="log", title="in-peak fragments", row=2, col=2)
        fig.update_yaxes(type="log", title="peaks", row=2, col=2)
    else:
        add_empty_panel(fig, 2, 2, "ATAC modality not present")

    # (3,1) UMAP: identity
    add_identity_umap(fig, qc, 3, 1)

    # (3,2) Protein: UMIs vs detected ADTs
    if "protein_umis" in qc:
        fig.add_trace(
            go.Scattergl(
                x=np.maximum(qc["protein_umis"], 1), y=np.maximum(qc["protein_features_detected"], 1), mode="markers",
                marker=dict(size=MS, color=qc["protein_top_feature_fraction"], colorscale="Teal", colorbar=cbar(3, "top fraction")),
                showlegend=False,
            ), 3, 2)
        fig.update_xaxes(type="log", title="protein UMIs", row=3, col=2)
        fig.update_yaxes(type="log", title="detected proteins", row=3, col=2)
    else:
        add_empty_panel(fig, 3, 2, "Protein modality not present")

    # (4,1) protein: signal vs isotype background
    if "protein_specific_total" in qc and "protein_isotype_total" in qc:
        present = qc["protein_present"].astype(bool)
        xx = np.maximum(qc.loc[present, "protein_isotype_total"].fillna(0), 1)
        yy = np.maximum(qc.loc[present, "protein_specific_total"].fillna(0), 1)
        fig.add_trace(go.Scattergl(x=xx, y=yy, mode="markers", marker=dict(size=MS, color="#54a24b"), showlegend=False), 4, 1)
        hi = float(max(yy.max() if len(yy) else 10, xx.max() if len(xx) else 10, 10))
        fig.add_trace(go.Scatter(x=[1, hi], y=[1, hi], mode="lines", line=dict(color="#888888", dash="dash"), showlegend=False), 4, 1)
        fig.update_xaxes(type="log", title="isotype-control UMIs (background)", row=4, col=1)
        fig.update_yaxes(type="log", title="specific-antibody UMIs", row=4, col=1)
    else:
        add_empty_panel(fig, 4, 1, "No protein modality for isotype QC")

    # (4,2) UMAP: protein marker
    protein_marker = summary["protein_marker"]
    marker_col = f"protein_marker_{protein_marker}" if protein_marker else ""
    sub = qc.dropna(subset=["umap1", "umap2"])
    if marker_col in qc:
        fig.add_trace(
            go.Scattergl(
                x=sub["umap1"], y=sub["umap2"], mode="markers",
                marker=dict(size=MS, color=qc.loc[sub.index, marker_col], colorscale="Magma", colorbar=cbar(4, protein_marker)),
                showlegend=False,
            ), 4, 2)
    else:
        add_empty_panel(fig, 4, 2, "No protein marker selected")

    overlap_rate = 100 * summary["n_four_factor_overlap"] / max(summary["n_rna_cells"], 1)
    assigned_rate = 100 * summary["n_four_factor_assigned"] / max(summary["n_rna_cells"], 1)
    fig.update_layout(
        height=1180, width=820, template="plotly_white",
        font=dict(size=15),
        legend=dict(itemsizing="constant", font=dict(size=12), orientation="h",
                    yanchor="top", y=-0.035, xanchor="center", x=0.5),
        margin=dict(l=78, r=112, t=92, b=72),
        title=dict(
            text=(f"{args.title}<br>cells={summary['n_rna_cells']}, "
                  f"4-factor={summary['n_four_factor_overlap']} ({overlap_rate:.1f}%), "
                  f"assigned={summary['n_four_factor_assigned']} ({assigned_rate:.1f}%)"),
            font=dict(size=15), x=0.5, xanchor="center", y=0.99, yanchor="top"),
    )
    fig.update_annotations(font_size=13)  # subplot titles

    out = Path(output_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(f"{out}.html", include_plotlyjs="cdn")
    Path(f"{out}.json").write_text(fig.to_json(), encoding="utf-8")
    Path(f"{out}.summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if not args.no_png:
        fig.write_image(f"{out}.png", scale=2)
    print(f">> QC summary: {summary}")
    print(f">> wrote {out}.html, {out}.json, {out}.summary.json" + ("" if args.no_png else f", {out}.png"))
    return summary


def main() -> None:
    args = parse_args()
    print(f">> reading MuData: {args.mudata}")
    mdata = md.read(args.mudata)
    render_qc(mdata, args.output_prefix, args)


if __name__ == "__main__":
    main()
