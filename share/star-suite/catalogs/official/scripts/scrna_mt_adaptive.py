#!/usr/bin/env python3
"""Shared helpers for the adaptive mitochondrial-percentage (mt%) QC guard.

A cell passes the mt% guard when

    mt_pct <= mt_threshold ,   with
    mt_threshold = max(mt_floor, median(mt_pct) + n_mad * MAD(mt_pct))

`mt_floor` (default 5.0%) is a strict guard: any cell at or below it is always
kept, so the rule never tightens below the historical fixed cutoff. The
MAD-based term is a per-sample soft guard that adapts to populations with a
biologically elevated baseline mt% (e.g. MSK ES undifferentiated cells, whose
median mt% already sits at ~5%). Because the floor is folded into the single
`mt_threshold` value, downstream code only needs `mt_pct <= mt_threshold`.

MAD is the raw median absolute deviation, matching the n_genes convention in
compute_adaptive_qc_threshold.py (median + n_mad * raw_MAD). Raw MAD x 3 is
already on the strict side of the field convention (3 scaled MADs), so no
scaling factor is applied.

This module is imported by both the new-implementation script
(apply_adaptive_mt_filter.py) and the conversion script
(convert_h5ad_mt_adaptive.py) so the two paths compute identical thresholds.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

MT_FLOOR_DEFAULT = 5.0
MT_N_MAD_DEFAULT = 3.0
# Sample-level review tripwire (never a per-cell cut): flag a sample when more
# than MT_FLAG_HIGH_FRAC of its cells sit above MT_FLAG_HIGH_PCT mt%.
MT_FLAG_HIGH_PCT_DEFAULT = 20.0
MT_FLAG_HIGH_FRAC_DEFAULT = 0.10


def finite_mt(mt_pct) -> np.ndarray:
    """Return the finite mt_pct values as a flat float array."""
    values = np.asarray(mt_pct, dtype=float).ravel()
    return values[np.isfinite(values)]


def compute_mt_threshold(
    mt_pct,
    n_mad: float = MT_N_MAD_DEFAULT,
    mt_floor: float = MT_FLOOR_DEFAULT,
) -> dict:
    """Combined strict-floor + MAD soft-guard mt% threshold for one sample."""
    values = finite_mt(mt_pct)
    if values.size == 0:
        raise ValueError("No finite mt_pct values for adaptive thresholding")
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    raw_threshold = float(median + n_mad * mad)
    threshold = float(max(mt_floor, raw_threshold))
    return {
        "mt_pct_median": median,
        "mt_pct_mad": mad,
        "mt_pct_n_mad": float(n_mad),
        "mt_pct_floor": float(mt_floor),
        "mt_pct_raw_threshold": raw_threshold,
        "mt_pct_threshold": threshold,
        "mt_pct_threshold_was_floored": bool(threshold > raw_threshold),
        "mt_pct_cells": int(values.size),
    }


def mt_sample_flag(
    mt_pct,
    high_pct: float = MT_FLAG_HIGH_PCT_DEFAULT,
    high_frac: float = MT_FLAG_HIGH_FRAC_DEFAULT,
) -> dict:
    """Sample-level review tripwire. Not a per-cell cut.

    Flags a sample when a large fraction of cells sit far into the high-mt
    tail, which is where the MAD soft guard could otherwise drift too high. A
    flagged sample should go to human QC review rather than be auto-rescued.
    """
    values = finite_mt(mt_pct)
    frac_high = float((values > high_pct).mean()) if values.size else 0.0
    return {
        "mt_pct_flag": bool(frac_high > high_frac),
        "mt_pct_flag_high_pct": float(high_pct),
        "mt_pct_flag_high_fraction": frac_high,
        "mt_pct_flag_high_fraction_limit": float(high_frac),
    }


def qc_filter_mask(
    n_genes,
    mt_pct,
    *,
    min_genes: int,
    max_genes: int,
    mt_threshold: float,
) -> np.ndarray:
    """Per-cell QC mask: n_genes in [min, max] AND mt_pct <= mt_threshold.

    `mt_threshold` already folds in the strict floor (see compute_mt_threshold),
    so `mt_pct <= mt_threshold` IS the combined strict + MAD guard. The n_genes
    bounds are unchanged from the existing adaptive filter.
    """
    n_genes = np.asarray(n_genes, dtype=float).ravel()
    mt_pct = np.asarray(mt_pct, dtype=float).ravel()
    return (
        (n_genes >= min_genes)
        & (n_genes <= max_genes)
        & (mt_pct <= mt_threshold)
    )


def apply_mt_adaptive_filter(
    adata,
    *,
    min_genes: int,
    max_genes: int,
    n_mad: float = MT_N_MAD_DEFAULT,
    mt_floor: float = MT_FLOOR_DEFAULT,
) -> dict:
    """Recompute `filter` / `singlet_filtered` obs using the adaptive mt guard.

    The threshold is computed over STAR-called singlets (obs['singlet']),
    matching compute_adaptive_qc_threshold.py, then applied to every barcode.
    The pre-existing strict-5% columns are preserved as `filter_strict_mt5`
    and `singlet_filtered_strict_mt5` so the conversion is auditable and
    reversible. Returns the threshold record (extended with a summary).
    """
    for col in ("mt_pct", "n_genes", "singlet"):
        if col not in adata.obs.columns:
            raise KeyError(f"input h5ad missing required obs column: {col}")

    singlet = adata.obs["singlet"].astype(bool).to_numpy()
    mt_all = adata.obs["mt_pct"].astype(float).to_numpy()
    n_genes = adata.obs["n_genes"].astype(float).to_numpy()
    if int(singlet.sum()) == 0:
        raise ValueError("No singlet cells available for adaptive mt thresholding")

    record = compute_mt_threshold(mt_all[singlet], n_mad=n_mad, mt_floor=mt_floor)
    record.update(mt_sample_flag(mt_all[singlet]))
    record["mt_pct_source"] = "singlets"

    new_filter = qc_filter_mask(
        n_genes, mt_all,
        min_genes=min_genes, max_genes=max_genes,
        mt_threshold=record["mt_pct_threshold"],
    )

    # Preserve the legacy strict-5% columns once (idempotent).
    strict_filter = None
    if "filter_strict_mt5" in adata.obs.columns:
        strict_filter = adata.obs["filter_strict_mt5"].astype(bool).to_numpy()
    elif "filter" in adata.obs.columns:
        strict_filter = adata.obs["filter"].astype(bool).to_numpy()
        adata.obs["filter_strict_mt5"] = strict_filter

    strict_singlet_filter = None
    if "singlet_filtered_strict_mt5" in adata.obs.columns:
        strict_singlet_filter = (
            adata.obs["singlet_filtered_strict_mt5"].astype(bool).to_numpy()
        )
    elif "singlet_filtered" in adata.obs.columns:
        strict_singlet_filter = adata.obs["singlet_filtered"].astype(bool).to_numpy()
        adata.obs["singlet_filtered_strict_mt5"] = strict_singlet_filter

    adata.obs["filter"] = new_filter
    adata.obs["singlet_filtered"] = singlet & new_filter

    record["min_genes"] = int(min_genes)
    record["max_genes"] = int(max_genes)
    record["filter_cells_strict_mt5"] = (
        int(strict_filter.sum()) if strict_filter is not None else None
    )
    record["singlet_filtered_cells_strict_mt5"] = (
        int(strict_singlet_filter.sum()) if strict_singlet_filter is not None else None
    )
    record["filter_cells_mt_adaptive"] = int(new_filter.sum())
    record["singlet_filtered_cells_mt_adaptive"] = int((singlet & new_filter).sum())
    adata.uns["mt_adaptive_filter"] = dict(record)
    return record


def merge_threshold_json(path, updates: dict) -> dict:
    """Merge `updates` into an adaptive_qc_threshold.json, preserving keys."""
    path = Path(path)
    record: dict = {}
    if path.exists():
        with open(path, "r", encoding="utf-8") as handle:
            record = json.load(handle)
    record.update(updates)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
    return record
