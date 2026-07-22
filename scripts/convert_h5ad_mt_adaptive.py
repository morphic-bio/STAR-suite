#!/usr/bin/env python3
"""Retrofit an already-built downstream sample to the adaptive mt% guard.

CONVERSION script. Operates on an existing downstream sample directory that was
produced with the legacy strict-5% mt cutoff. It:

  1. reads final_counts.h5ad (all barcodes, with mt_pct / n_genes / singlet obs),
  2. recomputes the adaptive mt% threshold and the `filter` / `singlet_filtered`
     obs columns (legacy columns preserved as *_strict_mt5),
  3. rewrites final_counts.h5ad,
  4. rebuilds filtered_counts.h5ad and default_singlet_filtered_counts.h5ad
     (same subsetting rule as postprocess_downstream_filters.py),
  5. merges the mt% keys into adaptive_qc_threshold.json.

The rebuilt filtered_counts.h5ad / default_singlet_filtered_counts.h5ad are the
analyst-facing objects and should be re-staged to the delivery tree afterwards.

Note: final_counts.h5ad (all barcodes) is required — the delivery-tree
filtered_counts.h5ad only holds survivors and cannot be un-filtered. It lives
in the production-sample archive. See
docs/RUNBOOK_SCRNA_MT_ADAPTIVE_FILTER_20260518.md.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import anndata as ad

from scrna_mt_adaptive import (
    MT_FLOOR_DEFAULT,
    MT_N_MAD_DEFAULT,
    apply_mt_adaptive_filter,
    merge_threshold_json,
)

FINAL_NAME = "final_counts.h5ad"
FILTERED_NAME = "filtered_counts.h5ad"
DEFAULT_SINGLET_NAME = "default_singlet_filtered_counts.h5ad"
THRESHOLD_NAME = "adaptive_qc_threshold.json"


def sample_label(sample_dir: Path) -> str:
    """Human-readable sample name. The downstream dir itself is always named
    `downstream_genefull_velocyto_cellbender`, so use its parent (ES, PP1, ...)."""
    if sample_dir.name == "downstream_genefull_velocyto_cellbender":
        return sample_dir.parent.name
    return sample_dir.name


def convert_sample(
    sample_dir: Path,
    output_dir: Path,
    *,
    n_mad: float,
    mt_floor: float,
    backup: bool,
    dry_run: bool,
) -> dict:
    label = sample_label(sample_dir)
    final_h5ad = sample_dir / FINAL_NAME
    threshold_json = sample_dir / THRESHOLD_NAME
    if not final_h5ad.exists():
        raise FileNotFoundError(
            f"{final_h5ad} not found — conversion needs the all-barcodes "
            f"final_counts.h5ad (production-sample archive), not the delivery tree."
        )
    if not threshold_json.exists():
        raise FileNotFoundError(f"{threshold_json} not found")

    with open(threshold_json, "r", encoding="utf-8") as handle:
        thresholds = json.load(handle)
    try:
        min_genes = int(thresholds["min_genes"])
        max_genes = int(thresholds["effective_max_genes"])
    except KeyError as exc:
        raise KeyError(
            f"{threshold_json} missing min_genes / effective_max_genes"
        ) from exc

    adata = ad.read_h5ad(final_h5ad)
    record = apply_mt_adaptive_filter(
        adata,
        min_genes=min_genes,
        max_genes=max_genes,
        n_mad=n_mad,
        mt_floor=mt_floor,
    )
    print(f"[{label}] " + json.dumps(
        {k: record[k] for k in (
            "mt_pct_median", "mt_pct_mad", "mt_pct_raw_threshold",
            "mt_pct_threshold", "mt_pct_threshold_was_floored",
            "filter_cells_strict_mt5", "filter_cells_mt_adaptive",
            "mt_pct_flag",
        )},
        sort_keys=True,
    ))
    if record["mt_pct_flag"]:
        print(
            f"[{label}] WARNING: flagged for mt% review — "
            f"{record['mt_pct_flag_high_fraction']:.1%} of singlets exceed "
            f"{record['mt_pct_flag_high_pct']:.0f}% mt."
        )

    if dry_run:
        print(f"[{label}] dry-run: no files written")
        return record

    output_dir.mkdir(parents=True, exist_ok=True)
    out_final = output_dir / FINAL_NAME
    out_filtered = output_dir / FILTERED_NAME
    out_singlet = output_dir / DEFAULT_SINGLET_NAME
    out_json = output_dir / THRESHOLD_NAME

    if backup:
        for src in (final_h5ad, sample_dir / FILTERED_NAME,
                    sample_dir / DEFAULT_SINGLET_NAME, threshold_json):
            if src.exists():
                bak = src.with_suffix(src.suffix + ".pre_mt_adaptive")
                if not bak.exists():
                    shutil.copy2(src, bak)
                    print(f"[{label}] backed up {src.name} -> {bak.name}")

    # Rewrite final_counts and rebuild the analyst-facing views. The view rule
    # matches postprocess_downstream_filters.py.
    adata.write_h5ad(out_final)

    qc_mask = adata.obs["filter"].astype(bool).to_numpy()
    non_empty = adata.obs["non_empty"].astype(bool).to_numpy()
    singlet_filtered = adata.obs["singlet_filtered"].astype(bool).to_numpy()
    adata[qc_mask & non_empty].copy().write_h5ad(out_filtered)
    adata[singlet_filtered].copy().write_h5ad(out_singlet)

    if out_json != threshold_json and threshold_json.exists():
        shutil.copy2(threshold_json, out_json)
    merge_threshold_json(out_json, record)

    print(f"[{label}] wrote {out_final}")
    print(f"[{label}] wrote {out_filtered}  "
          f"({int((qc_mask & non_empty).sum())} cells)")
    print(f"[{label}] wrote {out_singlet}  "
          f"({int(singlet_filtered.sum())} cells)")
    print(f"[{label}] updated {out_json}")
    return record


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Retrofit downstream sample(s) to the adaptive mt% guard."
    )
    parser.add_argument(
        "--sample-dir", required=True, nargs="+",
        help="One or more downstream_genefull_velocyto_cellbender directories.",
    )
    parser.add_argument(
        "--output-dir",
        help="Write converted files here (per sample). Default: in place.",
    )
    parser.add_argument("--mt-floor", type=float, default=MT_FLOOR_DEFAULT)
    parser.add_argument("--n-mad", type=float, default=MT_N_MAD_DEFAULT)
    parser.add_argument(
        "--no-backup", action="store_true",
        help="Skip writing *.pre_mt_adaptive backups of the originals.",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Compute and report per sample without writing any file.",
    )
    args = parser.parse_args()

    summaries = {}
    for raw in args.sample_dir:
        sample_dir = Path(raw).resolve()
        label = sample_label(sample_dir)
        output_dir = (
            Path(args.output_dir).resolve() / label
            if args.output_dir else sample_dir
        )
        summaries[label] = convert_sample(
            sample_dir,
            output_dir,
            n_mad=args.n_mad,
            mt_floor=args.mt_floor,
            backup=not args.no_backup,
            dry_run=args.dry_run,
        )

    print("\n=== conversion summary ===")
    for name, rec in summaries.items():
        delta = (
            rec["filter_cells_mt_adaptive"] - rec["filter_cells_strict_mt5"]
            if rec.get("filter_cells_strict_mt5") is not None else None
        )
        print(
            f"{name:10}  mt_threshold={rec['mt_pct_threshold']:6.2f}  "
            f"strict5%={rec.get('filter_cells_strict_mt5')}  "
            f"adaptive={rec['filter_cells_mt_adaptive']}  "
            f"delta={'+' + str(delta) if delta and delta > 0 else delta}  "
            f"flag={rec['mt_pct_flag']}"
        )


if __name__ == "__main__":
    main()
