#!/usr/bin/env python3
"""Apply the adaptive mt% QC guard to a downstream AnnData (new implementation).

NEW-IMPLEMENTATION pipeline step. Runs after combineFilters.py (which computes
`mt_pct`, `n_genes`, `singlet`, and the legacy strict-5% `filter`) and before
postprocess_downstream_filters.py. It:

  1. computes the per-sample adaptive mt% threshold over STAR-called singlets,
  2. recomputes the `filter` and `singlet_filtered` obs columns with the guard
     `mt_pct <= max(5%, median + n_mad*MAD)` (n_genes bounds unchanged),
  3. preserves the legacy columns as `filter_strict_mt5` /
     `singlet_filtered_strict_mt5`,
  4. merges the mt% keys into adaptive_qc_threshold.json.

postprocess_downstream_filters.py then builds filtered_counts.h5ad and
default_singlet_filtered_counts.h5ad from the recomputed `filter`, exactly as
today. The n_genes bounds (min_genes / effective_max_genes) are read from the
threshold JSON written earlier by compute_adaptive_qc_threshold.py.

See docs/RUNBOOK_SCRNA_MT_ADAPTIVE_FILTER_20260518.md for pipeline wiring.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad

from scrna_mt_adaptive import (
    MT_FLOOR_DEFAULT,
    MT_N_MAD_DEFAULT,
    apply_mt_adaptive_filter,
    merge_threshold_json,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Apply the adaptive mt% QC guard to a downstream AnnData."
    )
    parser.add_argument(
        "--input-h5ad", required=True,
        help="Unfiltered downstream h5ad from combineFilters.py (unfiltered_counts.h5ad).",
    )
    parser.add_argument(
        "--threshold-json", required=True,
        help="adaptive_qc_threshold.json (read for n_genes bounds, updated in place).",
    )
    parser.add_argument(
        "--output-h5ad",
        help="Where to write the updated h5ad (default: overwrite --input-h5ad).",
    )
    parser.add_argument("--mt-floor", type=float, default=MT_FLOOR_DEFAULT)
    parser.add_argument("--n-mad", type=float, default=MT_N_MAD_DEFAULT)
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Compute and print the threshold without writing any file.",
    )
    args = parser.parse_args()

    input_h5ad = Path(args.input_h5ad).resolve()
    threshold_json = Path(args.threshold_json).resolve()
    output_h5ad = Path(args.output_h5ad).resolve() if args.output_h5ad else input_h5ad

    if not threshold_json.exists():
        raise FileNotFoundError(f"threshold JSON not found: {threshold_json}")
    with open(threshold_json, "r", encoding="utf-8") as handle:
        thresholds = json.load(handle)
    try:
        min_genes = int(thresholds["min_genes"])
        max_genes = int(thresholds["effective_max_genes"])
    except KeyError as exc:
        raise KeyError(
            "threshold JSON missing n_genes bounds "
            "(min_genes / effective_max_genes) — run compute_adaptive_qc_threshold.py first"
        ) from exc

    adata = ad.read_h5ad(input_h5ad)
    record = apply_mt_adaptive_filter(
        adata,
        min_genes=min_genes,
        max_genes=max_genes,
        n_mad=args.n_mad,
        mt_floor=args.mt_floor,
    )

    print(json.dumps(record, indent=2, sort_keys=True))
    if record["mt_pct_flag"]:
        print(
            f"WARNING: sample flagged for mt% review — "
            f"{record['mt_pct_flag_high_fraction']:.1%} of singlets exceed "
            f"{record['mt_pct_flag_high_pct']:.0f}% mt."
        )

    if args.dry_run:
        print("dry-run: no files written")
        return

    adata.write_h5ad(output_h5ad)
    merge_threshold_json(threshold_json, record)
    print(f"Wrote {output_h5ad}")
    print(f"Updated {threshold_json}")


if __name__ == "__main__":
    main()
