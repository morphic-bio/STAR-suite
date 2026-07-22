#!/usr/bin/env python3
"""Repair/standardize the CellBender denoised layer across UCSF production h5ad files.

Three sample states are handled:

1. Layer missing entirely (e.g. EBs1_1 repair pilot whose downstream was rebuilt
   without --reuse-cellbender) -> copy the layer from the stale final_counts.h5ad
   which still has it, then propagate to filtered subsets.

2. Layer named 'cellbender' instead of 'denoised' (samples processed by the
   remote CellBender watcher, which defaults to --cellbender-layer cellbender)
   -> rename in-place via h5py (instant, no data copy).

3. Layer already named 'denoised' -> no-op.

After fixing unfiltered_counts.h5ad, the layer is propagated to:
  - filtered_counts.h5ad            (barcode subset)
  - default_singlet_filtered_counts.h5ad  (barcode subset)
  - counts.h5ad                     (same barcodes)
  - final_counts.h5ad               (rebuilt as copy of unfiltered)
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

LAYER_NAME = "denoised"
OLD_LAYER_NAME = "cellbender"

ALL_SAMPLES = [
    "EBs1_1", "EBs1_2", "EBs1_3", "EBs1_4", "EBs1_5",
    "EBs2_1", "EBs2_2", "EBs2_3", "EBs2_4", "EBs2_5",
    "iPSC1_1", "iPSC1_2", "iPSC1_3",
    "iPSC2_1", "iPSC2_2", "iPSC2_3",
]

H5AD_FILES = [
    "unfiltered_counts.h5ad",
    "filtered_counts.h5ad",
    "default_singlet_filtered_counts.h5ad",
    "final_counts.h5ad",
    "counts.h5ad",
]


def get_layer_names(h5ad_path: Path) -> list[str]:
    with h5py.File(h5ad_path, "r") as f:
        layers = f.get("layers")
        return list(layers.keys()) if layers else []


def rename_layer_inplace(h5ad_path: Path, old_name: str, new_name: str) -> bool:
    with h5py.File(h5ad_path, "a") as f:
        layers = f.get("layers")
        if layers is not None and old_name in layers:
            f.move(f"layers/{old_name}", f"layers/{new_name}")
            return True
    return False


def copy_layer_same_obs(src_path: Path, dst_path: Path,
                        src_layer: str, dst_layer: str) -> None:
    """Copy a sparse layer between h5ad files with identical obs order."""
    with h5py.File(src_path, "r") as src_f, h5py.File(dst_path, "a") as dst_f:
        src_obs = src_f["obs"]["_index"][:]
        dst_obs = dst_f["obs"]["_index"][:]
        if not np.array_equal(src_obs, dst_obs):
            raise ValueError(
                f"Barcode mismatch: {src_path.name} ({len(src_obs)}) "
                f"vs {dst_path.name} ({len(dst_obs)})"
            )
        dst_layers = dst_f.require_group("layers")
        if dst_layer in dst_layers:
            del dst_layers[dst_layer]
        src_f.copy(f"layers/{src_layer}", dst_layers, name=dst_layer)


def subset_layer_to_target(src_path: Path, target_path: Path,
                           layer_name: str) -> None:
    """Subset a sparse layer from source to target by barcode matching."""
    with h5py.File(src_path, "r") as src_f:
        src_obs = np.asarray(src_f["obs"]["_index"][:])
        grp = src_f[f"layers/{layer_name}"]
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
        n_vars = src_f["var"]["_index"].shape[0]
        src_matrix = sp.csr_matrix((data, indices, indptr),
                                   shape=(len(src_obs), n_vars))

    src_bc_to_idx = {
        (bc.decode() if isinstance(bc, bytes) else bc): i
        for i, bc in enumerate(src_obs)
    }

    with h5py.File(target_path, "a") as dst_f:
        dst_obs = np.asarray(dst_f["obs"]["_index"][:])
        row_indices = []
        for bc in dst_obs:
            key = bc.decode() if isinstance(bc, bytes) else bc
            idx = src_bc_to_idx.get(key)
            if idx is None:
                raise KeyError(f"Target barcode {key} not in source")
            row_indices.append(idx)

        subset = src_matrix[row_indices, :].tocsr()

        dst_layers = dst_f.require_group("layers")
        if layer_name in dst_layers:
            del dst_layers[layer_name]
        grp = dst_layers.create_group(layer_name)
        grp.create_dataset("data", data=subset.data)
        grp.create_dataset("indices", data=subset.indices)
        grp.create_dataset("indptr", data=subset.indptr)
        grp.attrs["encoding-type"] = "csr_matrix"
        grp.attrs["encoding-version"] = "0.1.0"
        grp.attrs["shape"] = np.array(subset.shape, dtype=np.int64)


def repair_sample(ds: Path, sample: str, *, dry_run: bool) -> bool:
    """Repair one sample. Returns True on success."""
    unfiltered = ds / "unfiltered_counts.h5ad"
    filtered_h5 = ds / "filtered_counts.h5ad"
    default_singlet = ds / "default_singlet_filtered_counts.h5ad"
    final = ds / "final_counts.h5ad"
    counts = ds / "counts.h5ad"

    if not unfiltered.exists():
        print(f"  ERROR: {unfiltered} not found", file=sys.stderr)
        return False

    layers = get_layer_names(unfiltered)

    # --- Case 1: already correct ---
    if LAYER_NAME in layers:
        print(f"  {sample}: unfiltered already has '{LAYER_NAME}'")
        # Still check other h5ads for consistency
        for h5ad in [filtered_h5, default_singlet, final, counts]:
            if h5ad.exists():
                h_layers = get_layer_names(h5ad)
                if OLD_LAYER_NAME in h_layers and LAYER_NAME not in h_layers:
                    print(f"    renaming in {h5ad.name}")
                    if not dry_run:
                        rename_layer_inplace(h5ad, OLD_LAYER_NAME, LAYER_NAME)
        return True

    # --- Case 2: wrong name ---
    if OLD_LAYER_NAME in layers:
        print(f"  {sample}: renaming '{OLD_LAYER_NAME}' -> '{LAYER_NAME}'")
        if not dry_run:
            for h5ad in [unfiltered, filtered_h5, default_singlet, final, counts]:
                if h5ad.exists() and OLD_LAYER_NAME in get_layer_names(h5ad):
                    rename_layer_inplace(h5ad, OLD_LAYER_NAME, LAYER_NAME)
                    print(f"    renamed in {h5ad.name}")
        return True

    # --- Case 3: layer missing -> source from final_counts.h5ad ---
    print(f"  {sample}: layer MISSING in unfiltered")

    source_file = None
    source_layer = None
    for candidate, lname in [(final, OLD_LAYER_NAME), (final, LAYER_NAME)]:
        if candidate.exists() and lname in get_layer_names(candidate):
            source_file = candidate
            source_layer = lname
            break

    if source_file is None:
        print(f"  ERROR: no CellBender layer found anywhere for {sample}",
              file=sys.stderr)
        return False

    print(f"  sourcing '{source_layer}' from {source_file.name}")
    if dry_run:
        return True

    # Copy to unfiltered (same barcodes)
    copy_layer_same_obs(source_file, unfiltered, source_layer, LAYER_NAME)
    print(f"    added '{LAYER_NAME}' to unfiltered_counts.h5ad")

    # Copy to counts.h5ad (should have same barcodes)
    if counts.exists():
        try:
            copy_layer_same_obs(source_file, counts, source_layer, LAYER_NAME)
            print(f"    added '{LAYER_NAME}' to counts.h5ad")
        except ValueError:
            subset_layer_to_target(unfiltered, counts, LAYER_NAME)
            print(f"    subset '{LAYER_NAME}' to counts.h5ad")

    # Propagate to filtered subsets
    for target in [filtered_h5, default_singlet]:
        if target.exists():
            subset_layer_to_target(unfiltered, target, LAYER_NAME)
            print(f"    subset '{LAYER_NAME}' to {target.name}")

    # Rebuild final_counts.h5ad from unfiltered
    shutil.copy2(str(unfiltered), str(final))
    print(f"    rebuilt final_counts.h5ad from unfiltered")

    return True


def verify(run_root: Path) -> bool:
    print("\n=== Verification ===")
    all_ok = True
    for sample in ALL_SAMPLES:
        ds = run_root / "samples" / sample / "downstream_genefull_velocyto_cellbender"
        for fname in H5AD_FILES:
            h5ad = ds / fname
            if h5ad.exists():
                layers = get_layer_names(h5ad)
                if LAYER_NAME not in layers:
                    print(f"  FAIL  {sample}/{fname}  layers={layers}")
                    all_ok = False
                elif OLD_LAYER_NAME in layers:
                    print(f"  WARN  {sample}/{fname}  still has '{OLD_LAYER_NAME}'")
    if all_ok:
        print("  ALL OK: every h5ad has 'denoised' layer")
    return all_ok


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Repair CellBender denoised layer in UCSF production h5ad files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--run-root", required=True,
                        help="UCSF production run root directory")
    parser.add_argument("--samples", default=None,
                        help="Comma-separated sample IDs (default: all 16)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print actions without modifying files")
    args = parser.parse_args()

    run_root = Path(args.run_root)
    if not run_root.is_dir():
        print(f"ERROR: {run_root} is not a directory", file=sys.stderr)
        sys.exit(1)

    samples = args.samples.split(",") if args.samples else ALL_SAMPLES

    ok = True
    for sample in samples:
        ds = run_root / "samples" / sample / "downstream_genefull_velocyto_cellbender"
        if not ds.is_dir():
            print(f"  ERROR: {ds} not found", file=sys.stderr)
            ok = False
            continue
        if not repair_sample(ds, sample, dry_run=args.dry_run):
            ok = False

    if not args.dry_run:
        if not verify(run_root):
            ok = False

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
