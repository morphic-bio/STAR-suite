#!/usr/bin/env python3
"""Rebuild filtered MEX surfaces after EmptyDrops refiltering.

Given a STAR run directory whose Solo.out/GeneFull/filtered/barcodes.tsv has
been updated (e.g. by --runMode soloCellFiltering), this script:

1. Subsets outs/raw_feature_bc_matrix → outs/filtered_feature_bc_matrix
2. Copies new GeneFull filtered barcodes → Solo.out/Velocyto/filtered/barcodes.tsv
3. Rebuilds Solo.out/Velocyto/filtered/ matrices from Velocyto raw
4. Rebuilds outs/filtered_velocyto_feature_bc_matrix via prepare_velocyto_mex logic
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

from scipy.io import mmread, mmwrite
from scipy.sparse import issparse


def read_barcodes(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return [line.strip() for line in f if line.strip()]


def write_barcodes_gz(path: Path, barcodes: list[str]) -> None:
    with gzip.open(path, "wt") as f:
        for bc in barcodes:
            f.write(bc + "\n")


def write_barcodes_plain(path: Path, barcodes: list[str]) -> None:
    with open(path, "w") as f:
        for bc in barcodes:
            f.write(bc + "\n")


def read_matrix(path: Path):
    mat = mmread(path)
    if not issparse(mat):
        from scipy.sparse import csc_matrix
        mat = csc_matrix(mat)
    return mat.tocsc()


def write_matrix_gz(path: Path, mat) -> None:
    with gzip.open(path, "wb") as f:
        mmwrite(f, mat)


def write_matrix_plain(path: Path, mat) -> None:
    mmwrite(str(path), mat)


def strip_suffix(bc: str) -> str:
    return bc[:-2] if bc.endswith("-1") else bc


def add_suffix(bc: str) -> str:
    return bc if bc.endswith("-1") else bc + "-1"


def build_index(barcodes: list[str]) -> dict[str, int]:
    return {bc: i for i, bc in enumerate(barcodes)}


def subset_genefull_mex(run_dir: Path, new_filtered_barcodes: list[str]) -> None:
    """Subset outs/raw_feature_bc_matrix → outs/filtered_feature_bc_matrix."""
    raw_dir = run_dir / "outs" / "raw_feature_bc_matrix"
    out_dir = run_dir / "outs" / "filtered_feature_bc_matrix"

    # Back up old filtered
    backup = run_dir / "outs" / "filtered_feature_bc_matrix_pre_guarded_backup"
    if not backup.exists():
        shutil.copytree(out_dir, backup)

    raw_barcodes = read_barcodes(raw_dir / "barcodes.tsv.gz")
    raw_index = build_index([strip_suffix(bc) for bc in raw_barcodes])

    # Map new filtered barcodes (no suffix) to raw indices
    indices = []
    filtered_out = []
    for bc in new_filtered_barcodes:
        bc_stripped = strip_suffix(bc)
        if bc_stripped in raw_index:
            indices.append(raw_index[bc_stripped])
            filtered_out.append(add_suffix(bc_stripped))

    print(f"  GeneFull MEX: {len(indices)}/{len(new_filtered_barcodes)} barcodes found in raw")

    mat = read_matrix(raw_dir / "matrix.mtx.gz")
    filtered_mat = mat[:, indices]

    out_dir.mkdir(parents=True, exist_ok=True)
    write_barcodes_gz(out_dir / "barcodes.tsv.gz", filtered_out)
    shutil.copy2(raw_dir / "features.tsv.gz", out_dir / "features.tsv.gz")
    write_matrix_gz(out_dir / "matrix.mtx.gz", filtered_mat)

    print(f"  Wrote {out_dir} ({len(filtered_out)} barcodes, {filtered_mat.nnz} nnz)")


def subset_velocyto_solo(run_dir: Path, new_filtered_barcodes: list[str]) -> None:
    """Update Solo.out/Velocyto/filtered/ with new barcodes."""
    velo_raw = run_dir / "Solo.out" / "Velocyto" / "raw"
    velo_filt = run_dir / "Solo.out" / "Velocyto" / "filtered"

    # Back up old filtered
    backup = run_dir / "Solo.out" / "Velocyto" / "filtered_pre_guarded_backup"
    if not backup.exists():
        shutil.copytree(velo_filt, backup)

    raw_barcodes = read_barcodes(velo_raw / "barcodes.tsv")
    raw_index = build_index(raw_barcodes)

    indices = []
    filtered_out = []
    for bc in new_filtered_barcodes:
        bc_stripped = strip_suffix(bc)
        if bc_stripped in raw_index:
            indices.append(raw_index[bc_stripped])
            filtered_out.append(bc_stripped)

    print(f"  Velocyto Solo: {len(indices)}/{len(new_filtered_barcodes)} barcodes found in raw")

    # Write new filtered barcodes
    write_barcodes_plain(velo_filt / "barcodes.tsv", filtered_out)
    shutil.copy2(velo_raw / "features.tsv", velo_filt / "features.tsv")

    # Subset each layer
    for layer in ("spliced.mtx", "unspliced.mtx", "ambiguous.mtx"):
        mat = read_matrix(velo_raw / layer)
        filtered_mat = mat[:, indices]
        write_matrix_plain(velo_filt / layer, filtered_mat)
        print(f"    {layer}: {filtered_mat.shape} nnz={filtered_mat.nnz}")


def subset_velocyto_outs(run_dir: Path, new_filtered_barcodes: list[str]) -> None:
    """Subset outs/raw_velocyto_feature_bc_matrix → outs/filtered_velocyto_feature_bc_matrix."""
    raw_dir = run_dir / "outs" / "raw_velocyto_feature_bc_matrix"
    out_dir = run_dir / "outs" / "filtered_velocyto_feature_bc_matrix"

    # Back up old filtered
    backup = run_dir / "outs" / "filtered_velocyto_feature_bc_matrix_pre_guarded_backup"
    if not backup.exists():
        shutil.copytree(out_dir, backup)

    raw_barcodes = read_barcodes(raw_dir / "barcodes.tsv.gz")
    raw_index = build_index([strip_suffix(bc) for bc in raw_barcodes])

    indices = []
    filtered_out = []
    for bc in new_filtered_barcodes:
        bc_stripped = strip_suffix(bc)
        if bc_stripped in raw_index:
            indices.append(raw_index[bc_stripped])
            filtered_out.append(add_suffix(bc_stripped))

    print(f"  Velocyto outs MEX: {len(indices)}/{len(new_filtered_barcodes)} barcodes found in raw")

    out_dir.mkdir(parents=True, exist_ok=True)
    write_barcodes_gz(out_dir / "barcodes.tsv.gz", filtered_out)
    shutil.copy2(raw_dir / "features.tsv.gz", out_dir / "features.tsv.gz")

    for layer in ("matrix.mtx.gz", "spliced.mtx.gz", "unspliced.mtx.gz", "ambiguous.mtx.gz"):
        mat = read_matrix(raw_dir / layer)
        filtered_mat = mat[:, indices]
        write_matrix_gz(out_dir / layer, filtered_mat)
        print(f"    {layer}: {filtered_mat.shape} nnz={filtered_mat.nnz}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, help="STAR run directory")
    parser.add_argument("--sample-name", default="", help="Sample name for logging")
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    label = args.sample_name or run_dir.name

    new_barcodes_path = run_dir / "Solo.out" / "GeneFull" / "filtered" / "barcodes.tsv"
    new_filtered_barcodes = read_barcodes(new_barcodes_path)
    print(f"[{label}] New filtered barcodes: {len(new_filtered_barcodes)}")

    print(f"[{label}] Subsetting GeneFull outs MEX...")
    subset_genefull_mex(run_dir, new_filtered_barcodes)

    print(f"[{label}] Updating Velocyto Solo.out/filtered...")
    subset_velocyto_solo(run_dir, new_filtered_barcodes)

    print(f"[{label}] Subsetting Velocyto outs MEX...")
    subset_velocyto_outs(run_dir, new_filtered_barcodes)

    print(f"[{label}] Done.")


if __name__ == "__main__":
    main()
