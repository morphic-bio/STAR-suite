#!/usr/bin/env python3
"""Extract a single feature_type surface from a Cell Ranger raw_feature_bc_matrix.

This rewrites the raw MEX into a smaller MEX containing only rows whose
features.tsv(.gz) third column matches the requested feature type.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return open(path, mode)


def resolve_required(base: Path, stem: str) -> Path:
    plain = base / stem
    gz = base / f"{stem}.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise SystemExit(f"Missing required file: {plain}(.gz)")


def count_lines(path: Path) -> int:
    n = 0
    with open_text(path, "r") as fh:
        for line in fh:
            if line.strip():
                n += 1
    return n


def write_lines_gz(src_path: Path, dst_path: Path) -> None:
    with open_text(src_path, "r") as src, gzip.open(dst_path, "wt") as dst:
        for line in src:
            dst.write(line)


def build_row_map(features_path: Path, feature_type: str, out_features_path: Path) -> dict[int, int]:
    row_map: dict[int, int] = {}
    selected = 0
    with open_text(features_path, "r") as src, gzip.open(out_features_path, "wt") as dst:
        for old_idx, raw in enumerate(src, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            if parts[2] != feature_type:
                continue
            selected += 1
            row_map[old_idx] = selected
            dst.write(line + "\n")
    if not row_map:
        raise SystemExit(f"No features matched feature_type={feature_type!r} in {features_path}")
    return row_map


def count_selected_nnz(matrix_path: Path, row_map: dict[int, int]) -> tuple[int, int, int]:
    rows = cols = nnz = 0
    dims_seen = False
    with open_text(matrix_path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            if not dims_seen:
                row_s, col_s, _ = line.split()
                rows = int(row_s)
                cols = int(col_s)
                dims_seen = True
                continue
            row_s, _, val_s = line.split()
            row = int(row_s)
            if row not in row_map:
                continue
            if int(val_s) != 0:
                nnz += 1
    if not dims_seen:
        raise SystemExit(f"Malformed Matrix Market file: {matrix_path}")
    return rows, cols, nnz


def write_filtered_matrix(matrix_path: Path, row_map: dict[int, int], n_rows: int, n_cols: int, nnz: int, out_matrix_path: Path) -> None:
    with gzip.open(out_matrix_path, "wt") as dst:
        dst.write("%%MatrixMarket matrix coordinate integer general\n")
        dst.write("%\n")
        dst.write(f"{n_rows} {n_cols} {nnz}\n")
        dims_seen = False
        with open_text(matrix_path, "r") as src:
            for raw in src:
                line = raw.strip()
                if not line or line.startswith("%"):
                    continue
                if not dims_seen:
                    dims_seen = True
                    continue
                row_s, col_s, val_s = line.split()
                row = int(row_s)
                new_row = row_map.get(row)
                if new_row is None:
                    continue
                value = int(val_s)
                if value == 0:
                    continue
                dst.write(f"{new_row} {col_s} {value}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-mex-dir", required=True, help="Cell Ranger raw_feature_bc_matrix directory")
    parser.add_argument("--feature-type", default="Gene Expression", help="Feature type to retain")
    parser.add_argument("--out-dir", required=True, help="Output MEX directory")
    args = parser.parse_args()

    input_dir = Path(args.input_mex_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    features_path = resolve_required(input_dir, "features.tsv")
    barcodes_path = resolve_required(input_dir, "barcodes.tsv")
    matrix_path = resolve_required(input_dir, "matrix.mtx")

    out_features = out_dir / "features.tsv.gz"
    out_barcodes = out_dir / "barcodes.tsv.gz"
    out_matrix = out_dir / "matrix.mtx.gz"

    row_map = build_row_map(features_path, args.feature_type, out_features)
    _, n_cols, nnz = count_selected_nnz(matrix_path, row_map)
    write_filtered_matrix(matrix_path, row_map, len(row_map), n_cols, nnz, out_matrix)
    write_lines_gz(barcodes_path, out_barcodes)
    barcode_count = count_lines(barcodes_path)

    if barcode_count != n_cols:
        raise SystemExit(
            f"Barcode count mismatch after extraction: barcodes={barcode_count} matrix_cols={n_cols}"
        )

    print(f"input_mex_dir={input_dir}")
    print(f"feature_type={args.feature_type}")
    print(f"selected_features={len(row_map)}")
    print(f"barcodes={barcode_count}")
    print(f"selected_nnz={nnz}")
    print(f"out_dir={out_dir}")


if __name__ == "__main__":
    main()
