#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from scipy.io import mmread, mmwrite
from scipy.sparse import issparse


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def read_rows(path: Path) -> list[list[str]]:
    with open_text(path, "rt") as handle:
        return [line.rstrip("\n").split("\t") for line in handle if line.strip()]


def read_single_column(path: Path) -> list[str]:
    return [row[0] for row in read_rows(path)]


def read_matrix(path: Path):
    matrix = mmread(path)
    if not issparse(matrix):
        raise TypeError(f"Expected sparse Matrix Market input: {path}")
    return matrix.tocsc()


def write_rows(path: Path, rows: list[list[str]]) -> None:
    with gzip.open(path, "wt") as handle:
        for row in rows:
            handle.write("\t".join(row))
            handle.write("\n")


def write_column(path: Path, values: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        for value in values:
            handle.write(value)
            handle.write("\n")


def write_matrix(path: Path, matrix) -> None:
    with gzip.open(path, "wb") as handle:
        mmwrite(handle, matrix)


def build_filtered_indices(raw_barcodes: list[str], filtered_barcodes: list[str]) -> list[int]:
    raw_index = {barcode: idx for idx, barcode in enumerate(raw_barcodes)}
    missing = [barcode for barcode in filtered_barcodes if barcode not in raw_index]
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(f"Filtered barcodes missing from raw barcode list: {preview}")
    return [raw_index[barcode] for barcode in filtered_barcodes]


def strip_barcode_suffix(barcode: str) -> str:
    if barcode.endswith("-1"):
        return barcode[:-2]
    return barcode


def ensure_parent(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_mex_dir(
    output_dir: Path,
    barcodes: list[str],
    features: list[list[str]],
    spliced,
    unspliced,
    ambiguous,
) -> dict[str, int]:
    ensure_parent(output_dir)

    total = (spliced + unspliced + ambiguous).tocsc()

    if total.shape[0] != len(features):
        raise ValueError(
            f"Feature count mismatch for {output_dir}: matrix rows={total.shape[0]} features={len(features)}"
        )
    if total.shape[1] != len(barcodes):
        raise ValueError(
            f"Barcode count mismatch for {output_dir}: matrix cols={total.shape[1]} barcodes={len(barcodes)}"
        )

    write_column(output_dir / "barcodes.tsv.gz", barcodes)
    write_rows(output_dir / "features.tsv.gz", features)
    write_matrix(output_dir / "matrix.mtx.gz", total)
    write_matrix(output_dir / "spliced.mtx.gz", spliced)
    write_matrix(output_dir / "unspliced.mtx.gz", unspliced)
    write_matrix(output_dir / "ambiguous.mtx.gz", ambiguous)

    return {
        "features": total.shape[0],
        "barcodes": total.shape[1],
        "nnz_total": int(total.nnz),
        "nnz_spliced": int(spliced.nnz),
        "nnz_unspliced": int(unspliced.nnz),
        "nnz_ambiguous": int(ambiguous.nnz),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Package STARsolo Velocyto raw outputs into stable raw/filtered MEX directories "
            "with a total matrix and per-layer matrices."
        )
    )
    parser.add_argument("--run-dir", required=True, help="STAR run directory containing Solo.out and outs")
    parser.add_argument(
        "--output-root",
        help="Directory that will receive raw_velocyto_feature_bc_matrix and filtered_velocyto_feature_bc_matrix "
        "(default: <run-dir>/outs)",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    output_root = Path(args.output_root).resolve() if args.output_root else run_dir / "outs"

    velocyto_raw_dir = run_dir / "Solo.out" / "Velocyto" / "raw"
    velocyto_filtered_dir = run_dir / "Solo.out" / "Velocyto" / "filtered"
    genefull_filtered_dir = run_dir / "Solo.out" / "GeneFull" / "filtered"

    raw_barcodes_path = velocyto_raw_dir / "barcodes.tsv"
    raw_features_path = velocyto_raw_dir / "features.tsv"
    filtered_barcodes_path = velocyto_filtered_dir / "barcodes.tsv"
    filtered_barcodes_source = "Solo.out/Velocyto/filtered/barcodes.tsv"
    if not filtered_barcodes_path.exists():
        genefull_filtered_barcodes_path = genefull_filtered_dir / "barcodes.tsv"
        if genefull_filtered_barcodes_path.exists():
            filtered_barcodes_path = genefull_filtered_barcodes_path
            filtered_barcodes_source = "Solo.out/GeneFull/filtered/barcodes.tsv"
    spliced_path = velocyto_raw_dir / "spliced.mtx"
    unspliced_path = velocyto_raw_dir / "unspliced.mtx"
    ambiguous_path = velocyto_raw_dir / "ambiguous.mtx"

    required_paths = [
        raw_barcodes_path,
        raw_features_path,
        filtered_barcodes_path,
        spliced_path,
        unspliced_path,
        ambiguous_path,
    ]
    missing_paths = [str(path) for path in required_paths if not path.exists()]
    if missing_paths:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing_paths))

    raw_barcodes = [strip_barcode_suffix(value) for value in read_single_column(raw_barcodes_path)]
    raw_features = read_rows(raw_features_path)
    filtered_barcodes = [strip_barcode_suffix(value) for value in read_single_column(filtered_barcodes_path)]

    spliced = read_matrix(spliced_path)
    unspliced = read_matrix(unspliced_path)
    ambiguous = read_matrix(ambiguous_path)

    shapes = {spliced.shape, unspliced.shape, ambiguous.shape}
    if len(shapes) != 1:
        raise ValueError(f"Velocyto layer shapes do not match: {shapes}")
    expected_shape = shapes.pop()
    if expected_shape != (len(raw_features), len(raw_barcodes)):
        raise ValueError(
            "Velocyto layer shape does not match Velocyto feature/barcode axes: "
            f"layers={expected_shape}, features={len(raw_features)}, barcodes={len(raw_barcodes)}"
        )

    filtered_indices = build_filtered_indices(raw_barcodes, filtered_barcodes)

    raw_dir = output_root / "raw_velocyto_feature_bc_matrix"
    filtered_dir = output_root / "filtered_velocyto_feature_bc_matrix"

    raw_summary = write_mex_dir(
        raw_dir,
        raw_barcodes,
        raw_features,
        spliced,
        unspliced,
        ambiguous,
    )

    filtered_summary = write_mex_dir(
        filtered_dir,
        filtered_barcodes,
        raw_features,
        spliced[:, filtered_indices],
        unspliced[:, filtered_indices],
        ambiguous[:, filtered_indices],
    )

    manifest = {
        "run_dir": str(run_dir),
        "output_root": str(output_root),
        "source": {
            "velocyto_raw_dir": str(velocyto_raw_dir),
            "velocyto_filtered_dir": str(velocyto_filtered_dir),
            "filtered_barcodes_source": filtered_barcodes_source,
        },
        "raw": raw_summary,
        "filtered": filtered_summary,
    }

    manifest_path = output_root / "velocyto_feature_bc_matrix_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Wrote {raw_dir}")
    print(f"Wrote {filtered_dir}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
