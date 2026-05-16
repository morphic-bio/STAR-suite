#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t")
    return open(path, mode, encoding="utf-8")


def resolve_required(base: Path, stem: str) -> Path:
    plain = base / stem
    gz = base / f"{stem}.gz"
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise FileNotFoundError(f"Missing required file: {plain}(.gz)")


def count_rows(path: Path) -> int:
    with open_text(path, "r") as handle:
        return sum(1 for line in handle if line.strip())


def read_barcodes(path: Path) -> list[str]:
    values: list[str] = []
    with open_text(path, "r") as handle:
        for raw in handle:
            line = raw.strip()
            if line:
                values.append(line.split("\t", 1)[0])
    return values


def ensure_unique(values: list[str], label: str) -> None:
    seen: set[str] = set()
    dupes: list[str] = []
    for value in values:
        if value in seen:
            dupes.append(value)
            if len(dupes) >= 5:
                break
        seen.add(value)
    if dupes:
        raise ValueError(f"Duplicate {label}: {', '.join(dupes)}")


def matrix_shape(path: Path) -> tuple[int, int, int]:
    with open_text(path, "r") as handle:
        header = handle.readline()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError(f"Not a Matrix Market file: {path}")
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) != 3:
                raise ValueError(f"Malformed Matrix Market size line in {path}: {line!r}")
            rows, cols, nnz = (int(part) for part in parts)
            return rows, cols, nnz
    raise ValueError(f"Missing Matrix Market size line in {path}")


def copy_text_gz(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    with open_text(source, "r") as src, gzip.open(target, "wt") as dst:
        shutil.copyfileobj(src, dst)


def package_mex(source_dir: Path, output_dir: Path, label: str) -> dict[str, int | str]:
    barcodes_path = resolve_required(source_dir, "barcodes.tsv")
    features_path = resolve_required(source_dir, "features.tsv")
    matrix_path = resolve_required(source_dir, "matrix.mtx")

    barcodes = read_barcodes(barcodes_path)
    ensure_unique(barcodes, f"{label} barcodes")
    n_barcodes = len(barcodes)
    n_features = count_rows(features_path)
    matrix_rows, matrix_cols, matrix_nnz = matrix_shape(matrix_path)

    if (matrix_rows, matrix_cols) != (n_features, n_barcodes):
        raise ValueError(
            f"{label} matrix shape mismatch: matrix={(matrix_rows, matrix_cols)} "
            f"features/barcodes={(n_features, n_barcodes)}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    copy_text_gz(barcodes_path, output_dir / "barcodes.tsv.gz")
    copy_text_gz(features_path, output_dir / "features.tsv.gz")
    copy_text_gz(matrix_path, output_dir / "matrix.mtx.gz")

    return {
        "source_dir": str(source_dir),
        "features": n_features,
        "barcodes": n_barcodes,
        "nnz": matrix_nnz,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Package STARsolo Solo.out/GeneFull raw and filtered MEX into "
            "CellRanger-style outs/raw_feature_bc_matrix and "
            "outs/filtered_feature_bc_matrix directories."
        )
    )
    parser.add_argument("--run-dir", required=True, help="STAR run directory containing Solo.out/GeneFull")
    parser.add_argument(
        "--output-root",
        help="Directory receiving raw_feature_bc_matrix and filtered_feature_bc_matrix (default: <run-dir>/outs)",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    output_root = Path(args.output_root).resolve() if args.output_root else run_dir / "outs"
    genefull_dir = run_dir / "Solo.out" / "GeneFull"
    raw_dir = genefull_dir / "raw"
    filtered_dir = genefull_dir / "filtered"

    if not raw_dir.exists():
        raise FileNotFoundError(f"Missing GeneFull raw directory: {raw_dir}")
    if not filtered_dir.exists():
        raise FileNotFoundError(f"Missing GeneFull filtered directory: {filtered_dir}")

    raw_summary = package_mex(raw_dir, output_root / "raw_feature_bc_matrix", "GeneFull raw")
    filtered_summary = package_mex(
        filtered_dir,
        output_root / "filtered_feature_bc_matrix",
        "GeneFull filtered",
    )

    raw_barcodes = set(read_barcodes(resolve_required(raw_dir, "barcodes.tsv")))
    filtered_barcodes = read_barcodes(resolve_required(filtered_dir, "barcodes.tsv"))
    missing = [barcode for barcode in filtered_barcodes if barcode not in raw_barcodes]
    if missing:
        raise ValueError(
            "GeneFull filtered barcodes missing from raw barcodes: "
            + ", ".join(missing[:5])
        )

    manifest = {
        "run_dir": str(run_dir),
        "output_root": str(output_root),
        "source": {
            "genefull_raw_dir": str(raw_dir),
            "genefull_filtered_dir": str(filtered_dir),
        },
        "raw": raw_summary,
        "filtered": filtered_summary,
    }
    manifest_path = output_root / "gene_full_feature_bc_matrix_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print(f"Wrote {output_root / 'raw_feature_bc_matrix'}")
    print(f"Wrote {output_root / 'filtered_feature_bc_matrix'}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
