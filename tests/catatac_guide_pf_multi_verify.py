#!/usr/bin/env python3
"""Validate STAR pf-multi CAT-ATAC guide-arm smoke outputs."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def find_mex_dir(base: Path) -> Path:
    if (base / "matrix.mtx").is_file() or (base / "matrix.mtx.gz").is_file():
        return base
    for child in sorted(base.iterdir()):
        if child.is_dir() and (
            (child / "matrix.mtx").is_file() or (child / "matrix.mtx.gz").is_file()
        ):
            return child
    return base


def load_lines(path: Path) -> list[str]:
    if not path.is_file():
        return []
    return [line.split()[0] for line in path.read_text().splitlines() if line.strip()]


def load_features_tsv(path: Path) -> list[str]:
    """Return feature ids from production features.tsv (first column)."""
    if not path.is_file():
        raise FileNotFoundError(f"missing production features.tsv: {path}")
    ids: list[str] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        ids.append(line.split("\t")[0])
    return ids


def load_assign_features(path: Path) -> list[str]:
    """Assign-side feature axis (features.txt only on native split-read output)."""
    txt_path = path / "features.txt"
    if not txt_path.is_file():
        raise FileNotFoundError(f"missing assign features.txt: {txt_path}")
    return load_lines(txt_path)


def parse_matrix_dims(matrix_path: Path) -> tuple[int, int, int]:
    if not matrix_path.is_file():
        raise FileNotFoundError(f"missing matrix.mtx: {matrix_path}")
    with matrix_path.open() as fh:
        for line in fh:
            if line.startswith("%") or not line.strip():
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            return int(parts[0]), int(parts[1]), int(parts[2])
    raise ValueError(f"matrix header missing in {matrix_path}")


def load_gex_whitelist(path: Path) -> set[str]:
    out: set[str] = set()
    for line in path.read_text().splitlines():
        token = line.split()[0].upper()
        if token:
            out.add(token)
    return out


def load_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            mapping[parts[0].upper()] = parts[1].upper()
    return mapping


def load_mex_counts(
    mex_dir: Path,
    *,
    features: list[str],
    barcodes: list[str],
) -> dict[tuple[str, str], int]:
    features_upper = [x.upper() for x in features]
    barcodes_upper = [x.upper() for x in barcodes]
    counts: dict[tuple[str, str], int] = {}
    with (mex_dir / "matrix.mtx").open() as fh:
        saw_dims = False
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            if not saw_dims:
                saw_dims = True
                continue
            feat_idx, bc_idx, value = int(parts[0]), int(parts[1]), int(float(parts[2]))
            key = (features_upper[feat_idx - 1], barcodes_upper[bc_idx - 1])
            counts[key] = counts.get(key, 0) + value
    return counts


def validate_production_mex(
    mex_dir: Path,
    *,
    expected_features: int | None = None,
) -> dict:
    features_tsv = mex_dir / "features.tsv"
    barcodes_tsv = mex_dir / "barcodes.tsv"
    matrix_path = mex_dir / "matrix.mtx"

    feature_ids = load_features_tsv(features_tsv)
    barcodes = load_lines(barcodes_tsv)
    matrix_rows, matrix_cols, matrix_nnz = parse_matrix_dims(matrix_path)

    return {
        "feature_ids": feature_ids,
        "barcodes": barcodes,
        "feature_count": len(feature_ids),
        "unique_feature_names": len(set(feature_ids)),
        "barcode_count": len(barcodes),
        "matrix_rows": matrix_rows,
        "matrix_cols": matrix_cols,
        "matrix_nnz": matrix_nnz,
        "matrix_rows_match_features": matrix_rows == len(feature_ids),
        "matrix_cols_match_barcodes": matrix_cols == len(barcodes),
        "expected_features": expected_features,
        "feature_count_matches_expected": (
            expected_features is None or len(feature_ids) == expected_features
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--native-dir", required=True)
    parser.add_argument("--pf-multi-assign-dir", required=True)
    parser.add_argument("--log-out", required=True)
    parser.add_argument("--native-log", required=True)
    parser.add_argument("--output-map", required=True)
    parser.add_argument("--gex-whitelist", required=True)
    parser.add_argument("--expected-features", type=int, default=54)
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    native_mex = find_mex_dir(Path(args.native_dir))
    pf_mex = find_mex_dir(Path(args.pf_multi_assign_dir))
    log_text = Path(args.log_out).read_text(errors="replace")
    native_log = Path(args.native_log).read_text(errors="replace")
    combined_log = log_text + "\n" + native_log

    pf_mex_info = validate_production_mex(
        pf_mex, expected_features=args.expected_features
    )
    pf_features = pf_mex_info["feature_ids"]
    pf_barcodes_tsv = [x.upper() for x in pf_mex_info["barcodes"]]
    pf_barcodes_txt = [x.upper() for x in load_lines(pf_mex / "barcodes.txt")]

    native_features = load_assign_features(native_mex)
    native_barcodes = load_lines(native_mex / "barcodes.txt")

    gex_whitelist = load_gex_whitelist(Path(args.gex_whitelist))

    dup_warnings = {
        "HIC2_1": "duplicate feature name 'HIC2_1'" in combined_log,
        "HIC2_2": "duplicate feature name 'HIC2_2'" in combined_log,
    }

    gex_hits = sum(1 for bc in pf_barcodes_tsv if bc in gex_whitelist)
    atac_direct_gex = sum(1 for bc in pf_barcodes_txt if bc in gex_whitelist)
    mapped_diff = sum(
        1
        for atac_bc, gex_bc in zip(pf_barcodes_txt, pf_barcodes_tsv)
        if atac_bc != gex_bc
    )

    native_counts = load_mex_counts(
        native_mex,
        features=native_features,
        barcodes=native_barcodes,
    )
    # matrix.mtx column indices follow assign-side barcodes.txt (ATAC), not GEX barcodes.tsv
    pf_matrix_barcodes = load_lines(pf_mex / "barcodes.txt")
    pf_counts = load_mex_counts(
        pf_mex,
        features=pf_features,
        barcodes=pf_matrix_barcodes,
    )

    report = {
        "native_mex_dir": str(native_mex),
        "pf_multi_mex_dir": str(pf_mex),
        "duplicate_warnings": dup_warnings,
        "native_feature_count": len(native_features),
        "pf_production_mex": pf_mex_info,
        "gex_whitelist_hits": gex_hits,
        "atac_direct_gex_hits": atac_direct_gex,
        "mapped_diff_from_atac": mapped_diff,
        "matrix_equal": native_counts == pf_counts,
        "matrix_only_native": len(set(native_counts) - set(pf_counts)),
        "matrix_only_pf": len(set(pf_counts) - set(native_counts)),
        "matrix_value_mismatches": sum(
            1
            for key in set(native_counts) & set(pf_counts)
            if native_counts[key] != pf_counts[key]
        ),
    }
    Path(args.report).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))

    errors: list[str] = []
    if not dup_warnings["HIC2_1"] or not dup_warnings["HIC2_2"]:
        errors.append("expected duplicate feature warnings for HIC2_1 and HIC2_2")
    if not pf_mex_info["feature_count_matches_expected"]:
        errors.append(
            f"features.tsv row count expected {args.expected_features}, "
            f"got {pf_mex_info['feature_count']}"
        )
    if pf_mex_info["feature_count"] != pf_mex_info["unique_feature_names"]:
        errors.append("duplicate feature names remain in production features.tsv")
    if not pf_mex_info["matrix_rows_match_features"]:
        errors.append(
            f"matrix rows ({pf_mex_info['matrix_rows']}) != features.tsv rows "
            f"({pf_mex_info['feature_count']})"
        )
    if pf_mex_info["matrix_cols"] != len(pf_matrix_barcodes):
        errors.append(
            f"matrix cols ({pf_mex_info['matrix_cols']}) != barcodes.txt cols "
            f"({len(pf_matrix_barcodes)})"
        )
    if pf_mex_info["barcode_count"] != len(pf_matrix_barcodes):
        errors.append(
            f"barcodes.tsv count ({pf_mex_info['barcode_count']}) != barcodes.txt count "
            f"({len(pf_matrix_barcodes)})"
        )
    if pf_barcodes_tsv and gex_hits != len(pf_barcodes_tsv):
        errors.append("barcodes.tsv is not entirely in GEX whitelist space")
    if pf_barcodes_tsv and atac_direct_gex > max(1, len(pf_barcodes_txt) // 100):
        errors.append("too many ATAC-space barcodes already match GEX whitelist directly")
    if pf_barcodes_tsv and mapped_diff == 0:
        errors.append("barcodes.tsv did not remap ATAC barcodes to GEX namespace")
    if not report["matrix_equal"]:
        errors.append("pf-multi matrix does not match native split-read baseline")

    if errors:
        for err in errors:
            print(f"VERIFY FAIL: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
