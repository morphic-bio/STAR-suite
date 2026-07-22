#!/usr/bin/env python3
"""Wiring checks for downsampled CAT-ATAC trimodal STAR smoke."""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from pathlib import Path


def existing_file(base: Path, names: list[str]) -> Path | None:
    for name in names:
        path = base / name
        if path.is_file():
            return path
    return None


def load_lines(path: Path) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [line.split()[0] for line in fh if line.strip()]


def find_mex_dir(base: Path) -> Path:
    if existing_file(base, ["matrix.mtx", "matrix.mtx.gz"]):
        return base
    if not base.is_dir():
        return base
    for child in sorted(base.iterdir()):
        if child.is_dir() and existing_file(child, ["matrix.mtx", "matrix.mtx.gz"]):
            return child
    return base


def parse_matrix_header(matrix_path: Path) -> tuple[int, int, int]:
    opener = gzip.open if str(matrix_path).endswith(".gz") else open
    with opener(matrix_path, "rt") as fh:
        for line in fh:
            if line.startswith("%") or not line.strip():
                continue
            parts = line.split()
            if len(parts) == 3:
                return int(parts[0]), int(parts[1]), int(parts[2])
    raise ValueError(f"no matrix header in {matrix_path}")


def read_mex_summary(
    mex_dir: Path,
    *,
    include_features: bool = False,
    include_barcodes: bool = False,
) -> dict:
    features_path = existing_file(mex_dir, ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"])
    barcodes_path = existing_file(mex_dir, ["barcodes.tsv", "barcodes.tsv.gz"])
    matrix_path = existing_file(mex_dir, ["matrix.mtx", "matrix.mtx.gz"])
    if features_path is None or barcodes_path is None or matrix_path is None:
        missing = []
        if features_path is None:
            missing.append("features.tsv[.gz]")
        if barcodes_path is None:
            missing.append("barcodes.tsv[.gz]")
        if matrix_path is None:
            missing.append("matrix.mtx[.gz]")
        raise FileNotFoundError(f"{mex_dir}: missing {', '.join(missing)}")

    feature_ids = load_lines(features_path)
    barcodes = [x.upper() for x in load_lines(barcodes_path)]
    rows, cols, nnz = parse_matrix_header(matrix_path)
    summary = {
        "path": str(mex_dir),
        "features_path": str(features_path),
        "barcodes_path": str(barcodes_path),
        "matrix_path": str(matrix_path),
        "feature_count": len(feature_ids),
        "unique_features": len(set(feature_ids)),
        "barcode_count": len(barcodes),
        "matrix_rows": rows,
        "matrix_cols": cols,
        "matrix_nnz": nnz,
        "rows_match_features": rows == len(feature_ids),
        "cols_match_barcodes": cols == len(barcodes),
    }
    if include_features:
        summary["feature_ids"] = feature_ids
    if include_barcodes:
        summary["barcodes"] = barcodes
    return summary


def verify_guide_mex(guide_mex: Path, gex_whitelist: Path, expected_features: int) -> dict:
    info = read_mex_summary(guide_mex, include_features=True, include_barcodes=True)
    barcodes = info["barcodes"]
    gex_wl = {ln.split()[0].upper() for ln in gex_whitelist.read_text().splitlines() if ln.split()}
    info.update({
        "gex_whitelist_hits": sum(1 for bc in barcodes if bc in gex_wl),
        "expected_features": expected_features,
        "feature_count_ok": info["feature_count"] == expected_features,
        "no_duplicate_features": info["feature_count"] == info["unique_features"],
    })
    return info


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--star-run", required=True)
    parser.add_argument("--gex-whitelist", required=True)
    parser.add_argument("--guide-library-id", default="catatac_guide_fixture")
    parser.add_argument("--expected-features", type=int, default=54)
    parser.add_argument("--expect-atac-peak-mex", choices=["yes", "no"], default="no")
    parser.add_argument("--report", required=True)
    args = parser.parse_args()

    star_run = Path(args.star_run)
    errors: list[str] = []
    report: dict = {"star_run": str(star_run)}

    log_final = star_run / "Log.final.out"
    log_out = star_run / "Log.out"
    report["log_final_exists"] = log_final.is_file()
    if not log_final.is_file():
        errors.append("missing Log.final.out")
    else:
        final_text = log_final.read_text(errors="replace")
        report["log_final_has_input_reads"] = "Number of input reads" in final_text
        if not report["log_final_has_input_reads"]:
            errors.append("Log.final.out does not contain final mapping metrics")

    report["log_out_exists"] = log_out.is_file()
    if not log_out.is_file():
        errors.append("missing Log.out")
    else:
        log_tail = log_out.read_text(errors="replace")[-16000:]
        report["finished_successfully"] = "ALL DONE!" in log_tail
        if not report["finished_successfully"]:
            errors.append("Log.out does not report ALL DONE")
        telemetry_matches = re.findall(
            r"ATAC permit telemetry: enabled=(yes|no) acquireCalls=(\d+)",
            log_tail,
        )
        if telemetry_matches:
            enabled, acquires = telemetry_matches[-1]
            report["atac_permit_telemetry_enabled"] = enabled == "yes"
            report["atac_permit_acquire_calls"] = int(acquires)
            if enabled == "yes" and int(acquires) <= 0:
                errors.append("ATAC permit telemetry enabled but acquireCalls is zero")

    gex_filtered = existing_file(
        star_run / "Solo.out/GeneFull/filtered", ["barcodes.tsv", "barcodes.tsv.gz"]
    )
    gex_raw = existing_file(
        star_run / "Solo.out/GeneFull/raw", ["barcodes.tsv", "barcodes.tsv.gz"]
    )
    gex_barcode_set: set[str] = set()
    gex_filtered_set = {x.upper() for x in load_lines(gex_filtered)} if gex_filtered else set()
    gex_raw_set = {x.upper() for x in load_lines(gex_raw)} if gex_raw else set()
    report["gex_filtered_barcodes_path"] = str(gex_filtered) if gex_filtered else None
    report["gex_raw_barcodes_path"] = str(gex_raw) if gex_raw else None
    report["gex_filtered_barcode_count"] = len(gex_filtered_set)
    report["gex_raw_barcode_count"] = len(gex_raw_set)
    gex_barcode_set = gex_filtered_set or gex_raw_set
    if not gex_barcode_set:
        errors.append("missing GEX GeneFull barcodes.tsv")

    atac_peaks = star_run / "atac/atac_peaks.narrowPeak"
    atac_bam = star_run / "atac_possorted.bam"
    atac_sidecar = star_run / "atac_fragments.bin"
    atac_sidecar_chroms = Path(str(atac_sidecar) + ".chroms.tsv")
    chromap_summary = star_run / "chromap_summary.csv"
    peak_mex = star_run / "atac/peak_mex"
    report["atac_peaks_exists"] = atac_peaks.is_file()
    report["atac_bam_exists"] = atac_bam.is_file()
    report["atac_bam_size"] = atac_bam.stat().st_size if atac_bam.is_file() else 0
    report["atac_sidecar_exists"] = atac_sidecar.is_file()
    report["atac_sidecar_size"] = atac_sidecar.stat().st_size if atac_sidecar.is_file() else 0
    report["atac_sidecar_chroms_exists"] = atac_sidecar_chroms.is_file()
    report["chromap_summary_exists"] = chromap_summary.is_file()
    report["atac_peak_mex_dir"] = str(peak_mex) if peak_mex.is_dir() else None
    if peak_mex.is_dir():
        try:
            report["atac_peak_mex"] = read_mex_summary(peak_mex, include_barcodes=True)
        except Exception as exc:
            errors.append(f"ATAC peak MEX unreadable: {exc}")
    if report["atac_bam_size"] <= 0:
        errors.append("missing or empty Chromap ATAC BAM")
    if report["atac_sidecar_size"] <= 0:
        errors.append("missing or empty Chromap ATAC sidecar")
    if not report["chromap_summary_exists"]:
        errors.append("missing Chromap summary")
    if args.expect_atac_peak_mex == "yes" and "atac_peak_mex" not in report:
        errors.append("missing requested ATAC peak MEX outputs")

    guide_base = star_run / "cr_assign/CRISPR_Guide_Capture" / args.guide_library_id
    guide_mex = find_mex_dir(guide_base)
    report["guide_mex_dir"] = str(guide_mex)
    guide_feature_ids: set[str] = set()
    guide_barcodes: set[str] = set()
    if existing_file(guide_mex, ["features.tsv", "features.tsv.gz"]) is None:
        errors.append(f"missing guide MEX under {guide_base}")
    else:
        try:
            guide_info = verify_guide_mex(
                guide_mex, Path(args.gex_whitelist), args.expected_features
            )
        except Exception as exc:
            errors.append(f"guide MEX unreadable: {exc}")
            guide_info = {}
        if guide_info:
            guide_feature_ids = set(guide_info.pop("feature_ids", []))
            guide_barcodes = set(guide_info.pop("barcodes", []))
            report["guide"] = guide_info
            if not guide_info["feature_count_ok"]:
                errors.append(
                    f"guide features.tsv count {guide_info['feature_count']} != {args.expected_features}"
                )
            if not guide_info["no_duplicate_features"]:
                errors.append("duplicate feature names in guide features.tsv")
            if not guide_info["rows_match_features"] or not guide_info["cols_match_barcodes"]:
                errors.append("guide matrix header does not match production axes")
            if guide_info["barcode_count"] and guide_info["gex_whitelist_hits"] != guide_info["barcode_count"]:
                errors.append("guide barcodes.tsv not entirely in GEX whitelist space")
            if gex_barcode_set:
                report["guide_gex_barcode_overlap"] = len(guide_barcodes & gex_barcode_set)

    raw_merged_dir = star_run / "outs/raw_feature_bc_matrix"
    filtered_merged_dir = star_run / "outs/filtered_feature_bc_matrix"
    if raw_merged_dir.is_dir():
        try:
            raw_merged = read_mex_summary(raw_merged_dir, include_features=True)
            raw_feature_ids = set(raw_merged.pop("feature_ids", []))
            report["raw_merged_mex"] = raw_merged
        except Exception as exc:
            errors.append(f"raw merged MEX unreadable: {exc}")
            raw_feature_ids = set()
    else:
        errors.append("missing raw merged feature_bc_matrix")
        raw_feature_ids = set()
    if filtered_merged_dir.is_dir():
        try:
            report["filtered_merged_mex"] = read_mex_summary(filtered_merged_dir)
        except Exception as exc:
            errors.append(f"filtered merged MEX unreadable: {exc}")

    if raw_feature_ids and guide_feature_ids:
        report["guide_features_in_raw_merged"] = len(guide_feature_ids & raw_feature_ids)
        if guide_feature_ids - raw_feature_ids:
            errors.append("not all guide features are present in raw merged MEX")

    if "atac_peak_mex" in report and guide_barcodes:
        atac_bcs = set(report["atac_peak_mex"]["barcodes"])
        report["guide_atac_peak_barcode_overlap"] = len(guide_barcodes & atac_bcs)
        report["atac_peak_mex"].pop("barcodes", None)

    Path(args.report).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if errors:
        for err in errors:
            print(f"VERIFY FAIL: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
