#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import tarfile
from pathlib import Path
from typing import Iterable


DEFAULT_ORACLE_DIR = Path("/mnt/pikachu/JAX_scRNAseq02/cellranger-logs")
REQUIRED_MEX_FILES = ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")
OCM_IDS = ("OB1", "OB2", "OB3", "OB4")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt", encoding="utf-8")


def count_lines(path: Path) -> int:
    with open_text(path) as handle:
        return sum(1 for line in handle if line.strip())


def read_single_column(path: Path) -> list[str]:
    with open_text(path) as handle:
        values = [line.rstrip("\n").split(",")[0] for line in handle if line.strip()]
    if values and values[0].lower() in {"barcode", "barcodes"}:
        values = values[1:]
    return values


def normalize_barcode(barcode: str) -> str:
    return barcode[:-2] if barcode.endswith("-1") else barcode


def load_cells_per_tag(path: Path) -> dict[str, list[str]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return {key: [normalize_barcode(str(value)) for value in values] for key, values in data.items()}


def parse_multi_config(path: Path) -> dict[str, list[str]]:
    section = ""
    headers: list[str] = []
    samples: dict[str, list[str]] = {}

    with path.open("r", encoding="utf-8", newline="") as handle:
        for raw_row in csv.reader(handle):
            row = [cell.strip() for cell in raw_row]
            if not row or not any(row):
                continue
            first = row[0]
            if first.startswith("[") and first.endswith("]"):
                section = first.strip("[]")
                headers = []
                continue
            if section != "samples":
                continue
            if not headers:
                headers = row
                continue
            record = dict(zip(headers, row))
            sample_id = record.get("sample_id", "")
            ocm_ids = record.get("ocm_barcode_ids", "")
            if sample_id and ocm_ids:
                samples[sample_id] = [value.strip() for value in ocm_ids.split("|") if value.strip()]

    if not samples:
        raise ValueError(f"No [samples] entries found in {path}")
    return samples


def list_mri_paths(path: Path) -> set[str]:
    with tarfile.open(path, "r:gz") as tar:
        return set(tar.getnames())


def require_path(path: Path, errors: list[str]) -> None:
    if not path.exists():
        errors.append(f"missing: {path}")


def require_mex_dir(path: Path, errors: list[str]) -> dict[str, int]:
    stats: dict[str, int] = {}
    for filename in REQUIRED_MEX_FILES:
        require_path(path / filename, errors)
    barcode_path = path / "barcodes.tsv.gz"
    if barcode_path.exists():
        stats["barcodes"] = count_lines(barcode_path)
    return stats


def sample_filtered_count(path: Path) -> int:
    return len(read_single_column(path))


def proportions(counts: dict[str, int]) -> dict[str, float]:
    total = sum(counts.values())
    if total == 0:
        return {key: 0.0 for key in counts}
    return {key: value / total for key, value in counts.items()}


def max_abs_delta(left: dict[str, float], right: dict[str, float], keys: Iterable[str]) -> float:
    return max(abs(left.get(key, 0.0) - right.get(key, 0.0)) for key in keys)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Validate STAR JAX scRNAseq02 OCM output against the 25E32-L3 "
            "Cell Ranger multi oracle structure and OCM barcode counts. This "
            "does not assert matrix-value parity because the oracle uses "
            "exonic Gene counts while the production STAR surface uses GeneFull."
        )
    )
    parser.add_argument("--star-run-dir", required=True, help="STAR smoke run directory containing outs/")
    parser.add_argument("--oracle-dir", default=str(DEFAULT_ORACLE_DIR), help="Cell Ranger oracle/log directory")
    parser.add_argument("--report-json", help="Write validation report JSON")
    parser.add_argument(
        "--min-total-cells",
        type=int,
        default=1,
        help="Minimum total cells across OCM tags for a smoke run (default: 1)",
    )
    parser.add_argument(
        "--max-tag-proportion-delta",
        type=float,
        default=0.25,
        help=(
            "Maximum absolute OCM tag proportion delta before warning "
            "(default: 0.25). Use --strict-tag-proportion-delta to fail on this."
        ),
    )
    parser.add_argument(
        "--strict-tag-proportion-delta",
        action="store_true",
        help=(
            "Fail when --max-tag-proportion-delta is exceeded. By default this "
            "is a warning because the oracle is exonic Gene while production "
            "STAR outputs use GeneFull/EmptyDrops_CR cell calls."
        ),
    )
    parser.add_argument(
        "--min-oracle-overlap",
        type=float,
        default=0.0,
        help=(
            "Optional minimum fraction of STAR-called OCM barcodes present in "
            "Cell Ranger oracle tags. This is precision against the oracle and "
            "is usually left unset for GeneFull production outputs."
        ),
    )
    parser.add_argument(
        "--min-oracle-recall",
        type=float,
        default=0.90,
        help=(
            "Minimum fraction of Cell Ranger oracle barcodes recovered in the "
            "matching STAR OCM tag before failing (default: 0.90)."
        ),
    )
    args = parser.parse_args()

    star_run_dir = Path(args.star_run_dir).resolve()
    oracle_dir = Path(args.oracle_dir).resolve()
    outs_dir = star_run_dir / "outs"
    report_path = Path(args.report_json).resolve() if args.report_json else star_run_dir / "ocm_oracle_validation.json"

    config_path = oracle_dir / "config.csv"
    oracle_cells_path = oracle_dir / "cells_per_tag.json"
    mri_path = oracle_dir / "25E32-L3_Day4-pool-1.mri.tgz"

    errors: list[str] = []
    warnings: list[str] = []

    for path in [outs_dir, config_path, oracle_cells_path, mri_path]:
        require_path(path, errors)
    if errors:
        report = {"status": "FAIL", "errors": errors, "warnings": warnings}
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        raise SystemExit("FAIL: required inputs are missing; see " + str(report_path))

    samples = parse_multi_config(config_path)
    oracle_cells = load_cells_per_tag(oracle_cells_path)
    mri_paths = list_mri_paths(mri_path)

    expected_archive_markers = [
        "25E32-L3_Day4-pool-1/outs/multi/multiplexing_analysis/cells_per_tag.json",
    ]
    for sample in samples:
        expected_archive_markers.extend(
            [
                f"25E32-L3_Day4-pool-1/outs/per_sample_outs/{sample}/count/sample_filtered_feature_bc_matrix/barcodes.tsv.gz",
                f"25E32-L3_Day4-pool-1/outs/per_sample_outs/{sample}/count/sample_raw_feature_bc_matrix/features.tsv.gz",
                f"25E32-L3_Day4-pool-1/outs/per_sample_outs/{sample}/count/sample_filtered_barcodes.csv",
            ]
        )
    missing_archive_markers = [path for path in expected_archive_markers if path not in mri_paths]
    if missing_archive_markers:
        warnings.append(
            "oracle MRI archive is missing expected layout markers: "
            + ", ".join(missing_archive_markers[:5])
        )

    multi_count_dir = outs_dir / "multi" / "count"
    multi_raw_dir = multi_count_dir / "raw_feature_bc_matrix"
    multi_mux_dir = outs_dir / "multi" / "multiplexing_analysis"
    output_cells_path = multi_mux_dir / "cells_per_tag.json"

    multi_raw_stats = require_mex_dir(multi_raw_dir, errors)
    require_path(output_cells_path, errors)

    sample_reports: dict[str, dict[str, object]] = {}
    for sample, tag_ids in samples.items():
        sample_count_dir = outs_dir / "per_sample_outs" / sample / "count"
        raw_dir = sample_count_dir / "sample_raw_feature_bc_matrix"
        filtered_dir = sample_count_dir / "sample_filtered_feature_bc_matrix"
        filtered_csv = sample_count_dir / "sample_filtered_barcodes.csv"

        raw_stats = require_mex_dir(raw_dir, errors)
        filtered_stats = require_mex_dir(filtered_dir, errors)
        require_path(filtered_csv, errors)

        filtered_csv_count = sample_filtered_count(filtered_csv) if filtered_csv.exists() else 0
        filtered_mex_count = filtered_stats.get("barcodes", 0)
        if filtered_csv.exists() and filtered_mex_count and filtered_csv_count != filtered_mex_count:
            errors.append(
                f"{sample}: sample_filtered_barcodes.csv count {filtered_csv_count} "
                f"!= filtered MEX barcode count {filtered_mex_count}"
            )

        sample_reports[sample] = {
            "ocm_ids": tag_ids,
            "raw_barcodes": raw_stats.get("barcodes", 0),
            "filtered_barcodes": filtered_mex_count,
            "sample_filtered_barcodes_csv": filtered_csv_count,
        }

    output_cells = load_cells_per_tag(output_cells_path) if output_cells_path.exists() else {}
    oracle_counts = {key: len(oracle_cells.get(key, [])) for key in OCM_IDS}
    output_counts = {key: len(output_cells.get(key, [])) for key in OCM_IDS}
    output_total = sum(output_counts.values())
    if output_total < args.min_total_cells:
        errors.append(f"total OCM cells {output_total} is below minimum {args.min_total_cells}")

    oracle_props = proportions(oracle_counts)
    output_props = proportions(output_counts)
    tag_delta = max_abs_delta(oracle_props, output_props, OCM_IDS)
    if output_total > 0 and tag_delta > args.max_tag_proportion_delta:
        message = (
            f"max OCM tag proportion delta {tag_delta:.3f} exceeds "
            f"{args.max_tag_proportion_delta:.3f}; STAR GeneFull/EmptyDrops_CR "
            "may call additional cells relative to the exonic Cell Ranger oracle"
        )
        if args.strict_tag_proportion_delta:
            errors.append(message)
        else:
            warnings.append(message)

    overlap_by_tag: dict[str, dict[str, float | int]] = {}
    output_overlap_total = 0
    output_cell_total_for_overlap = 0
    oracle_overlap_total = 0
    oracle_cell_total_for_overlap = 0
    for key in OCM_IDS:
        oracle_set = set(oracle_cells.get(key, []))
        output_set = set(output_cells.get(key, []))
        overlap = len(oracle_set & output_set)
        output_denom = len(output_set)
        oracle_denom = len(oracle_set)
        output_frac = overlap / output_denom if output_denom else 0.0
        oracle_recall = overlap / oracle_denom if oracle_denom else 0.0
        output_overlap_total += overlap
        output_cell_total_for_overlap += output_denom
        oracle_overlap_total += overlap
        oracle_cell_total_for_overlap += oracle_denom
        overlap_by_tag[key] = {
            "output": output_denom,
            "oracle": oracle_denom,
            "overlap": overlap,
            "overlap_fraction_of_output": output_frac,
            "oracle_recall": oracle_recall,
            "star_extra": len(output_set - oracle_set),
            "oracle_missed": len(oracle_set - output_set),
        }
        if oracle_denom and oracle_recall < args.min_oracle_recall:
            errors.append(
                f"{key}: oracle recall {oracle_recall:.3f} is below minimum "
                f"{args.min_oracle_recall:.3f}"
            )
    overall_overlap = (
        output_overlap_total / output_cell_total_for_overlap if output_cell_total_for_overlap else 0.0
    )
    overall_oracle_recall = (
        oracle_overlap_total / oracle_cell_total_for_overlap if oracle_cell_total_for_overlap else 0.0
    )
    if args.min_oracle_overlap and overall_overlap < args.min_oracle_overlap:
        errors.append(
            f"oracle barcode overlap {overall_overlap:.3f} is below minimum {args.min_oracle_overlap:.3f}"
        )
    if oracle_cell_total_for_overlap and overall_oracle_recall < args.min_oracle_recall:
        errors.append(
            f"overall oracle recall {overall_oracle_recall:.3f} is below minimum "
            f"{args.min_oracle_recall:.3f}"
        )

    report = {
        "status": "FAIL" if errors else "PASS",
        "star_run_dir": str(star_run_dir),
        "oracle_dir": str(oracle_dir),
        "gene_model_caveat": (
            "Cell Ranger oracle config has include-introns=false, so it is a Gene/exonic "
            "count oracle. This validator checks OCM multi layout and cell assignment, "
            "not GeneFull MEX matrix-value parity."
        ),
        "samples": sample_reports,
        "multi_raw": multi_raw_stats,
        "oracle_counts_per_tag": oracle_counts,
        "output_counts_per_tag": output_counts,
        "oracle_tag_proportions": oracle_props,
        "output_tag_proportions": output_props,
        "max_tag_proportion_delta": tag_delta,
        "strict_tag_proportion_delta": bool(args.strict_tag_proportion_delta),
        "oracle_overlap": {
            "overall_fraction_of_output": overall_overlap,
            "overall_fraction_of_oracle": overall_oracle_recall,
            "by_tag": overlap_by_tag,
        },
        "warnings": warnings,
        "errors": errors,
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"STAR run: {star_run_dir}")
    print(f"Status: {report['status']}")
    print(f"Output OCM cells: {output_counts} total={output_total}")
    print(f"Oracle OCM cells: {oracle_counts} total={sum(oracle_counts.values())}")
    print(f"Max tag proportion delta: {tag_delta:.3f}")
    print(f"Oracle overlap fraction of output: {overall_overlap:.3f}")
    print(f"Oracle recall fraction: {overall_oracle_recall:.3f}")
    print("Gene model caveat: oracle is Gene/exonic; matrix-value parity is not asserted.")
    print(f"Report: {report_path}")
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"  - {warning}")
    if errors:
        print("Errors:")
        for error in errors:
            print(f"  - {error}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
