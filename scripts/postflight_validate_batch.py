#!/usr/bin/env python3
"""Batch post-flight validation for STAR perturb-seq samples.

Auto-discovers samples under prod-root/samples/*/ and runs
postflight_validate.py checks on each, printing a summary table.

Exit codes:
  0  All samples PASS
  1  At least one sample has ERRORs
  2  Strict mode: at least one sample has WARNINGs
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional, Tuple

from postflight_validate import (
    ValidationReport,
    format_report,
    report_to_dict,
    validate_sample,
)


DEFAULT_DOWNSTREAM_SUFFIX = "downstream_genefull_velocyto_cellbender"
DEFAULT_RUN_SUBDIR = "run"


def discover_samples(prod_root: Path) -> List[str]:
    """Find sample names under prod-root/samples/."""
    samples_dir = prod_root / "samples"
    if not samples_dir.is_dir():
        print(f"ERROR: samples directory not found: {samples_dir}",
              file=sys.stderr)
        return []
    names = sorted(
        d.name for d in samples_dir.iterdir()
        if d.is_dir() and not d.name.startswith(".")
        and (d / "run").is_dir()
    )
    return names


def print_summary_table(reports: List[ValidationReport]) -> None:
    """Print a compact summary table across all samples."""
    print()
    print("=" * 100)
    print("BATCH VALIDATION SUMMARY")
    print("=" * 100)

    # Collect per-sample stats
    header = (f"{'Sample':<25} {'Verdict':<22} {'Errors':>6} {'Warnings':>8} "
              f"{'OK':>4} {'Info':>5}")
    print(header)
    print("-" * 100)

    for report in reports:
        n_err = sum(1 for c in report.checks if c.severity == "ERROR")
        n_warn = sum(1 for c in report.checks if c.severity == "WARN")
        n_ok = sum(1 for c in report.checks if c.severity == "OK")
        n_info = sum(1 for c in report.checks if c.severity == "INFO")
        print(f"{report.sample_name:<25} {report.verdict:<22} {n_err:>6} "
              f"{n_warn:>8} {n_ok:>4} {n_info:>5}")

    print("-" * 100)

    any_errors = any(r.has_errors for r in reports)
    any_warnings = any(r.has_warnings for r in reports)
    if any_errors:
        print("OVERALL: FAIL (errors found)")
    elif any_warnings:
        print("OVERALL: PASS (with warnings)")
    else:
        print("OVERALL: PASS")
    print("=" * 100)


def parse_cells_range(value: str) -> Tuple[int, int]:
    """Parse LOW-HIGH string into (low, high) tuple."""
    parts = value.split("-")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"Expected format LOW-HIGH (e.g. 5000-25000), got: {value}")
    try:
        low, high = int(parts[0]), int(parts[1])
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Expected integers in LOW-HIGH, got: {value}")
    if low > high:
        raise argparse.ArgumentTypeError(
            f"LOW ({low}) must be <= HIGH ({high})")
    return low, high


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--prod-root", required=True, type=Path,
                        help="Production root directory containing samples/")
    parser.add_argument("--downstream-suffix", default=DEFAULT_DOWNSTREAM_SUFFIX,
                        help=f"Downstream directory name (default: {DEFAULT_DOWNSTREAM_SUFFIX})")
    parser.add_argument("--run-subdir", default=DEFAULT_RUN_SUBDIR,
                        help=f"Run subdirectory name within each sample (default: {DEFAULT_RUN_SUBDIR})")
    parser.add_argument("--samples", nargs="*", default=None,
                        help="Specific sample names (default: auto-discover)")
    parser.add_argument("--reference-downstream-suffix", default=None,
                        help="Old downstream suffix for comparison (e.g. downstream_old)")
    parser.add_argument("--expected-cells-range", type=parse_cells_range, default=None,
                        metavar="LOW-HIGH",
                        help="Warn if filtered cell count outside range")
    parser.add_argument("--require-layers", nargs="+", default=["denoised"],
                        metavar="LAYER",
                        help="Layer names required in h5ads (default: denoised)")
    parser.add_argument("--json-output", type=Path, default=None,
                        help="Write batch JSON report")
    parser.add_argument("--strict", action="store_true",
                        help="Exit non-zero on any WARNING")
    parser.add_argument("--verbose", action="store_true",
                        help="Print per-sample detail reports (default: summary only)")

    args = parser.parse_args()
    prod_root = args.prod_root.resolve()

    # Discover or use provided samples
    if args.samples:
        sample_names = args.samples
    else:
        sample_names = discover_samples(prod_root)

    if not sample_names:
        print("No samples found.", file=sys.stderr)
        return 1

    print(f"Validating {len(sample_names)} samples under {prod_root}",
          file=sys.stderr)

    reports: List[ValidationReport] = []

    for sample_name in sample_names:
        sample_dir = prod_root / "samples" / sample_name
        run_dir = sample_dir / args.run_subdir
        downstream_dir = sample_dir / args.downstream_suffix

        # Reference downstream
        ref_ds: Optional[Path] = None
        if args.reference_downstream_suffix:
            ref_ds = sample_dir / args.reference_downstream_suffix
            if not ref_ds.exists():
                ref_ds = None

        print(f"  Validating {sample_name}...", file=sys.stderr)

        if not run_dir.exists():
            report = ValidationReport(sample_name=sample_name)
            report.add("Setup", "ERROR",
                       f"Run directory not found: {run_dir}")
            reports.append(report)
            continue

        if not downstream_dir.exists():
            report = ValidationReport(sample_name=sample_name)
            report.add("Setup", "ERROR",
                       f"Downstream directory not found: {downstream_dir}")
            reports.append(report)
            continue

        report = validate_sample(
            run_dir=run_dir,
            downstream_dir=downstream_dir,
            ref_downstream_dir=ref_ds,
            expected_cells_range=args.expected_cells_range,
            require_layers=args.require_layers,
            sample_name=sample_name,
        )
        reports.append(report)

        if args.verbose:
            print(format_report(report))

    # Summary
    print_summary_table(reports)

    # JSON output
    if args.json_output:
        out_path = args.json_output.resolve()
        batch_result = {
            "prod_root": str(prod_root),
            "downstream_suffix": args.downstream_suffix,
            "n_samples": len(reports),
            "samples": [report_to_dict(r) for r in reports],
        }
        with open(out_path, "w") as fh:
            json.dump(batch_result, fh, indent=2)
        print(f"\nBatch JSON report written to {out_path}", file=sys.stderr)

    # Exit code
    any_errors = any(r.has_errors for r in reports)
    any_warnings = any(r.has_warnings for r in reports)
    if any_errors:
        return 1
    if args.strict and any_warnings:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
