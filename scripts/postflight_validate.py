#!/usr/bin/env python3
"""Post-flight validation for STAR perturb-seq sample outputs.

Validates a single STAR perturb-seq sample by checking barcode consistency
across surfaces, h5ad integrity and layer presence, CRISPR feature integration,
QC stat sanity, matrix shape agreement, and optional reference comparison.

Usable standalone (one sample) or called in a batch loop via
postflight_validate_batch.py.

Exit codes:
  0  PASS
  1  FAIL (errors found)
  2  FAIL (strict mode, warnings treated as errors)
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

class Severity:
    OK = "OK"
    WARN = "WARN"
    ERROR = "ERROR"
    INFO = "INFO"


@dataclass
class CheckResult:
    section: str
    severity: str
    message: str
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ValidationReport:
    sample_name: str
    checks: List[CheckResult] = field(default_factory=list)

    def add(self, section: str, severity: str, message: str,
            details: Optional[Dict[str, Any]] = None) -> None:
        self.checks.append(CheckResult(
            section=section,
            severity=severity,
            message=message,
            details=details or {},
        ))

    @property
    def has_errors(self) -> bool:
        return any(c.severity == Severity.ERROR for c in self.checks)

    @property
    def has_warnings(self) -> bool:
        return any(c.severity == Severity.WARN for c in self.checks)

    @property
    def verdict(self) -> str:
        if self.has_errors:
            return "FAIL"
        if self.has_warnings:
            return "PASS (with warnings)"
        return "PASS"


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def _open_text(path: Path):
    """Open a plain or gzipped text file."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def _read_barcodes(path: Path) -> List[str]:
    """Read barcodes from a plain or gzipped file, returning raw strings."""
    if not path.exists():
        return []
    with _open_text(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def _strip_suffix(bc: str) -> str:
    """Strip trailing -1 (or -N) suffix common in 10x barcodes."""
    if "-" in bc:
        base, suffix = bc.rsplit("-", 1)
        if suffix.isdigit():
            return base
    return bc


def _barcode_set(path: Path) -> Set[str]:
    """Load barcodes and strip suffix, returning a set."""
    return {_strip_suffix(bc) for bc in _read_barcodes(path)}


def _jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def _count_lines_gz(path: Path) -> int:
    """Count non-empty lines in a (possibly gzipped) file."""
    count = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.strip():
                count += 1
    return count


def _safe_read_h5ad(path: Path):
    """Open h5ad in backed mode, returning (adata, error_string)."""
    import anndata as ad
    try:
        adata = ad.read_h5ad(path, backed="r")
        return adata, None
    except Exception as exc:
        return None, str(exc)


def _close_h5ad(adata) -> None:
    """Safely close a backed AnnData."""
    try:
        adata.file.close()
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_barcode_consistency(report: ValidationReport, run_dir: Path) -> int:
    """Check 1: Surface barcode consistency. Returns GeneFull filtered count."""
    section = "Barcode consistency"

    sources = {
        "GeneFull/filtered": run_dir / "Solo.out/GeneFull/filtered/barcodes.tsv",
        "Velocyto/filtered": run_dir / "Solo.out/Velocyto/filtered/barcodes.tsv",
        "outs/filtered_feature_bc_matrix": run_dir / "outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        "outs/filtered_velocyto_feature_bc_matrix": run_dir / "outs/filtered_velocyto_feature_bc_matrix/barcodes.tsv.gz",
    }

    sets: Dict[str, Set[str]] = {}
    counts: Dict[str, int] = {}
    missing: List[str] = []

    for label, path in sources.items():
        if not path.exists():
            missing.append(label)
            report.add(section, Severity.ERROR,
                       f"Missing barcode file: {label} ({path})")
            continue
        bc_set = _barcode_set(path)
        sets[label] = bc_set
        counts[label] = len(bc_set)

    if missing:
        # Still report what we found
        for label, cnt in counts.items():
            report.add(section, Severity.INFO,
                       f"{label}: {cnt} barcodes", {"count": cnt})
        return counts.get("GeneFull/filtered", 0)

    # Report counts
    for label, cnt in counts.items():
        report.add(section, Severity.INFO,
                   f"{label}: {cnt} barcodes", {"count": cnt})

    # Check all sets match
    reference_label = "GeneFull/filtered"
    reference_set = sets[reference_label]
    all_match = True
    for label, bc_set in sets.items():
        if label == reference_label:
            continue
        if bc_set != reference_set:
            all_match = False
            only_ref = len(reference_set - bc_set)
            only_other = len(bc_set - reference_set)
            report.add(section, Severity.ERROR,
                       f"Mismatch: {reference_label} vs {label} "
                       f"({only_ref} only in ref, {only_other} only in other)",
                       {"ref_only": only_ref, "other_only": only_other})

    if all_match:
        report.add(section, Severity.OK,
                   f"All 4 surfaces match ({counts[reference_label]} barcodes)")

    return counts.get("GeneFull/filtered", 0)


H5AD_NAMES = [
    "counts.h5ad",
    "filtered_counts.h5ad",
    "default_singlet_filtered_counts.h5ad",
    "unfiltered_counts.h5ad",
]


def check_h5ad_integrity(report: ValidationReport, downstream_dir: Path,
                         genefull_count: int) -> Dict[str, Tuple[int, int]]:
    """Check 2: h5ad integrity. Returns dict of {filename: (n_obs, n_vars)}."""
    section = "h5ad integrity"
    shapes: Dict[str, Tuple[int, int]] = {}

    for name in H5AD_NAMES:
        path = downstream_dir / name
        if not path.exists():
            report.add(section, Severity.ERROR,
                       f"Missing h5ad: {name}", {"path": str(path)})
            continue

        adata, err = _safe_read_h5ad(path)
        if err:
            report.add(section, Severity.ERROR,
                       f"Cannot open {name}: {err}")
            continue

        n_obs, n_vars = adata.n_obs, adata.n_vars
        shapes[name] = (n_obs, n_vars)
        report.add(section, Severity.INFO,
                   f"{name}: n_obs={n_obs}, n_vars={n_vars}",
                   {"n_obs": n_obs, "n_vars": n_vars})
        _close_h5ad(adata)

    # counts.h5ad is built from the raw MEX (all non-empty barcodes), so its
    # n_obs will be much larger than the filtered count.  The meaningful check
    # is that filtered_counts.h5ad n_obs <= GeneFull filtered count (QC may
    # remove some cells, but cannot add any).
    if "filtered_counts.h5ad" in shapes and genefull_count > 0:
        fc_obs = shapes["filtered_counts.h5ad"][0]
        if fc_obs <= genefull_count:
            report.add(section, Severity.OK,
                       f"filtered_counts.h5ad n_obs ({fc_obs}) <= GeneFull filtered count ({genefull_count})")
        else:
            report.add(section, Severity.ERROR,
                       f"filtered_counts.h5ad n_obs ({fc_obs}) > GeneFull filtered count ({genefull_count})")

    # filtered_counts.h5ad n_obs <= counts.h5ad n_obs
    if "filtered_counts.h5ad" in shapes and "counts.h5ad" in shapes:
        fc_obs = shapes["filtered_counts.h5ad"][0]
        c_obs = shapes["counts.h5ad"][0]
        if fc_obs <= c_obs:
            report.add(section, Severity.OK,
                       f"filtered_counts.h5ad n_obs ({fc_obs}) <= counts.h5ad n_obs ({c_obs})")
        else:
            report.add(section, Severity.ERROR,
                       f"filtered_counts.h5ad n_obs ({fc_obs}) > counts.h5ad n_obs ({c_obs})")

    # n_vars consistency
    vars_values = [v for _, v in shapes.values()]
    if vars_values:
        if len(set(vars_values)) == 1:
            report.add(section, Severity.OK,
                       f"n_vars consistent across all h5ads ({vars_values[0]})")
        else:
            var_map = {name: v for name, (_, v) in shapes.items()}
            report.add(section, Severity.ERROR,
                       f"n_vars inconsistent across h5ads: {var_map}")

    return shapes


def check_required_layers(report: ValidationReport, downstream_dir: Path,
                          require_layers: List[str]) -> None:
    """Check 3: Required layers present and non-empty in each h5ad."""
    section = "Required layers"

    if not require_layers:
        report.add(section, Severity.OK, "No required layers specified")
        return

    for name in H5AD_NAMES:
        path = downstream_dir / name
        if not path.exists():
            # Already reported in h5ad integrity
            continue

        adata, err = _safe_read_h5ad(path)
        if err:
            continue

        available_layers = list(adata.layers.keys())
        for layer_name in require_layers:
            if layer_name not in available_layers:
                report.add(section, Severity.ERROR,
                           f"{name}: missing required layer '{layer_name}'",
                           {"available_layers": available_layers})
            else:
                # Check nnz > 0
                try:
                    import scipy.sparse as sp
                    import numpy as np
                    layer_data = adata.layers[layer_name]
                    if sp.issparse(layer_data):
                        nnz = layer_data.nnz
                    else:
                        nnz = int(np.count_nonzero(layer_data))
                    if nnz > 0:
                        report.add(section, Severity.OK,
                                   f"{name}: layer '{layer_name}' present (nnz={nnz})",
                                   {"nnz": nnz})
                    else:
                        report.add(section, Severity.ERROR,
                                   f"{name}: layer '{layer_name}' present but empty (nnz=0)")
                except Exception as exc:
                    report.add(section, Severity.WARN,
                               f"{name}: could not check nnz for layer '{layer_name}': {exc}")

        _close_h5ad(adata)


def check_feature_libraries(report: ValidationReport,
                            downstream_dir: Path) -> None:
    """Check 4: Feature library / CRISPR integration."""
    section = "Feature library / CRISPR"

    feat_dir = downstream_dir / "feature_libraries"
    if feat_dir.exists() and feat_dir.is_dir():
        report.add(section, Severity.OK,
                   f"feature_libraries/ directory exists")
    else:
        report.add(section, Severity.WARN,
                   "feature_libraries/ directory missing")

    crispr_obs_cols = [
        "guide_calls", "num_features", "feature_call",
        "guide_identity", "num_guides",
    ]
    crispr_obsm_prefixes = ["crispr", "guide"]

    for name in ("counts.h5ad", "filtered_counts.h5ad"):
        path = downstream_dir / name
        if not path.exists():
            continue

        adata, err = _safe_read_h5ad(path)
        if err:
            continue

        found_obs = [c for c in adata.obs.columns
                     if any(c.lower().startswith(p) or p in c.lower()
                            for p in crispr_obs_cols)]
        found_obsm = [k for k in adata.obsm.keys()
                      if any(p in k.lower() for p in crispr_obsm_prefixes)]

        if found_obs or found_obsm:
            report.add(section, Severity.OK,
                       f"{name}: CRISPR annotations found",
                       {"obs_columns": found_obs, "obsm_keys": found_obsm})
        else:
            report.add(section, Severity.WARN,
                       f"{name}: no CRISPR annotations found in obs/obsm",
                       {"obs_columns_checked": crispr_obs_cols,
                        "obsm_prefixes_checked": crispr_obsm_prefixes})

        _close_h5ad(adata)


def check_qc_stats(report: ValidationReport, downstream_dir: Path,
                    expected_cells_range: Optional[Tuple[int, int]],
                    h5ad_shapes: Dict[str, Tuple[int, int]]) -> None:
    """Check 5: QC stats sanity."""
    section = "QC stats"

    # Read adaptive_qc_threshold.json
    aqc_path = downstream_dir / "adaptive_qc_threshold.json"
    if aqc_path.exists():
        try:
            with open(aqc_path, "r") as fh:
                aqc = json.load(fh)
            effective_max = aqc.get("effective_max_genes")
            report.add(section, Severity.INFO,
                       f"adaptive_qc_threshold.json: effective_max_genes={effective_max}",
                       {"adaptive_qc": aqc})
        except Exception as exc:
            report.add(section, Severity.WARN,
                       f"Cannot parse adaptive_qc_threshold.json: {exc}")
    else:
        report.add(section, Severity.INFO,
                   "adaptive_qc_threshold.json not found (may not be adaptive)")

    # Expected cell range check using filtered_counts.h5ad
    if expected_cells_range is not None:
        low, high = expected_cells_range
        if "filtered_counts.h5ad" in h5ad_shapes:
            n_cells = h5ad_shapes["filtered_counts.h5ad"][0]
            if low <= n_cells <= high:
                report.add(section, Severity.OK,
                           f"Filtered cell count ({n_cells}) within expected range [{low}, {high}]")
            else:
                report.add(section, Severity.WARN,
                           f"Filtered cell count ({n_cells}) outside expected range [{low}, {high}]",
                           {"n_cells": n_cells, "expected_low": low, "expected_high": high})
        else:
            report.add(section, Severity.WARN,
                       "Cannot check cell range: filtered_counts.h5ad not available")

    # Read summary.txt
    summary_path = downstream_dir / "summary.txt"
    if summary_path.exists():
        try:
            text = summary_path.read_text()
            report.add(section, Severity.INFO,
                       f"summary.txt content:\n{text.rstrip()}",
                       {"summary_text": text})
        except Exception as exc:
            report.add(section, Severity.WARN,
                       f"Cannot read summary.txt: {exc}")
    else:
        report.add(section, Severity.WARN, "summary.txt not found")


def check_matrix_shape(report: ValidationReport, run_dir: Path,
                       h5ad_shapes: Dict[str, Tuple[int, int]]) -> None:
    """Check 6: Matrix shape agreement — n_vars vs features.tsv line count.

    The h5ads are built from GeneFull counts, so compare against
    Solo.out/GeneFull/filtered/features.tsv (not the outs/ features file which
    may include additional CRISPR/feature library rows).
    """
    section = "Matrix shape agreement"

    # Prefer GeneFull features file; fall back to outs/ if missing
    genefull_features = run_dir / "Solo.out/GeneFull/filtered/features.tsv"
    outs_features = run_dir / "outs/filtered_feature_bc_matrix/features.tsv.gz"
    features_path = genefull_features if genefull_features.exists() else outs_features
    if not features_path.exists():
        report.add(section, Severity.ERROR,
                   f"Missing features file: {features_path}")
        return

    try:
        n_features = _count_lines_gz(features_path)
    except Exception as exc:
        report.add(section, Severity.ERROR,
                   f"Cannot read features file: {exc}")
        return

    feat_label = features_path.name
    report.add(section, Severity.INFO,
               f"{feat_label} line count: {n_features}")

    mismatch = False
    for name, (_, n_vars) in h5ad_shapes.items():
        if n_vars != n_features:
            report.add(section, Severity.ERROR,
                       f"{name} n_vars ({n_vars}) != {feat_label} ({n_features})")
            mismatch = True

    if not mismatch and h5ad_shapes:
        report.add(section, Severity.OK,
                   f"All h5ad n_vars match {feat_label} ({n_features})")


def check_reference_comparison(report: ValidationReport, run_dir: Path,
                               downstream_dir: Path,
                               ref_downstream_dir: Optional[Path],
                               ref_filtered_barcodes: Optional[Path]) -> None:
    """Check 7: Reference comparison (Jaccard, stat diff)."""
    section = "Reference comparison"

    # Determine new filtered barcodes
    new_bc_path = run_dir / "Solo.out/GeneFull/filtered/barcodes.tsv"
    new_bcs = _barcode_set(new_bc_path)
    if not new_bcs:
        report.add(section, Severity.WARN,
                   "Cannot load new filtered barcodes for comparison")
        return

    # Get old barcodes
    old_bcs: Set[str] = set()
    if ref_filtered_barcodes:
        old_bcs = _barcode_set(ref_filtered_barcodes)
    elif ref_downstream_dir:
        # Try to find filtered barcodes in reference downstream
        # The canonical location depends on what the reference was
        # Try counts.h5ad obs_names
        ref_counts = ref_downstream_dir / "counts.h5ad"
        if ref_counts.exists():
            adata, err = _safe_read_h5ad(ref_counts)
            if adata is not None:
                old_bcs = {_strip_suffix(str(bc)) for bc in adata.obs_names}
                _close_h5ad(adata)

    if old_bcs:
        intersection = len(new_bcs & old_bcs)
        j = _jaccard(new_bcs, old_bcs)
        report.add(section, Severity.INFO,
                   f"Barcode Jaccard: {j:.6f} "
                   f"(old={len(old_bcs)}, new={len(new_bcs)}, "
                   f"intersection={intersection}, "
                   f"old_only={len(old_bcs - new_bcs)}, "
                   f"new_only={len(new_bcs - old_bcs)})",
                   {"jaccard": j,
                    "old_count": len(old_bcs),
                    "new_count": len(new_bcs),
                    "intersection": intersection,
                    "old_only": len(old_bcs - new_bcs),
                    "new_only": len(new_bcs - old_bcs)})
    else:
        report.add(section, Severity.WARN,
                   "Could not load reference barcodes for Jaccard comparison")

    # Compare summary.txt side by side
    if ref_downstream_dir:
        old_summary = ref_downstream_dir / "summary.txt"
        new_summary = downstream_dir / "summary.txt"
        if old_summary.exists() and new_summary.exists():
            try:
                old_text = old_summary.read_text().strip()
                new_text = new_summary.read_text().strip()
                old_lines = old_text.splitlines()
                new_lines = new_text.splitlines()
                comparison_lines = []
                max_lines = max(len(old_lines), len(new_lines))
                for i in range(max_lines):
                    old_l = old_lines[i] if i < len(old_lines) else ""
                    new_l = new_lines[i] if i < len(new_lines) else ""
                    marker = " " if old_l == new_l else "*"
                    comparison_lines.append(
                        f"  {marker} OLD: {old_l}\n  {marker} NEW: {new_l}")
                comparison = "\n".join(comparison_lines)
                report.add(section, Severity.INFO,
                           f"summary.txt comparison:\n{comparison}",
                           {"old_summary": old_text, "new_summary": new_text})
            except Exception as exc:
                report.add(section, Severity.WARN,
                           f"Could not compare summary.txt files: {exc}")
        else:
            if not old_summary.exists():
                report.add(section, Severity.WARN,
                           "Reference summary.txt not found")


# ---------------------------------------------------------------------------
# Main validation orchestrator
# ---------------------------------------------------------------------------

def validate_sample(run_dir: Path, downstream_dir: Path, *,
                    ref_downstream_dir: Optional[Path] = None,
                    ref_filtered_barcodes: Optional[Path] = None,
                    expected_cells_range: Optional[Tuple[int, int]] = None,
                    require_layers: Optional[List[str]] = None,
                    sample_name: Optional[str] = None) -> ValidationReport:
    """Run all validation checks on a single sample."""
    if sample_name is None:
        sample_name = run_dir.parent.name

    if require_layers is None:
        require_layers = ["denoised"]

    report = ValidationReport(sample_name=sample_name)

    # 1. Barcode consistency
    genefull_count = check_barcode_consistency(report, run_dir)

    # 2. h5ad integrity
    h5ad_shapes = check_h5ad_integrity(report, downstream_dir, genefull_count)

    # 3. Required layers
    check_required_layers(report, downstream_dir, require_layers)

    # 4. Feature library / CRISPR
    check_feature_libraries(report, downstream_dir)

    # 5. QC stats
    check_qc_stats(report, downstream_dir, expected_cells_range, h5ad_shapes)

    # 6. Matrix shape agreement
    check_matrix_shape(report, run_dir, h5ad_shapes)

    # 7. Reference comparison (only if requested)
    if ref_downstream_dir or ref_filtered_barcodes:
        check_reference_comparison(report, run_dir, downstream_dir,
                                   ref_downstream_dir, ref_filtered_barcodes)

    return report


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

def format_report(report: ValidationReport) -> str:
    """Format a ValidationReport as structured text."""
    lines: List[str] = []
    lines.append("=" * 80)
    lines.append(f"POST-FLIGHT VALIDATION: {report.sample_name}")
    lines.append("=" * 80)

    current_section = None
    for check in report.checks:
        if check.section != current_section:
            current_section = check.section
            lines.append("")
            lines.append(f"--- {current_section} ---")

        tag = f"[{check.severity}]"
        # Indent multi-line messages
        msg_lines = check.message.splitlines()
        lines.append(f"  {tag:8s} {msg_lines[0]}")
        for ml in msg_lines[1:]:
            lines.append(f"           {ml}")

    lines.append("")
    lines.append("-" * 80)
    lines.append(f"VERDICT: {report.verdict}")
    lines.append("-" * 80)

    return "\n".join(lines)


def report_to_dict(report: ValidationReport) -> Dict[str, Any]:
    """Convert a ValidationReport to a JSON-serializable dict."""
    return {
        "sample_name": report.sample_name,
        "verdict": report.verdict,
        "has_errors": report.has_errors,
        "has_warnings": report.has_warnings,
        "checks": [asdict(c) for c in report.checks],
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # Required
    parser.add_argument("--run-dir", required=True, type=Path,
                        help="STAR run directory containing Solo.out/ and outs/")
    parser.add_argument("--downstream-dir", required=True, type=Path,
                        help="Downstream output directory containing h5ads, summary.txt, etc.")

    # Optional
    parser.add_argument("--reference-downstream-dir", type=Path, default=None,
                        help="Old downstream dir for comparison (Jaccard, stat diff)")
    parser.add_argument("--reference-filtered-barcodes", type=Path, default=None,
                        help="Old filtered barcodes.tsv for Jaccard comparison")
    parser.add_argument("--expected-cells-range", type=parse_cells_range, default=None,
                        metavar="LOW-HIGH",
                        help="Warn if filtered cell count outside range (e.g. 5000-25000)")
    parser.add_argument("--require-layers", nargs="+", default=["denoised"],
                        metavar="LAYER",
                        help="Layer names that must exist in h5ads (default: denoised)")
    parser.add_argument("--sample-name", default=None,
                        help="Label for output (default: basename of run-dir parent)")
    parser.add_argument("--json-output", type=Path, default=None,
                        help="Write detailed results as JSON")
    parser.add_argument("--strict", action="store_true",
                        help="Exit non-zero on any WARNING (default: only errors)")

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    downstream_dir = args.downstream_dir.resolve()
    ref_ds = args.reference_downstream_dir.resolve() if args.reference_downstream_dir else None
    ref_bc = args.reference_filtered_barcodes.resolve() if args.reference_filtered_barcodes else None

    report = validate_sample(
        run_dir=run_dir,
        downstream_dir=downstream_dir,
        ref_downstream_dir=ref_ds,
        ref_filtered_barcodes=ref_bc,
        expected_cells_range=args.expected_cells_range,
        require_layers=args.require_layers,
        sample_name=args.sample_name,
    )

    print(format_report(report))

    if args.json_output:
        out_path = args.json_output.resolve()
        with open(out_path, "w") as fh:
            json.dump(report_to_dict(report), fh, indent=2)
        print(f"\nJSON report written to {out_path}", file=sys.stderr)

    # Exit code
    if report.has_errors:
        return 1
    if args.strict and report.has_warnings:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
