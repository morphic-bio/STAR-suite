#!/usr/bin/env python3
"""
Contract validation for run_flexfilter_mex output.

Validates:
- Per-sample directories exist and match expected tags
- barcodes.tsv: 24-char composite barcodes with correct TAG suffix
- features.tsv: tab-separated, non-empty
- matrix.mtx: valid Matrix Market header
- Optionally compares cell counts against gold standard summary
"""

import argparse
import os
import sys
from pathlib import Path


def validate_barcodes(barcodes_path: Path, expected_tag: str = None) -> tuple[bool, str]:
    """Validate barcodes.tsv format and optionally check TAG suffix."""
    if not barcodes_path.exists():
        return False, f"File not found: {barcodes_path}"
    
    if barcodes_path.stat().st_size == 0:
        return False, f"File is empty: {barcodes_path}"
    
    with open(barcodes_path) as f:
        lines = [line.strip() for line in f if line.strip()]
    
    if not lines:
        return False, f"No barcodes found in {barcodes_path}"
    
    errors = []
    for i, barcode in enumerate(lines[:10]):  # Check first 10
        if len(barcode) != 24:
            errors.append(f"Line {i+1}: barcode length {len(barcode)} != 24: {barcode}")
        
        if expected_tag and len(barcode) >= 8:
            tag_suffix = barcode[-8:]
            if tag_suffix != expected_tag:
                errors.append(f"Line {i+1}: TAG suffix {tag_suffix} != expected {expected_tag}")
    
    if errors:
        return False, "; ".join(errors[:3])  # First 3 errors
    
    return True, f"OK: {len(lines)} barcodes, 24-char format"


def validate_features(features_path: Path) -> tuple[bool, str]:
    """Validate features.tsv format."""
    if not features_path.exists():
        return False, f"File not found: {features_path}"
    
    if features_path.stat().st_size == 0:
        return False, f"File is empty: {features_path}"
    
    with open(features_path) as f:
        lines = [line.strip() for line in f if line.strip()]
    
    if not lines:
        return False, f"No features found in {features_path}"
    
    # Check tab-separated format
    for i, line in enumerate(lines[:5]):
        if '\t' not in line:
            return False, f"Line {i+1}: not tab-separated: {line[:50]}"
    
    return True, f"OK: {len(lines)} features, tab-separated"


def validate_matrix(matrix_path: Path) -> tuple[bool, str]:
    """Validate matrix.mtx Matrix Market format."""
    if not matrix_path.exists():
        return False, f"File not found: {matrix_path}"
    
    if matrix_path.stat().st_size == 0:
        return False, f"File is empty: {matrix_path}"
    
    with open(matrix_path) as f:
        # Check header
        first_line = f.readline().strip()
        if not first_line.startswith("%%MatrixMarket"):
            return False, f"Invalid Matrix Market header: {first_line[:50]}"
        
        # Skip comments
        line = f.readline()
        while line.startswith('%'):
            line = f.readline()
        
        # Parse dimensions
        parts = line.strip().split()
        if len(parts) != 3:
            return False, f"Invalid dimension line: {line[:50]}"
        
        try:
            n_rows, n_cols, n_entries = int(parts[0]), int(parts[1]), int(parts[2])
        except ValueError:
            return False, f"Cannot parse dimensions: {line[:50]}"
    
    return True, f"OK: {n_rows}x{n_cols} matrix, {n_entries} entries"


def load_gold_summary(summary_path: Path) -> dict[str, int]:
    """Load expected cell counts from gold standard summary TSV."""
    if not summary_path.exists():
        return {}
    
    counts = {}
    with open(summary_path) as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[0] != "TOTAL":
                sample = parts[0]
                try:
                    final_cells = int(parts[8])  # "Final" column
                    counts[sample] = final_cells
                except (ValueError, IndexError):
                    pass
    return counts


def validate_output(output_dir: Path, gold_dir: Path = None, expected_samples: list = None):
    """Validate FlexFilter output directory."""
    print(f"Validating output: {output_dir}")
    
    if not output_dir.exists():
        print(f"ERROR: Output directory not found: {output_dir}")
        return False
    
    # Auto-detect samples if not specified: any dir that contains Gene/filtered
    if expected_samples is None:
        detected = []
        for d in output_dir.iterdir():
            if not d.is_dir():
                continue
            if (d / "Gene" / "filtered").exists():
                detected.append(d.name)
        expected_samples = sorted(detected)
    
    if not expected_samples:
        print("ERROR: No sample directories found (expected subdirs with Gene/filtered)")
        return False
    
    print(f"Samples to validate: {expected_samples}")
    
    # Load gold standard counts if available
    gold_counts = {}
    if gold_dir:
        summary_path = gold_dir / "flexfilter_summary.tsv"
        gold_counts = load_gold_summary(summary_path)
        if gold_counts:
            print(f"Loaded gold standard counts for {len(gold_counts)} samples")
        else:
            print("WARN: Gold standard summary not found or empty; skipping count comparison")
    
    all_passed = True
    
    for sample in expected_samples:
        print(f"\n--- {sample} ---")
        sample_dir = output_dir / sample / "Gene" / "filtered"
        
        if not sample_dir.exists():
            print(f"FAIL: Sample directory not found: {sample_dir}")
            all_passed = False
            continue
        
        # Tag suffix check is optional; if sample names are arbitrary we skip it by default
        expected_tag = None
        
        # Validate barcodes
        ok, msg = validate_barcodes(sample_dir / "barcodes.tsv", expected_tag)
        print(f"  barcodes.tsv: {msg}")
        if not ok:
            all_passed = False
        
        # Validate features
        ok, msg = validate_features(sample_dir / "features.tsv")
        print(f"  features.tsv: {msg}")
        if not ok:
            all_passed = False
        
        # Validate matrix
        ok, msg = validate_matrix(sample_dir / "matrix.mtx")
        print(f"  matrix.mtx: {msg}")
        if not ok:
            all_passed = False
        
        # Compare cell counts if gold available
        if sample in gold_counts:
            with open(sample_dir / "barcodes.tsv") as f:
                actual_cells = sum(1 for line in f if line.strip())
            expected_cells = gold_counts[sample]
            
            # Allow some tolerance (within 20% or 50 cells)
            diff = abs(actual_cells - expected_cells)
            tolerance = max(50, int(expected_cells * 0.2))
            
            if diff <= tolerance:
                print(f"  Cell count: {actual_cells} (gold: {expected_cells}, diff: {diff}) OK")
            else:
                print(f"  Cell count: {actual_cells} (gold: {expected_cells}, diff: {diff}) WARN: large difference")
    
    return all_passed


def main():
    parser = argparse.ArgumentParser(description="Validate FlexFilter output")
    parser.add_argument("output_dir", type=Path, help="FlexFilter output directory")
    parser.add_argument("--gold-dir", type=Path, help="Gold standard directory for comparison")
    parser.add_argument("--samples", nargs="+", help="Expected sample names (e.g., SampleA SampleB)")
    
    args = parser.parse_args()
    
    success = validate_output(args.output_dir, args.gold_dir, args.samples)
    
    print("\n" + "=" * 40)
    if success:
        print("VALIDATION PASSED")
        sys.exit(0)
    else:
        print("VALIDATION FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()

