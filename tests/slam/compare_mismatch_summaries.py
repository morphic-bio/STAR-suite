#!/usr/bin/env python3
"""
Compare STAR-Slam mismatch summaries against GRAND-SLAM reference outputs.

Supports:
  - mismatches.tsv (Category, Condition, Orientation, Genomic, Read, Coverage, Mismatches)
  - mismatchdetails.tsv (Category, Genomic, Read, Position, Overlap, Opposite, Coverage, Mismatches)
"""

import argparse
import csv
import gzip
import sys
from collections import defaultdict


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_mismatches(path, match_condition):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = ["Category", "Condition", "Orientation", "Genomic", "Read", "Coverage", "Mismatches"]
        for col in required:
            if col not in reader.fieldnames:
                raise ValueError(f"{path} missing column: {col}")
        data = defaultdict(lambda: [0.0, 0.0])
        for row in reader:
            if not row:
                continue
            if match_condition:
                key = (row["Category"], row["Condition"], row["Orientation"], row["Genomic"], row["Read"])
            else:
                key = (row["Category"], row["Orientation"], row["Genomic"], row["Read"])
            cov = float(row["Coverage"])
            mis = float(row["Mismatches"])
            data[key][0] += cov
            data[key][1] += mis
        return data


def parse_mismatch_details(path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = ["Category", "Genomic", "Read", "Position", "Overlap", "Opposite", "Coverage", "Mismatches"]
        for col in required:
            if col not in reader.fieldnames:
                raise ValueError(f"{path} missing column: {col}")
        data = defaultdict(lambda: [0.0, 0.0])
        for row in reader:
            if not row:
                continue
            key = (
                row["Category"],
                row["Genomic"],
                row["Read"],
                int(row["Position"]),
                int(row["Overlap"]),
                int(row["Opposite"]),
            )
            cov = float(row["Coverage"])
            mis = float(row["Mismatches"])
            data[key][0] += cov
            data[key][1] += mis
        return data


def compare_tables(ref, test, tol, label, max_report):
    keys = set(ref) | set(test)
    missing_in_test = [k for k in keys if k not in test]
    missing_in_ref = [k for k in keys if k not in ref]
    mismatches = []
    for key in keys:
        r_cov, r_mis = ref.get(key, (0.0, 0.0))
        t_cov, t_mis = test.get(key, (0.0, 0.0))
        dcov = abs(r_cov - t_cov)
        dmis = abs(r_mis - t_mis)
        if dcov > tol or dmis > tol:
            mismatches.append((max(dcov, dmis), key, r_cov, t_cov, r_mis, t_mis))
    mismatches.sort(reverse=True)

    ok = True
    if missing_in_test:
        ok = False
        print(f"FAIL: {label} missing rows in test: {len(missing_in_test)}")
        for key in missing_in_test[:max_report]:
            print(f"  missing_in_test\t{key}")
    if missing_in_ref:
        ok = False
        print(f"FAIL: {label} missing rows in reference: {len(missing_in_ref)}")
        for key in missing_in_ref[:max_report]:
            print(f"  missing_in_ref\t{key}")
    if mismatches:
        ok = False
        print(f"FAIL: {label} mismatches: {len(mismatches)}")
        for delta, key, r_cov, t_cov, r_mis, t_mis in mismatches[:max_report]:
            print(f"  {key}\tref_cov={r_cov}\ttest_cov={t_cov}\tref_mis={r_mis}\ttest_mis={t_mis}\tdelta={delta}")
    if ok:
        print(f"PASS: {label} matched within tol={tol}")
    return ok


def main():
    parser = argparse.ArgumentParser(description="Compare SLAM mismatch summary outputs.")
    parser.add_argument("--reference-mismatches", help="Reference mismatches.tsv (GEDI)")
    parser.add_argument("--test-mismatches", help="Test mismatches.tsv (STAR)")
    parser.add_argument("--reference-details", help="Reference mismatchdetails.tsv (GEDI)")
    parser.add_argument("--test-details", help="Test mismatchdetails.tsv (STAR)")
    parser.add_argument("--match-condition", action="store_true",
                        help="Include Condition column in mismatches.tsv key")
    parser.add_argument("--tol", type=float, default=0.0,
                        help="Absolute tolerance for coverage/mismatch deltas")
    parser.add_argument("--max-report", type=int, default=10,
                        help="Maximum rows to report per section")
    args = parser.parse_args()

    ok = True
    if args.reference_mismatches and args.test_mismatches:
        ref = parse_mismatches(args.reference_mismatches, args.match_condition)
        test = parse_mismatches(args.test_mismatches, args.match_condition)
        ok &= compare_tables(ref, test, args.tol, "mismatches.tsv", args.max_report)
    elif args.reference_mismatches or args.test_mismatches:
        raise ValueError("Both --reference-mismatches and --test-mismatches are required")

    if args.reference_details and args.test_details:
        ref = parse_mismatch_details(args.reference_details)
        test = parse_mismatch_details(args.test_details)
        ok &= compare_tables(ref, test, args.tol, "mismatchdetails.tsv", args.max_report)
    elif args.reference_details or args.test_details:
        raise ValueError("Both --reference-details and --test-details are required")

    if not (args.reference_mismatches or args.reference_details):
        raise ValueError("No inputs provided for comparison")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
