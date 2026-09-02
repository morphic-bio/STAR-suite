#!/usr/bin/env python3
"""Compare Flex hash-screen and align-everything MEX outputs by biological keys."""

import argparse
from pathlib import Path


def read_lines(path):
    with open(path, "r", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


def read_mex(directory):
    directory = Path(directory)
    features = read_lines(directory / "features.tsv")
    barcodes = read_lines(directory / "barcodes.tsv")
    feature_keys = [line.split("\t", 1)[0] for line in features]
    if len(feature_keys) != len(set(feature_keys)):
        raise ValueError(f"duplicate feature IDs in {directory / 'features.tsv'}")
    if len(barcodes) != len(set(barcodes)):
        raise ValueError(f"duplicate barcodes in {directory / 'barcodes.tsv'}")

    counts = {}
    dimensions = None
    with open(directory / "matrix.mtx", "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if dimensions is None:
                dimensions = tuple(int(value) for value in fields)
                if len(dimensions) != 3:
                    raise ValueError(f"invalid MatrixMarket dimensions in {directory}")
                continue
            feature_index, barcode_index, count = (int(value) for value in fields)
            key = (feature_keys[feature_index - 1], barcodes[barcode_index - 1])
            if key in counts:
                raise ValueError(f"duplicate MatrixMarket coordinate in {directory}: {key}")
            counts[key] = count

    if dimensions is None:
        raise ValueError(f"missing MatrixMarket dimensions in {directory}")
    if dimensions[0] != len(features) or dimensions[1] != len(barcodes):
        raise ValueError(f"MatrixMarket dimensions do not match labels in {directory}")
    if dimensions[2] != len(counts):
        raise ValueError(f"MatrixMarket nnz does not match data rows in {directory}")
    return features, barcodes, counts


def fraction(numerator, denominator):
    return numerator / max(denominator, 1)


def compare_one(label, baseline_dir, current_dir, limits):
    baseline_features, baseline_barcodes, baseline = read_mex(baseline_dir)
    current_features, current_barcodes, current = read_mex(current_dir)

    baseline_feature_set = set(baseline_features)
    current_feature_set = set(current_features)
    if baseline_feature_set != current_feature_set:
        raise AssertionError(f"{label}: feature sets differ")

    baseline_barcode_set = set(baseline_barcodes)
    current_barcode_set = set(current_barcodes)
    barcode_symmetric_difference = len(baseline_barcode_set ^ current_barcode_set)
    barcode_difference_fraction = fraction(
        barcode_symmetric_difference,
        max(len(baseline_barcode_set), len(current_barcode_set)),
    )

    all_keys = baseline.keys() | current.keys()
    mismatches = [key for key in all_keys if baseline.get(key, 0) != current.get(key, 0)]
    mismatch_fraction = fraction(len(mismatches), max(len(baseline), len(current)))
    baseline_total = sum(baseline.values())
    current_total = sum(current.values())
    count_delta_fraction = fraction(abs(current_total - baseline_total), baseline_total)
    maximum_coordinate_delta = max(
        (abs(current.get(key, 0) - baseline.get(key, 0)) for key in mismatches),
        default=0,
    )
    only_baseline = sum(key not in current for key in baseline)
    only_current = sum(key not in baseline for key in current)

    print(
        f"{label}\tbaseline_nnz={len(baseline)}\tcurrent_nnz={len(current)}"
        f"\tmismatches={len(mismatches)}\tmismatch_fraction={mismatch_fraction:.9f}"
        f"\tbaseline_total={baseline_total}\tcurrent_total={current_total}"
        f"\tcount_delta_fraction={count_delta_fraction:.9f}"
        f"\tbarcode_symmetric_difference={barcode_symmetric_difference}"
        f"\tbarcode_difference_fraction={barcode_difference_fraction:.9f}"
        f"\tonly_baseline={only_baseline}\tonly_current={only_current}"
        f"\tmax_coordinate_delta={maximum_coordinate_delta}"
    )

    failures = []
    if mismatch_fraction > limits.max_mismatch_fraction:
        failures.append("matrix mismatch fraction")
    if count_delta_fraction > limits.max_count_delta_fraction:
        failures.append("total-count delta fraction")
    if barcode_difference_fraction > limits.max_barcode_difference_fraction:
        failures.append("barcode-set difference fraction")
    if maximum_coordinate_delta > limits.max_coordinate_delta:
        failures.append("maximum coordinate delta")
    if failures:
        raise AssertionError(f"{label}: exceeded {', '.join(failures)}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("baseline")
    parser.add_argument("current")
    parser.add_argument("--samples", default="BC004,BC006,BC007,BC008")
    parser.add_argument("--max-mismatch-fraction", type=float, default=0.002)
    parser.add_argument("--max-count-delta-fraction", type=float, default=0.002)
    parser.add_argument("--max-barcode-difference-fraction", type=float, default=0.001)
    parser.add_argument("--max-coordinate-delta", type=int, default=1)
    args = parser.parse_args()

    compare_one(
        "pooled_raw",
        Path(args.baseline) / "Solo.out/Gene/raw",
        Path(args.current) / "Solo.out/Gene/raw",
        args,
    )
    for sample in (value for value in args.samples.split(",") if value):
        compare_one(
            sample,
            Path(args.baseline) / f"per_sample/{sample}/Gene/filtered",
            Path(args.current) / f"per_sample/{sample}/Gene/filtered",
            args,
        )
    print("Flex hash-screen keyed MEX parity passed")


if __name__ == "__main__":
    main()
