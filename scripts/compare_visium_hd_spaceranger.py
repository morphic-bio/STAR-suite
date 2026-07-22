#!/usr/bin/env python3
"""Descriptively compare Visium HD policy MEX outputs with Space Ranger raw H5."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from pathlib import Path

import h5py
import numpy as np
from scipy.io import mmread
from scipy.stats import pearsonr, spearmanr


POLICIES = ("strict", "soft_expected", "hard", "gated_hard")
SCALES = {"square_002um": 3350, "square_008um": 838, "square_016um": 419}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy-mex-root", type=Path, required=True)
    parser.add_argument("--spaceranger-outs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def identity(path: Path) -> dict[str, object]:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256(path)}


def correlation(left: np.ndarray, right: np.ndarray, rank: bool = False) -> float | None:
    if len(left) < 2 or np.std(left) == 0 or np.std(right) == 0:
        return None
    statistic = spearmanr(left, right).statistic if rank else pearsonr(left, right).statistic
    return float(statistic)


def load_spaceranger(
    path: Path, grid_size: int
) -> tuple[list[str], np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(value) for value in matrix["shape"][:])
        if shape[1] != grid_size * grid_size:
            raise ValueError(f"unexpected Space Ranger grid shape in {path}: {shape}")
        feature_ids = [value.decode() for value in matrix["features"]["id"][:]]
        if len(feature_ids) != shape[0] or len(set(feature_ids)) != len(feature_ids):
            raise ValueError(f"invalid Space Ranger feature axis in {path}")

        data = matrix["data"]
        indices = matrix["indices"]
        feature_totals = np.zeros(shape[0], dtype=np.int64)
        for start in range(0, len(data), 2_000_000):
            end = min(len(data), start + 2_000_000)
            feature_totals += np.bincount(
                indices[start:end], weights=data[start:end], minlength=shape[0]
            ).astype(np.int64)

        indptr = matrix["indptr"][:]
        barcode_totals = np.zeros(shape[1], dtype=np.int64)
        for first_column in range(0, shape[1], 100_000):
            last_column = min(shape[1], first_column + 100_000)
            starts = indptr[first_column : last_column + 1]
            values = data[int(starts[0]) : int(starts[-1])]
            cumulative = np.empty(len(values) + 1, dtype=np.int64)
            cumulative[0] = 0
            np.cumsum(values, dtype=np.int64, out=cumulative[1:])
            relative = starts - starts[0]
            barcode_totals[first_column:last_column] = (
                cumulative[relative[1:]] - cumulative[relative[:-1]]
            )

        token = path.parent.name.split("_")[1]
        expected_first = f"s_{token}_00000_00000-1"
        expected_last = f"s_{token}_{grid_size - 1:05d}_{grid_size - 1:05d}-1"
        if matrix["barcodes"][0].decode() != expected_first:
            raise ValueError(f"unexpected first Space Ranger barcode in {path}")
        if matrix["barcodes"][-1].decode() != expected_last:
            raise ValueError(f"unexpected last Space Ranger barcode in {path}")

    if int(feature_totals.sum()) != int(barcode_totals.sum()):
        raise ValueError(f"Space Ranger row and column mass do not reconcile in {path}")
    return feature_ids, feature_totals, barcode_totals


def load_mex(
    directory: Path, scale: str, grid_size: int
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    feature_ids = [
        line.split("\t")[0]
        for line in (directory / "features.tsv").read_text(encoding="utf-8").splitlines()
    ]
    barcodes = (directory / "barcodes.tsv").read_text(encoding="utf-8").splitlines()
    matrix = mmread(directory / "matrix.mtx").tocsr()
    if matrix.shape != (len(feature_ids), len(barcodes)):
        raise ValueError(f"MEX axes do not match matrix dimensions in {directory}")
    feature_totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    barcode_totals = np.asarray(matrix.sum(axis=0)).ravel().astype(float)

    token = scale.split("_")[1]
    pattern = re.compile(rf"^s_{token}_(\d+)_(\d+)-1$")
    coordinate_indices = []
    for barcode in barcodes:
        match = pattern.match(barcode)
        if match is None:
            raise ValueError(f"invalid Visium HD barcode in {directory}: {barcode}")
        row, column = map(int, match.groups())
        if not (0 <= row < grid_size and 0 <= column < grid_size):
            raise ValueError(f"out-of-range Visium HD barcode in {directory}: {barcode}")
        coordinate_indices.append(row * grid_size + column)
    coordinates = np.asarray(coordinate_indices, dtype=np.int64)
    if len(np.unique(coordinates)) != len(coordinates):
        raise ValueError(f"duplicate Visium HD coordinates in {directory}")
    return feature_ids, feature_totals, coordinates, barcode_totals


def main() -> int:
    args = arguments()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite comparison output: {args.output}")

    metrics_path = args.spaceranger_outs / "metrics_summary.csv"
    with metrics_path.open(newline="", encoding="utf-8") as handle:
        metrics_rows = list(csv.DictReader(handle))
    if len(metrics_rows) != 1:
        raise ValueError("Space Ranger metrics summary must contain exactly one sample row")

    space_ranger: dict[str, object] = {
        "outs": str(args.spaceranger_outs.resolve()),
        "metrics_summary": identity(metrics_path),
        "number_of_reads": int(metrics_rows[0]["Number of Reads"]),
    }
    report: dict[str, object] = {
        "schema": "star_suite.visium_hd_spaceranger_sanity.v1",
        "purpose": "descriptive_sanity_only_not_computational_evidence",
        "space_ranger": space_ranger,
        "policy_mex_root": str(args.policy_mex_root.resolve()),
        "comparisons": [],
    }

    comparisons: list[dict[str, object]] = []
    for scale, grid_size in SCALES.items():
        h5_path = args.spaceranger_outs / "binned_outputs" / scale / "raw_feature_bc_matrix.h5"
        sr_ids, sr_gene_totals, sr_bin_totals = load_spaceranger(h5_path, grid_size)
        space_ranger[scale] = identity(h5_path)
        sr_gene_index = {feature: index for index, feature in enumerate(sr_ids)}
        sr_mass = int(sr_gene_totals.sum())
        sr_occupied = int(np.count_nonzero(sr_bin_totals))

        for policy in POLICIES:
            mex_dir = args.policy_mex_root / policy / scale
            our_ids, our_gene_totals, coordinates, our_bin_totals = load_mex(
                mex_dir, scale, grid_size
            )
            aligned = np.zeros(len(sr_ids), dtype=float)
            missing = []
            for feature, value in zip(our_ids, our_gene_totals):
                index = sr_gene_index.get(feature)
                if index is None:
                    missing.append(feature)
                else:
                    aligned[index] += value

            sr_at_our_bins = sr_bin_totals[coordinates].astype(float)
            our_occupied = our_bin_totals > 0
            occupied_overlap = our_occupied & (sr_at_our_bins > 0)
            common_detected = (aligned > 0) & (sr_gene_totals > 0)
            top_ours = set(np.argsort(aligned)[-100:])
            top_sr = set(np.argsort(sr_gene_totals)[-100:])
            our_mass = float(our_gene_totals.sum())
            comparisons.append(
                {
                    "policy": policy,
                    "scale": scale,
                    "our_mass": our_mass,
                    "space_ranger_full_mass": sr_mass,
                    "depth_fraction": our_mass / sr_mass,
                    "our_features": len(our_ids),
                    "missing_feature_ids": len(missing),
                    "common_detected_genes": int(common_detected.sum()),
                    "gene_raw_pearson_all": correlation(aligned, sr_gene_totals),
                    "gene_log1p_pearson_all": correlation(
                        np.log1p(aligned), np.log1p(sr_gene_totals)
                    ),
                    "gene_spearman_common_detected": correlation(
                        aligned[common_detected], sr_gene_totals[common_detected], rank=True
                    ),
                    "top100_gene_overlap": len(top_ours & top_sr),
                    "our_occupied_bins": int(our_occupied.sum()),
                    "space_ranger_full_occupied_bins": sr_occupied,
                    "occupied_bin_overlap": int(occupied_overlap.sum()),
                    "our_occupied_bins_present_in_space_ranger_fraction": float(
                        occupied_overlap.sum() / max(1, our_occupied.sum())
                    ),
                    "bin_raw_pearson_on_our_occupied": correlation(
                        our_bin_totals[our_occupied], sr_at_our_bins[our_occupied]
                    ),
                    "bin_log1p_pearson_on_our_occupied": correlation(
                        np.log1p(our_bin_totals[our_occupied]),
                        np.log1p(sr_at_our_bins[our_occupied]),
                    ),
                }
            )

    report["comparisons"] = comparisons
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(args.output)
    print(json.dumps({"comparisons": len(comparisons), "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
