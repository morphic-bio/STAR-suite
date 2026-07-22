#!/usr/bin/env python3
"""Run native 100K conformance against sealed HD or live frozen Python summaries."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("hd", "scrna"), required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--resolver", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--processing-repo", type=Path, default=Path("/mnt/pikachu/visium-hd-processing"))
    parser.add_argument("--h0-prior", type=Path)
    parser.add_argument("--sealed-summary", type=Path)
    parser.add_argument("--condition", default="challenged")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def opener(path: Path, mode: str):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(mode, encoding="utf-8", newline="")


def load_h0_prior(path: Path) -> dict[tuple[str, int], int]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            (row["barcode_half"], int(row["oligo_index"])): int(row["exact_h0_read_count"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def normalize_hd(source: Path, destination: Path, prior_path: Path) -> int:
    prior = load_h0_prior(prior_path)
    rows = 0
    with opener(source, "rt") as input_handle, destination.open("w", encoding="utf-8", newline="") as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("read_id", "feature_id", "raw_umi", "candidate", "log_sequence_likelihood", "exact_read_count"))
        for row in reader:
            row2, col2 = int(row["row2"]), int(row["col2"])
            synthetic_count = (prior[("BC1", col2)] + 1) * (prior[("BC2", row2)] + 1) - 1
            writer.writerow((row["read_id"], row["feature_id"], row["raw_umi"],
                             f"s_002um_{row2}_{col2}", row["log_sequence_likelihood"], synthetic_count))
            rows += 1
    return rows


def normalize_scrna(source: Path, destination: Path, condition: str) -> tuple[int, list[dict]]:
    rows, records = 0, []
    with opener(source, "rt") as input_handle, destination.open("w", encoding="utf-8", newline="") as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("read_id", "feature_id", "raw_umi", "candidate", "log_sequence_likelihood", "exact_read_count"))
        for row in reader:
            if row["condition"] != condition:
                continue
            candidates = row["candidates"].split(",")
            likelihoods = [float(value) for value in row["log_likelihoods"].split(",")]
            if not candidates or len(candidates) != len(likelihoods):
                raise ValueError(f"invalid candidate row: {row['read_id']}")
            records.append({
                "read_id": row["read_id"], "feature": row["feature_signature"],
                "umi": row["raw_umi"], "scores": dict(zip(candidates, likelihoods, strict=True)),
            })
            for candidate, likelihood in zip(candidates, likelihoods, strict=True):
                writer.writerow((row["read_id"], row["feature_signature"], row["raw_umi"],
                                 candidate, repr(likelihood), 0))
                rows += 1
    return rows, records


def parse_summary(path: Path) -> dict[str, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        output = {}
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                output[row["metric"]] = float(row["value"])
            except ValueError:
                continue
        return output


def molecule_count(cliques, corrections, calls, product: str) -> int:
    call_by_id = {call.clique_id: call for call in calls}
    groups = set()
    for clique in cliques:
        if product == "strict":
            if len(clique.candidates) != 1:
                continue
            candidate = clique.candidates[0]
        elif product == "hard":
            candidate = min(clique.probabilities(), key=lambda value: (-clique.probabilities()[value], value))
        else:
            candidate = call_by_id[clique.clique_id].candidate
            if candidate is None:
                continue
        groups.add((clique.feature_id, corrections[(clique.feature_id, candidate, clique.umi)], candidate))
    return len(groups)


def expected_scrna(records: list[dict], processing_repo: Path) -> dict[str, float]:
    sys.path.insert(0, str(processing_repo / "src"))
    from star_spatial.hd_probabilistic_umi import (  # pylint: disable=import-error
        CandidateRead, build_read_cliques, corrected_umi_maps, gated_hard_calls,
        weighted_occupancies,
    )
    reads = [CandidateRead(row["read_id"], row["feature"], row["umi"], row["scores"]) for row in records]
    cliques = build_read_cliques(reads)
    calls = gated_hard_calls(cliques, min_posterior=0.95, min_margin=0.90)
    output: dict[str, float] = {
        "candidate_reads": len(reads), "read_cliques": len(cliques),
        "gated_assigned": sum(call.candidate is not None for call in calls),
        "gated_deferred": sum(call.candidate is None for call in calls),
    }
    for mode in ("1mm_cr", "exact"):
        corrections = corrected_umi_maps(cliques, mode=mode)
        occupancies = weighted_occupancies(cliques, umi_mode=mode)
        output[f"soft_rows_{mode}"] = len(occupancies)
        output[f"soft_mass_{mode}"] = math.fsum(row.occupancy for row in occupancies)
        for product in ("strict", "hard", "gated_hard"):
            output[f"molecules_{mode}_{product}"] = molecule_count(cliques, corrections, calls, product)
    return output


def expected_hd(path: Path) -> dict[str, float]:
    data = json.loads(path.read_text(encoding="utf-8"))["counts"]
    return {key: float(value) for key, value in data.items()}


def compare(native: dict[str, float], expected: dict[str, float]) -> dict[str, float]:
    mapping = {
        "candidate_reads": "candidate_reads",
        "read_cliques": "read_cliques",
        "gated_assigned": "gated_assigned_cliques",
        "gated_deferred": "gated_deferred_cliques",
    }
    for mode in ("1mm_cr", "exact"):
        mapping[f"soft_rows_{mode}"] = f"{mode}.soft_rows_count"
        mapping[f"soft_mass_{mode}"] = f"{mode}.soft_occupancy_mass"
        mapping[f"molecules_{mode}_strict"] = f"{mode}.strict_count"
        mapping[f"molecules_{mode}_hard"] = f"{mode}.hard_count"
        mapping[f"molecules_{mode}_postcollapse_hard"] = f"{mode}.hard_count"
        mapping[f"molecules_{mode}_gated_hard"] = f"{mode}.gated_hard_count"
    residuals = {}
    for expected_key, native_key in mapping.items():
        if expected_key not in expected:
            continue
        if native_key not in native:
            raise AssertionError(f"native summary missing {native_key}")
        residual = native[native_key] - expected[expected_key]
        tolerance = 1e-8 * max(1.0, abs(expected[expected_key])) if "mass" in expected_key else 1e-9
        if abs(residual) > tolerance:
            raise AssertionError((expected_key, native_key, expected[expected_key], native[native_key], residual))
        residuals[expected_key] = residual
    return residuals


def main() -> int:
    args = parse_args()
    if args.kind == "hd" and (args.h0_prior is None or args.sealed_summary is None):
        raise SystemExit("--h0-prior and --sealed-summary are required for HD")
    args.out_dir.mkdir(parents=True, exist_ok=False)
    normalized = args.out_dir / "candidate_reads.normalized.tsv"
    records = []
    if args.kind == "hd":
        candidate_rows = normalize_hd(args.input, normalized, args.h0_prior)
        expected = expected_hd(args.sealed_summary)
    else:
        candidate_rows, records = normalize_scrna(args.input, normalized, args.condition)
        expected = expected_scrna(records, args.processing_repo)
    native_dir = args.out_dir / "native"
    subprocess.run([str(args.resolver), "--input", str(normalized), "--out-dir", str(native_dir)], check=True)
    native = parse_summary(native_dir / "summary.tsv")
    residuals = compare(native, expected)
    result = {
        "schema": "star_suite.molecule_first.100k_conformance.v1",
        "name": args.name,
        "kind": args.kind,
        "status": "pass",
        "input": {"path": str(args.input.resolve()), "sha256": sha256(args.input)},
        "normalized": {"path": str(normalized.resolve()), "sha256": sha256(normalized), "candidate_rows": candidate_rows},
        "h0_prior": None if args.h0_prior is None else {"path": str(args.h0_prior.resolve()), "sha256": sha256(args.h0_prior)},
        "sealed_summary": None if args.sealed_summary is None else {"path": str(args.sealed_summary.resolve()), "sha256": sha256(args.sealed_summary)},
        "native_summary": native,
        "expected_summary": expected,
        "residuals": residuals,
    }
    (args.out_dir / "conformance.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"name": args.name, "status": "pass", "residuals": residuals}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
