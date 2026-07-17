#!/usr/bin/env python3
"""Compare native STAR Suite products with the frozen Python reference."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolver", type=Path, required=True)
    parser.add_argument("--processing-repo", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    return parser.parse_args()


def close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1e-10, abs_tol=1e-12)


def load_fixture(path: Path, candidate_read):
    grouped: dict[str, dict] = {}
    counts: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            read = grouped.setdefault(row["read_id"], {
                "feature": row["feature_id"], "umi": row["raw_umi"], "scores": {},
            })
            assert (read["feature"], read["umi"]) == (row["feature_id"], row["raw_umi"])
            read["scores"][row["candidate"]] = float(row["log_sequence_likelihood"])
            count = int(row["exact_read_count"])
            assert row["candidate"] not in counts or counts[row["candidate"]] == count
            counts[row["candidate"]] = count
    reads = [
        candidate_read(read_id, row["feature"], row["umi"], row["scores"])
        for read_id, row in sorted(grouped.items())
    ]
    return reads, counts


def reference_clique_key(clique) -> tuple:
    return (clique.feature_id, clique.umi, clique.read_ids, clique.candidates)


def native_cliques(path: Path):
    by_id: dict[str, list[dict]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_id[row["clique_id"]].append(row)
    output, id_to_key = {}, {}
    for clique_id, rows in by_id.items():
        first = rows[0]
        members = tuple(first["member_read_ids"].split(";"))
        ordered = sorted(rows, key=lambda row: row["candidate"])
        candidates = tuple(row["candidate"] for row in ordered)
        key = (first["feature_id"], first["raw_umi"], members, candidates)
        id_to_key[clique_id] = key
        output[key] = {
            row["candidate"]: tuple(float(row[name]) for name in (
                "log_sequence_likelihood_sum", "log_exact_read_prior", "log_evidence", "posterior",
            ))
            for row in ordered
        }
    return output, id_to_key


def compare_cliques(reference, native) -> None:
    expected = {
        reference_clique_key(clique): {
            candidate: values
            for candidate, values in zip(clique.candidates, zip(
                clique.log_likelihood_sums, clique.log_read_priors,
                clique.log_evidence, clique.posterior, strict=True,
            ), strict=True)
        }
        for clique in reference
    }
    assert expected.keys() == native.keys()
    for key in expected:
        assert expected[key].keys() == native[key].keys()
        for candidate in expected[key]:
            assert all(close(left, right) for left, right in zip(
                expected[key][candidate], native[key][candidate], strict=True,
            )), (key, candidate, expected[key][candidate], native[key][candidate])


def native_soft(path: Path) -> dict[tuple, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            (row["umi_mode"], row["feature_id"], row["corrected_umi"], row["candidate"]):
                float(row["expected_count"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def native_molecules(path: Path) -> set[tuple]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            (row["umi_mode"], row["product"], row["feature_id"], row["corrected_umi"],
             row["candidate"], tuple(row["member_read_ids"].split(";")))
            for row in csv.DictReader(handle, delimiter="\t")
        }


def reference_molecules(cliques, corrections, calls, mode: str, product: str) -> set[tuple]:
    grouped: dict[tuple[str, str, str], list] = defaultdict(list)
    calls_by_id = {call.clique_id: call for call in calls}
    for clique in cliques:
        if product == "strict":
            if len(clique.candidates) != 1:
                continue
            candidate = clique.candidates[0]
        elif product == "hard":
            candidate = min(clique.probabilities(), key=lambda value: (-clique.probabilities()[value], value))
        else:
            candidate = calls_by_id[clique.clique_id].candidate
            if candidate is None:
                continue
        corrected = corrections[(clique.feature_id, candidate, clique.umi)]
        grouped[(clique.feature_id, corrected, candidate)].append(clique)
    output = set()
    for (feature, corrected, candidate), members in grouped.items():
        read_ids = tuple(sorted({read_id for clique in members for read_id in clique.read_ids}))
        output.add((mode, product, feature, corrected, candidate, read_ids))
    return output


def main() -> int:
    args = parse_args()
    sys.path.insert(0, str(args.processing_repo / "src"))
    from star_spatial.hd_probabilistic_umi import (  # pylint: disable=import-error
        CandidateRead, build_read_cliques, corrected_umi_maps, gated_hard_calls,
        weighted_occupancies,
    )

    reads, counts = load_fixture(args.fixture, CandidateRead)
    prior = {candidate: math.log(count + 1.0) for candidate, count in counts.items()}
    cliques = build_read_cliques(reads, log_read_prior=prior, temperature=1.0, prior_beta=1.0)
    calls = gated_hard_calls(cliques, min_posterior=0.95, min_margin=0.90)

    with tempfile.TemporaryDirectory(prefix="star_molecule_first_conformance.") as temp:
        out = Path(temp) / "native"
        subprocess.run([
            str(args.resolver), "--input", str(args.fixture), "--out-dir", str(out),
        ], check=True)
        native, native_ids = native_cliques(out / "read_cliques.tsv")
        compare_cliques(cliques, native)

        expected_soft = {}
        expected_molecules = set()
        for mode in ("1mm_cr", "exact"):
            corrections = corrected_umi_maps(cliques, mode=mode)
            for row in weighted_occupancies(cliques, umi_mode=mode):
                expected_soft[(mode, row.feature_id, row.corrected_umi, row.candidate)] = row.occupancy
            for product in ("strict", "hard", "gated_hard"):
                expected_molecules |= reference_molecules(cliques, corrections, calls, mode, product)

        observed_soft = native_soft(out / "soft_expected_molecules.tsv")
        assert expected_soft.keys() == observed_soft.keys()
        assert all(close(expected_soft[key], observed_soft[key]) for key in expected_soft)
        observed_molecules = set()
        for name in ("strict_molecules.tsv", "hard_molecules.tsv", "gated_hard_molecules.tsv"):
            observed_molecules |= native_molecules(out / name)
        assert expected_molecules == observed_molecules

        reference_by_key = {reference_clique_key(clique): call for clique, call in zip(cliques, calls, strict=True)}
        with (out / "hard_call_audit.tsv").open(newline="", encoding="utf-8") as handle:
            audit_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(audit_rows) == len(reference_by_key)
        for row in audit_rows:
            call = reference_by_key[native_ids[row["clique_id"]]]
            assert row["status"] == call.status
            assert row["candidate"] == (call.candidate or "")
            assert row["reason"] == call.reason
            assert close(float(row["posterior"]), call.posterior)
            assert close(float(row["margin"]), call.margin)

    print("native molecule-first products conform to frozen Python reference")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
