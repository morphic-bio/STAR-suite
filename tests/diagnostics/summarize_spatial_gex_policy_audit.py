#!/usr/bin/env python3
"""Summarize retained native stage accounting and bounded resolver diagnostics.

This reads existing run metadata, not large matrices, and executes no STAR job.
It is diagnostic evidence, not a generator of replacement publication timings.
"""
import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_primary(root):
    native = root / "starSpatialGex.out/run_summary.tsv"
    summary = root / "starSpatialGex.out/summary.tsv"
    timing = root / "logs/star.time.txt"
    manifest = root / "mex.relative.sha256"
    values = dict(line.split("\t", 1) for line in native.read_text().splitlines())
    lines = summary.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("product\tscale\t"))
    rows = [dict(zip(lines[start].split("\t"), line.split("\t")))
            for line in lines[start + 1:]]
    duration = re.search(r"Elapsed .*?:\s*([\d:.]+)\s*$", timing.read_text(), re.M).group(1)
    elapsed = 0.0
    for part in duration.split(":"):
        elapsed = elapsed * 60 + float(part)
    hashes = {}
    for line in manifest.read_text().splitlines():
        digest, name = line.split(None, 1)
        hashes[name.removeprefix("./")] = digest
    expected_matrix_bytes = sum(int(row["nnz"]) for row in rows) * 16 + int(values["downstream_matrix_runs"]) * 104
    expected_contribution_bytes = int(values["downstream_contribution_records"]) * 32 + int(values["downstream_contribution_runs"]) * 104
    assert expected_matrix_bytes == int(values["downstream_matrix_bytes"])
    assert expected_contribution_bytes == int(values["downstream_contribution_bytes"])
    return {
        "root": str(root), "source_sha256": {str(p): sha(p) for p in [native, summary, timing, manifest]},
        "native": values, "policy_rows": rows, "elapsed_seconds": elapsed,
        "mex_hashes": hashes,
        "matrix_spool_bytes_exactly_accounted_for": expected_matrix_bytes,
        "contribution_spool_bytes_exactly_accounted_for": expected_contribution_bytes,
        "scope": "Sealed historical MEX manifests compared; no fresh large-payload hashing or execution.",
    }


def read_diagnostic(path):
    probes = []
    isolation_pass = False
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if parts[:2] == ["policy_isolation", "PASS"]:
            isolation_pass = True
        if parts[0] == "soft_order_probe":
            probes.append(dict(zip(parts[1::2], parts[2::2])))
    assert isolation_pass and len(probes) == 65
    return {"path": str(path), "sha256": sha(path), "policy_isolation_passed": True,
            "probes": probes, "distinct_target_masses": sorted({float(p["target_mass"]) for p in probes})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v161", type=Path, required=True)
    parser.add_argument("--v171", type=Path, required=True)
    parser.add_argument("--v195", type=Path, required=True)
    parser.add_argument("--diagnostic-logs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    primary = {version: load_primary(getattr(args, version)) for version in ["v161", "v171", "v195"]}
    current = primary["v195"]
    keys = ["reads_decoded", "reads_with_candidates", "unique_gene_reads", "joined_reads",
            "candidate_rows", "exact_h0_reads", "barcode_reads_with_n", "barcode_n_bases",
            "barcode_dp_recovered_reads", "barcode_dp_ambiguous_reads", "barcode_dp_unassigned_reads",
            "barcode_unsupported_reads", "umi_reads_with_n", "umi_reads_with_invalid_base",
            "read_cliques", "strict_molecules", "hard_molecules", "downstream_contribution_records",
            "downstream_contribution_runs", "downstream_contribution_bytes", "downstream_largest_shard_records"]
    comparisons = {}
    for version in ["v161", "v171"]:
        old = primary[version]
        assert all(old["native"][k] == current["native"][k] for k in keys)
        assert all(current["mex_hashes"].get(k) == v for k, v in old["mex_hashes"].items())
        comparisons[version] = {"identical_native_fields": keys,
                                "matching_sealed_mex_components": len(old["mex_hashes"])}
    stages = ["spill_merge_seconds", "downstream_resolve_seconds", "downstream_materialize_seconds"]
    deltas = {k: float(current["native"][k]) - float(primary["v171"]["native"][k]) for k in stages}
    deltas["total_elapsed_seconds"] = current["elapsed_seconds"] - primary["v171"]["elapsed_seconds"]
    deltas["all_other_work_seconds"] = deltas["total_elapsed_seconds"] - sum(deltas[k] for k in stages)
    diagnostics = {v: read_diagnostic(args.diagnostic_logs / ("policy_audit_" + v + ".tsv"))
                   for v in ["v171", "v195_valid", "v195_strict_sort_control"]}
    assert diagnostics["v171"]["probes"] == diagnostics["v195_valid"]["probes"]
    assert len(diagnostics["v195_valid"]["distinct_target_masses"]) == 2
    assert len(diagnostics["v195_strict_sort_control"]["distinct_target_masses"]) == 1
    result = {"schema": "star_suite.spatial_gex_policy_audit.v1",
              "created_at": datetime.now(timezone.utc).isoformat(),
              "primaries": primary, "comparisons_to_v195": comparisons,
              "v195_minus_v171_seconds": deltas, "diagnostics": diagnostics,
              "conclusion": "No changed upstream counts or duplicated spool records found. Extra policies explain additional work; an old soft-policy ordering defect is independently reproduced.",
              "limits": ["No matched full-slide repeat performed.",
                         "Constructed near-tie defect is not shown to occur in ovarian data.",
                         "Strict-sort control changes only an external diagnostic source copy; no production fix applied."]}
    with args.out.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"output": str(args.out), "comparisons": comparisons,
                      "time_deltas": deltas}, indent=2))


if __name__ == "__main__":
    main()
