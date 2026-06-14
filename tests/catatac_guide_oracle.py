#!/usr/bin/env python3
"""Oracle and native split-read parity checks for CAT-ATAC guide capture."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def read_fastq_records(path: Path, limit: int | None = None):
    n = 0
    with gzip.open(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            fh.readline()
            qual = fh.readline().strip()
            yield seq, qual
            n += 1
            if limit is not None and n >= limit:
                break


def load_whitelist(path: Path) -> set[str]:
    out: set[str] = set()
    for line in path.read_text().splitlines():
        token = line.split()[0].upper()
        if token:
            out.add(token)
    return out


def load_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            mapping[parts[0].upper()] = parts[1].upper()
    return mapping


def barcode_preflight(r2: Path, whitelist: set[str], limit: int = 200000) -> dict:
    checks = {(8, 16): 0, (9, 16): 0, (8, 24): 0}
    n = 0
    for seq, _ in read_fastq_records(r2, limit=limit):
        n += 1
        for (off, length), hits in checks.items():
            sub = seq[off : off + length]
            if len(sub) == length and revcomp(sub).upper() in whitelist:
                checks[(off, length)] = hits + 1
    return {
        "n": n,
        "checks": {f"offset={k[0]}_length={k[1]}": v / n if n else 0.0 for k, v in checks.items()},
    }


def capture_preflight(r1: Path, limit: int = 200000) -> dict:
    cs1 = "CAAGTTGATAACGGACTAGCC"
    cs2 = "CAAGTTGTAAACGGACTAGCC"
    n = cs1_hits = cs2_hits = either_hits = 0
    for seq, _ in read_fastq_records(r1, limit=limit):
        n += 1
        a = cs1 in seq
        b = cs2 in seq
        cs1_hits += int(a)
        cs2_hits += int(b)
        either_hits += int(a or b)
    return {
        "n": n,
        "cs1_rate": cs1_hits / n if n else 0.0,
        "cs2_rate": cs2_hits / n if n else 0.0,
        "either_rate": either_hits / n if n else 0.0,
    }


def materialize_oracle_fastqs(fixture_dir: Path, out_dir: Path, limit: int | None = None) -> dict:
    """Materialize CB+UMI oracle FASTQs.

    ``limit`` caps emitted captured reads (matching native ``max_reads``), not raw triplets.
    """
    r1 = fixture_dir / "guide_R1.fastq.gz"
    r2 = fixture_dir / "guide_R2.fastq.gz"
    r3 = fixture_dir / "guide_R3.fastq.gz"
    cs1 = "CAAGTTGATAACGGACTAGCC"
    cs2 = "CAAGTTGTAAACGGACTAGCC"
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics = {
        "total_reads": 0,
        "capture_either_hits": 0,
        "barcode_synth_ok": 0,
    }
    with gzip.open(out_dir / "oracle_R1_.fastq.gz", "wt") as bc_out, gzip.open(
        out_dir / "oracle_R2_.fastq.gz", "wt"
    ) as feat_out:
        for (r1_seq, r1_qual), (r2_seq, r2_qual), (r3_seq, r3_qual) in zip(
            read_fastq_records(r1), read_fastq_records(r2), read_fastq_records(r3)
        ):
            metrics["total_reads"] += 1
            if cs1 not in r1_seq and cs2 not in r1_seq:
                continue
            metrics["capture_either_hits"] += 1
            cb = revcomp(r2_seq[8:24])
            umi = r1_seq[0:12]
            if len(cb) != 16 or len(umi) != 12:
                continue
            bc_qual = r2_qual[8:24][::-1]
            umi_qual = r1_qual[0:12]
            synth = cb + umi
            synth_qual = bc_qual + umi_qual
            metrics["barcode_synth_ok"] += 1
            idx = metrics["barcode_synth_ok"]
            bc_out.write(f"@oracle_{idx}\n{synth}\n+\n{synth_qual}\n")
            feat_out.write(f"@oracle_{idx}\n{r3_seq}\n+\n{r3_qual}\n")
            if limit is not None and metrics["barcode_synth_ok"] >= limit:
                break
    return metrics


def run_assign(assign_bin: Path, whitelist: Path, feature_ref: Path, fastq_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(assign_bin),
        "-w",
        str(whitelist),
        "-f",
        str(feature_ref),
        "-d",
        str(out_dir),
        str(fastq_dir),
        "--maxHammingDistance",
        "1",
        "--limit_search",
        "-1",
        "--min_counts",
        "0",
        "--skip_heatmaps",
        "--consumer_threads_per_set",
        "1",
        "-S",
        "1",
        "--barcode_fastq_pattern",
        "_R1_",
        "--forward_fastq_pattern",
        "_R2_",
    ]
    subprocess.run(cmd, check=True)


def run_split_native(split_bin: Path, fixture_dir: Path, out_dir: Path, whitelist: Path, feature_ref: Path, limit: int) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["CATATAC_GUIDE_MAX_READS"] = str(limit)
    subprocess.run(
        [str(split_bin), str(fixture_dir), str(out_dir), str(whitelist), str(feature_ref)],
        check=True,
        env=env,
    )


def _read_name_list(mex_dir: Path, stem: str) -> list[str]:
    for suffix in (".tsv", ".txt"):
        path = mex_dir / f"{stem}{suffix}"
        if path.is_file():
            return [line.split()[0] for line in path.read_text().splitlines() if line.strip()]
    raise FileNotFoundError(f"missing {stem}.tsv or {stem}.txt under {mex_dir}")


def load_mex_triplet(mex_dir: Path):
    features = _read_name_list(mex_dir, "features")
    barcodes = _read_name_list(mex_dir, "barcodes")
    matrix_path = mex_dir / "matrix.mtx"
    return features, barcodes, matrix_path.read_bytes()


def apply_barcode_output_map(barcodes_txt: Path, output_map: dict[str, str]) -> list[str]:
    """Mirror PfMultiMexStub::copyBarcodesTsv col1->col2 remap."""
    mapped: list[str] = []
    for line in barcodes_txt.read_text().splitlines():
        token = line.split()[0].upper()
        if not token:
            continue
        mapped.append(output_map.get(token, token))
    return mapped


def verify_gex_barcodes_tsv(
    assign_dir: Path,
    output_map: dict[str, str],
    gex_whitelist: set[str],
) -> dict:
    """Assert stub-mapped barcodes.tsv entries land in GEX whitelist space."""
    barcodes_txt = assign_dir / "barcodes.txt"
    if not barcodes_txt.is_file():
        raise FileNotFoundError(f"missing {barcodes_txt}")
    atac_barcodes = [line.split()[0].upper() for line in barcodes_txt.read_text().splitlines() if line.strip()]
    gex_barcodes = apply_barcode_output_map(barcodes_txt, output_map)
    gex_hits = sum(1 for bc in gex_barcodes if bc in gex_whitelist)
    atac_in_gex = sum(1 for bc in atac_barcodes if bc in gex_whitelist)
    mapped_diff = sum(1 for a, g in zip(atac_barcodes, gex_barcodes) if a != g)
    barcodes_tsv = assign_dir / "barcodes.tsv"
    barcodes_tsv.write_text("\n".join(gex_barcodes) + "\n")
    return {
        "barcode_count": len(gex_barcodes),
        "gex_whitelist_hits": gex_hits,
        "atac_direct_gex_hits": atac_in_gex,
        "mapped_diff_from_atac": mapped_diff,
        "barcodes_tsv": str(barcodes_tsv),
    }


def load_mex_counts(mex_dir: Path) -> dict[tuple[str, str], int]:
    features, barcodes, _ = load_mex_triplet(mex_dir)
    counts: dict[tuple[str, str], int] = {}
    with (mex_dir / "matrix.mtx").open() as fh:
        saw_dims = False
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            if not saw_dims:
                saw_dims = True
                continue
            feat_idx, bc_idx, value = int(parts[0]), int(parts[1]), int(float(parts[2]))
            key = (features[feat_idx - 1], barcodes[bc_idx - 1])
            counts[key] = counts.get(key, 0) + value
    return counts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture-dir", required=True)
    parser.add_argument("--feature-ref", required=True)
    parser.add_argument("--whitelist", required=True)
    parser.add_argument("--output-map", required=True)
    parser.add_argument("--oracle-out", required=True)
    parser.add_argument("--native-out", required=True)
    parser.add_argument("--assign-bin", default="core/features/process_features/assignBarcodes")
    parser.add_argument("--split-bin", default="core/features/process_features/tests/test_catatac_split_read")
    parser.add_argument("--limit", type=int, default=200000)
    parser.add_argument(
        "--gex-whitelist",
        default="/mnt/pikachu/GEX_whitelist/737K-arc-v1.txt",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    assign_bin = (repo_root / args.assign_bin).resolve()
    split_bin = (repo_root / args.split_bin).resolve()
    fixture_dir = Path(args.fixture_dir)
    feature_ref = Path(args.feature_ref)
    whitelist = Path(args.whitelist)
    output_map = Path(args.output_map)
    oracle_out = Path(args.oracle_out)
    native_out = Path(args.native_out)

    wl = load_whitelist(whitelist)
    preflight = {
        "barcode_window": barcode_preflight(fixture_dir / "guide_R2.fastq.gz", wl, limit=args.limit),
        "capture": capture_preflight(fixture_dir / "guide_R1.fastq.gz", limit=args.limit),
    }
    print(json.dumps(preflight, indent=2))

    off816 = preflight["barcode_window"]["checks"]["offset=8_length=16"]
    off916 = preflight["barcode_window"]["checks"]["offset=9_length=16"]
    off824 = preflight["barcode_window"]["checks"]["offset=8_length=24"]
    if off816 < 0.5:
        raise SystemExit(f"barcode preflight failed: offset=8 length=16 rate={off816}")
    if off916 > 0.01:
        raise SystemExit(f"barcode preflight failed: offset=9 length=16 rate={off916}")
    if off824 > 0.01:
        raise SystemExit(f"barcode preflight failed: offset=8 length=24 rate={off824}")
    if preflight["capture"]["cs2_rate"] <= 0:
        raise SystemExit("capture preflight failed: CS2 not accepted")

    with tempfile.TemporaryDirectory(prefix="catatac_oracle_") as tmp:
        oracle_fastqs = Path(tmp) / "oracle_fastqs"
        oracle_metrics = materialize_oracle_fastqs(fixture_dir, oracle_fastqs, limit=args.limit)
        run_assign(assign_bin, whitelist, feature_ref, oracle_fastqs, oracle_out / "assign")

    run_split_native(split_bin, fixture_dir, native_out / "split_assign", whitelist, feature_ref, args.limit)
    print("native split-read assignment complete", flush=True)

    oracle_mex = oracle_out / "assign" / "oracle_fastqs"
    native_mex = native_out / "split_assign" / "sample"
    oracle_features, oracle_barcodes, _ = load_mex_triplet(oracle_mex)
    native_features, native_barcodes, _ = load_mex_triplet(native_mex)
    oracle_counts = load_mex_counts(oracle_mex)
    native_counts = load_mex_counts(native_mex)

    gex_map = load_map(output_map)
    gex_whitelist = load_whitelist(Path(args.gex_whitelist))
    gex_stub = verify_gex_barcodes_tsv(native_mex, gex_map, gex_whitelist)
    gex_tsv_barcodes = [
        line.split()[0].upper()
        for line in (native_mex / "barcodes.tsv").read_text().splitlines()
        if line.strip()
    ]
    gex_tsv_hits = sum(1 for bc in gex_tsv_barcodes if bc in gex_whitelist)

    only_oracle = len(set(oracle_counts) - set(native_counts))
    only_native = len(set(native_counts) - set(oracle_counts))
    value_mismatches = sum(
        1
        for key in set(oracle_counts) & set(native_counts)
        if oracle_counts[key] != native_counts[key]
    )

    parity = {
        "oracle_metrics": oracle_metrics,
        "features_equal": oracle_features == native_features,
        "barcodes_equal": sorted(oracle_barcodes) == sorted(native_barcodes),
        "matrix_equal": oracle_counts == native_counts,
        "matrix_only_oracle": only_oracle,
        "matrix_only_native": only_native,
        "matrix_value_mismatches": value_mismatches,
        "native_barcode_count": len(native_barcodes),
        "native_feature_count": len(native_features),
        "oracle_barcode_count": len(oracle_barcodes),
        "oracle_feature_count": len(oracle_features),
        "gex_stub": gex_stub,
        "gex_tsv_whitelist_hits": gex_tsv_hits,
    }
    (native_out / "parity.json").write_text(json.dumps(parity, indent=2) + "\n")
    print(json.dumps(parity, indent=2))

    if not native_features or not native_barcodes:
        raise SystemExit("native split-read output is empty")
    if gex_tsv_hits != len(gex_tsv_barcodes):
        raise SystemExit("barcodes.tsv stub remap did not land entirely in GEX whitelist")
    if gex_stub["mapped_diff_from_atac"] == 0:
        raise SystemExit("barcodes.tsv stub remap did not translate ATAC barcodes to GEX")
    if not parity["features_equal"] or not parity["matrix_equal"]:
        raise SystemExit("oracle/native MEX parity failed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
