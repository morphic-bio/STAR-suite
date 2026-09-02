#!/usr/bin/env python3
"""Independent FASTQ scanner, bucket partitioner, and exact-dedup counter."""

import argparse
import hashlib
import json
import struct
from collections import defaultdict
from pathlib import Path


def read_table(path, key_column, value_builder):
    result = {}
    with Path(path).open(encoding="ascii") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            try:
                result[fields[key_column]] = value_builder(fields)
            except (IndexError, ValueError) as exc:
                raise SystemExit(f"{path}:{line_number}: malformed row: {exc}") from exc
    return result


def fastq_records(path):
    with Path(path).open(encoding="ascii", newline="") as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not name.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                raise SystemExit(f"{path}: malformed FASTQ record at {name.rstrip()!r}")
            yield name[1:].split()[0], sequence


def pack_umi(sequence):
    value = 0
    for base in sequence:
        value = (value << 2) | "ACGT".index(base)
    return value


def pack_key(cb_idx, umi24, gene_idx, tag_idx):
    if not (0 <= cb_idx < (1 << 20)):
        raise ValueError("CB index does not fit CB20")
    if not (0 <= umi24 < (1 << 24)):
        raise ValueError("UMI does not fit UMI24")
    if not (0 <= gene_idx < (1 << 15)):
        raise ValueError("gene index does not fit gene15")
    if not (0 <= tag_idx < (1 << 5)):
        raise ValueError("tag index does not fit tag5")
    return (cb_idx << 44) | (umi24 << 20) | (gene_idx << 5) | tag_idx


def bucket_for_cb(cb_idx, bucket_count, whitelist_size):
    # Fixed-point high-bit partition over the active whitelist domain. This is
    # equivalent to taking the high log2(P) bits after normalizing CB index to
    # [0, 2^20), and keeps contiguous, balanced CB ranges for any WL size.
    return min(bucket_count - 1, (cb_idx * bucket_count) // whitelist_size)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--tags", required=True)
    parser.add_argument("--probes", required=True)
    parser.add_argument("--bucket-count", required=True, type=int)
    args = parser.parse_args()
    if args.bucket_count <= 0 or args.bucket_count & (args.bucket_count - 1):
        raise SystemExit("bucket count must be a power of two")

    barcode_to_idx = read_table(args.barcodes, 1, lambda row: int(row[0]))
    tag_to_idx = read_table(args.tags, 2, lambda row: int(row[0]))
    probe_by_prefix = read_table(
        args.probes, 2, lambda row: (int(row[0]), int(row[3])))
    whitelist_size = len(barcode_to_idx)

    raw_by_bucket = [[] for _ in range(args.bucket_count)]
    aggregate = {}
    total_records = 0
    r1_iter = fastq_records(args.r1)
    r2_iter = fastq_records(args.r2)
    for left, right in zip(r1_iter, r2_iter):
        r1_name, r1_sequence = left
        r2_name, r2_sequence = right
        if r1_name.replace("/1", "") != r2_name.replace("/2", ""):
            raise SystemExit(f"mate-name mismatch: {r1_name} != {r2_name}")
        cb_idx = barcode_to_idx[r1_sequence[:16]]
        umi24 = pack_umi(r1_sequence[16:28])
        gene_idx, region = probe_by_prefix[r2_sequence[:32]]
        tag_idx = tag_to_idx[r2_sequence[68:76]]
        key = pack_key(cb_idx, umi24, gene_idx, tag_idx)
        encoded = 1 | (region << 30)
        bucket = bucket_for_cb(cb_idx, args.bucket_count, whitelist_size)
        raw_by_bucket[bucket].append(struct.pack("<QI", key, encoded))
        old_count, old_region = aggregate.get(key, (0, 0))
        merged_region = region if old_region == 0 else old_region if region == 0 or region == old_region else 3
        aggregate[key] = (old_count + 1, merged_region)
        total_records += 1

    try:
        next(r1_iter)
        raise SystemExit("R1 has more records than R2")
    except StopIteration:
        pass
    try:
        next(r2_iter)
        raise SystemExit("R2 has more records than R1")
    except StopIteration:
        pass

    # Count distinct UMI keys per final (CB, tag, gene) group.
    umi_sets = defaultdict(set)
    for key in aggregate:
        umi_sets[((key >> 44) & 0xFFFFF, key & 0x1F, (key >> 5) & 0x7FFF)].add(
            (key >> 20) & 0xFFFFFF)

    bucket_rows = []
    for index, chunks in enumerate(raw_by_bucket):
        payload = b"".join(chunks)
        if payload:
            bucket_rows.append({
                "bucket": index,
                "records": len(chunks),
                "sha256": hashlib.sha256(payload).hexdigest(),
            })
    aggregate_payload = b"".join(
        struct.pack("<QI", key, count | (region << 30))
        for key, (count, region) in sorted(aggregate.items())
    )
    matrix_rows = [
        {"cb_idx": cb, "tag_idx": tag, "gene_idx": gene, "umis": len(umis)}
        for (cb, tag, gene), umis in sorted(umi_sets.items())
    ]
    result = {
        "format": "STAR-CB-BUCKET-REFERENCE-v1",
        "bucket_count": args.bucket_count,
        "whitelist_size": whitelist_size,
        "total_records": total_records,
        "unique_records": len(aggregate),
        "aggregate_sha256": hashlib.sha256(aggregate_payload).hexdigest(),
        "buckets": bucket_rows,
        "matrix": matrix_rows,
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
