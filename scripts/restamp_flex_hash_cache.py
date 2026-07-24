#!/usr/bin/env python3
"""Convert a legacy FH01SEQ1 v1 cache to sample-aware v2.

The v1 cache stores zero in the final uint16 field of every record.  A cache
calibrated with one sample tag must stamp its exact-match H0 rows with that
sample's sequential index when it is used with a larger whitelist.  H1/H2 and
negative rows remain global (sample index zero), matching the production cache
generator contract.  Rows for other samples then fall back to normal alignment
instead of consuming H0 evidence calibrated on the source sample.
"""

import argparse
import hashlib
import os
import struct
import tempfile


MAGIC = b"FH01SEQ1"
HEADER = struct.Struct("<8sHHIQ")
RECORD = struct.Struct("<QQIBBH")
KMER_LENGTH = 50
VERSION_LEGACY = 1
VERSION_SAMPLE_AWARE = 2


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def restamp(
    source,
    output,
    sample_index,
    expected_source_sha256=None,
    include_classes=None,
    stamp_classes=(0,),
):
    sample_index = int(sample_index)
    if not 1 <= sample_index <= 0xFFFF:
        raise ValueError("sample index must be in [1, 65535]")
    included = None if include_classes is None else {int(value) for value in include_classes}
    if included is not None and not included:
        raise ValueError("include_classes must not be empty")
    stamped = {int(value) for value in stamp_classes}
    if not stamped:
        raise ValueError("stamp_classes must not be empty")
    if 0 not in stamped:
        raise ValueError("H0 cache class 0 must be sample-stamped")
    if os.path.abspath(source) == os.path.abspath(output):
        raise ValueError("source and output paths must differ")
    if os.path.exists(output):
        raise FileExistsError(f"refusing to overwrite existing output: {output}")

    source_sha256 = sha256_file(source)
    if expected_source_sha256 and source_sha256 != expected_source_sha256.lower():
        raise ValueError(
            f"source SHA-256 mismatch: expected {expected_source_sha256.lower()}, got {source_sha256}"
        )

    output_dir = os.path.dirname(os.path.abspath(output))
    os.makedirs(output_dir, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".restamp_flex_hash_cache.", dir=output_dir)
    os.close(fd)

    class_counts = {}
    stamped_class_counts = {}
    try:
        with open(source, "rb") as source_handle, open(temporary, "wb") as output_handle:
            header = source_handle.read(HEADER.size)
            if len(header) != HEADER.size:
                raise ValueError("truncated FH01SEQ1 header")
            magic, version, kmer_length, record_size, record_count = HEADER.unpack(header)
            if magic != MAGIC:
                raise ValueError("invalid FH01SEQ1 magic")
            if version != VERSION_LEGACY:
                raise ValueError(f"expected FH01SEQ1 v1 source, found v{version}")
            if kmer_length != KMER_LENGTH or record_size != RECORD.size:
                raise ValueError(
                    f"unsupported FH01SEQ1 layout: kmer={kmer_length}, record_size={record_size}"
                )

            expected_size = HEADER.size + record_count * RECORD.size
            actual_size = os.fstat(source_handle.fileno()).st_size
            if actual_size != expected_size:
                raise ValueError(
                    f"FH01SEQ1 size mismatch: header expects {expected_size} bytes, file has {actual_size}"
                )

            output_handle.write(b"\0" * HEADER.size)
            for record_number in range(record_count):
                raw = source_handle.read(RECORD.size)
                if len(raw) != RECORD.size:
                    raise ValueError(f"truncated FH01SEQ1 record {record_number}")
                seq_lo, seq_hi, gene_idx15, cache_class, negative_code, reserved = RECORD.unpack(raw)
                if reserved != 0:
                    raise ValueError(
                        f"v1 reserved field is nonzero at record {record_number}: {reserved}"
                    )
                if included is not None and cache_class not in included:
                    continue
                class_counts[cache_class] = class_counts.get(cache_class, 0) + 1
                output_sample_index = sample_index if cache_class in stamped else 0
                if output_sample_index:
                    stamped_class_counts[cache_class] = stamped_class_counts.get(cache_class, 0) + 1
                output_handle.write(
                    RECORD.pack(
                        seq_lo,
                        seq_hi,
                        gene_idx15,
                        cache_class,
                        negative_code,
                        output_sample_index,
                    )
                )
        if class_counts.get(0, 0) == 0:
            raise ValueError("source cache contains no H0 rows to restamp")
        output_record_count = sum(class_counts.values())
        with open(temporary, "r+b") as output_handle:
            output_handle.write(
                HEADER.pack(
                    MAGIC,
                    VERSION_SAMPLE_AWARE,
                    KMER_LENGTH,
                    RECORD.size,
                    output_record_count,
                )
            )
        os.replace(temporary, output)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise

    return {
        "source_sha256": source_sha256,
        "output_sha256": sha256_file(output),
        "record_count": sum(class_counts.values()),
        "class_counts": ",".join(
            f"{cache_class}:{count}" for cache_class, count in sorted(class_counts.items())
        ),
        "stamped_class_counts": ",".join(
            f"{cache_class}:{count}"
            for cache_class, count in sorted(stamped_class_counts.items())
        ),
        "sample_index": sample_index,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Restamp all FH01SEQ1 v1 rows with one sample index and write v2."
    )
    parser.add_argument("--input", required=True, help="legacy FH01SEQ1 v1 cache")
    parser.add_argument("--output", required=True, help="new FH01SEQ1 v2 cache")
    parser.add_argument("--sample-index", required=True, type=int)
    parser.add_argument(
        "--include-classes",
        help="optional comma-separated cache classes to retain (for example, 0 for H0-only)",
    )
    parser.add_argument(
        "--stamp-classes",
        default="0",
        help="comma-separated cache classes to sample-stamp (default: 0/H0 only)",
    )
    parser.add_argument("--expected-input-sha256")
    args = parser.parse_args()

    result = restamp(
        args.input,
        args.output,
        args.sample_index,
        args.expected_input_sha256,
        None if args.include_classes is None else [
            int(value) for value in args.include_classes.split(",") if value
        ],
        [int(value) for value in args.stamp_classes.split(",") if value],
    )
    for key in (
        "source_sha256",
        "output_sha256",
        "record_count",
        "class_counts",
        "stamped_class_counts",
        "sample_index",
    ):
        print(f"{key}={result[key]}")


if __name__ == "__main__":
    main()
