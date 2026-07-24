#!/usr/bin/env python3

import hashlib
import importlib.util
import os
import struct
import tempfile


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "scripts", "restamp_flex_hash_cache.py")
SPEC = importlib.util.spec_from_file_location("restamp_flex_hash_cache", TOOL)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def write_cache(path, version, records):
    with open(path, "wb") as handle:
        handle.write(
            MODULE.HEADER.pack(
                MODULE.MAGIC,
                version,
                MODULE.KMER_LENGTH,
                MODULE.RECORD.size,
                len(records),
            )
        )
        for record in records:
            handle.write(MODULE.RECORD.pack(*record))


def read_cache(path):
    with open(path, "rb") as handle:
        header = MODULE.HEADER.unpack(handle.read(MODULE.HEADER.size))
        records = []
        while True:
            raw = handle.read(MODULE.RECORD.size)
            if not raw:
                break
            records.append(MODULE.RECORD.unpack(raw))
    return header, records


with tempfile.TemporaryDirectory(prefix="test_restamp_flex_hash_cache.") as tmp:
    source = os.path.join(tmp, "source.bin")
    output = os.path.join(tmp, "output.bin")
    records = [
        (1, 10, 101, 0, 0, 0),
        (2, 20, 102, 1, 0, 0),
        (3, 30, 0, 2, 1, 0),
        (4, 40, 104, 0, 0, 0),
    ]
    write_cache(source, MODULE.VERSION_LEGACY, records)
    source_sha = hashlib.sha256(open(source, "rb").read()).hexdigest()

    result = MODULE.restamp(source, output, 4, source_sha)
    header, observed = read_cache(output)
    assert header[1] == MODULE.VERSION_SAMPLE_AWARE
    assert header[4] == 4
    assert [record[-1] for record in observed] == [4, 0, 0, 4]
    assert result["class_counts"] == "0:2,1:1,2:1"
    assert result["stamped_class_counts"] == "0:2"
    assert result["record_count"] == 4

    h0_output = os.path.join(tmp, "h0-only.bin")
    h0_result = MODULE.restamp(source, h0_output, 4, source_sha, [0])
    h0_header, h0_observed = read_cache(h0_output)
    assert h0_header[4] == 2
    assert [record[3] for record in h0_observed] == [0, 0]
    assert [record[-1] for record in h0_observed] == [4, 4]
    assert h0_result["class_counts"] == "0:2"

    all_stamped_output = os.path.join(tmp, "all-stamped.bin")
    all_stamped_result = MODULE.restamp(
        source, all_stamped_output, 4, source_sha, stamp_classes=[0, 1, 2]
    )
    _, all_stamped_observed = read_cache(all_stamped_output)
    assert [record[-1] for record in all_stamped_observed] == [4, 4, 4, 4]
    assert all_stamped_result["stamped_class_counts"] == "0:2,1:1,2:1"

    try:
        MODULE.restamp(source, os.path.join(tmp, "bad-sha.bin"), 4, "0" * 64)
        raise AssertionError("expected checksum mismatch")
    except ValueError as error:
        assert "SHA-256 mismatch" in str(error)

    version2 = os.path.join(tmp, "version2.bin")
    write_cache(version2, MODULE.VERSION_SAMPLE_AWARE, records)
    try:
        MODULE.restamp(version2, os.path.join(tmp, "double.bin"), 4)
        raise AssertionError("expected v2 source rejection")
    except ValueError as error:
        assert "expected FH01SEQ1 v1" in str(error)

print("Flex hash-cache restamp tests passed")
