#!/usr/bin/env python3
"""Independent BGZF/FASTQ scanner used as the gold oracle for ingest tests."""

import argparse
import bisect
import json
import struct
import zlib
from pathlib import Path

BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)
FNV_OFFSET = 14695981039346656037
FNV_PRIME = 1099511628211
MASK64 = (1 << 64) - 1


def fnv1a(payload):
    value = FNV_OFFSET
    for byte in payload:
        value ^= byte
        value = (value * FNV_PRIME) & MASK64
    return value


def find_bsize(header):
    if len(header) < 12 or header[:3] != b"\x1f\x8b\x08" or not header[3] & 4:
        return None, None
    xlen = struct.unpack_from("<H", header, 10)[0]
    if len(header) < 12 + xlen:
        return None, None
    pos = 12
    end = 12 + xlen
    while pos + 4 <= end:
        si1, si2, slen = struct.unpack_from("<BBH", header, pos)
        pos += 4
        if pos + slen > end:
            return None, None
        if si1 == ord("B") and si2 == ord("C") and slen == 2:
            return struct.unpack_from("<H", header, pos)[0] + 1, end
        pos += slen
    return None, None


def parse_fastq(payload):
    lines = payload.splitlines()
    if not lines:
        return []
    if len(lines) % 4:
        raise ValueError(f"FASTQ line count {len(lines)} is not divisible by four")
    records = []
    global_offset = 0
    for idx in range(0, len(lines), 4):
        name, seq, plus, qual = lines[idx:idx + 4]
        if not name.startswith(b"@"):
            raise ValueError(f"FASTQ record {idx // 4} does not start with @")
        if not plus.startswith(b"+"):
            raise ValueError(f"FASTQ record {idx // 4} has an invalid plus line")
        if len(seq) != len(qual):
            raise ValueError(f"FASTQ record {idx // 4} has unequal sequence/quality lengths")
        records.append((global_offset, name[1:], seq, qual))
        global_offset += sum(len(line) + 1 for line in (name, seq, plus, qual))
    if not payload.endswith(b"\n"):
        global_offset -= 1
    return records


def scan(path, verify_crc=True):
    data = Path(path).read_bytes()
    detected = False
    if len(data) >= 18:
        detected = find_bsize(data[:18])[0] is not None
    if not detected:
        return {
            "path": str(Path(path)), "bgzf": False,
            "has_eof": False, "file_size": len(data), "blocks": [],
            "record_count": 0, "checksum64": "0000000000000000",
        }

    blocks = []
    chunks = []
    uncompressed_starts = []
    offset = 0
    while offset < len(data):
        if len(data) - offset < 18:
            raise ValueError(f"truncated BGZF header at block offset {offset}")
        xlen = struct.unpack_from("<H", data, offset + 10)[0]
        header_end = offset + 12 + xlen
        if header_end > len(data):
            raise ValueError(f"truncated BGZF extra field at block offset {offset}")
        size, relative_header_end = find_bsize(data[offset:header_end])
        if size is None:
            raise ValueError(f"missing BGZF BC field at block offset {offset}")
        if offset + size > len(data):
            raise ValueError(f"truncated BGZF member at block offset {offset}")
        member = data[offset:offset + size]
        if size < relative_header_end + 8:
            raise ValueError(f"invalid BGZF member size at block offset {offset}")
        expected_crc, isize = struct.unpack_from("<II", member, size - 8)
        try:
            inflated = zlib.decompress(member[relative_header_end:size - 8], -15)
        except zlib.error as exc:
            raise ValueError(f"inflate failed at block offset {offset}: {exc}") from exc
        if len(inflated) != isize:
            raise ValueError(f"ISIZE mismatch at block offset {offset}")
        if verify_crc and zlib.crc32(inflated) != expected_crc:
            raise ValueError(f"CRC mismatch at block offset {offset}")
        if member == BGZF_EOF:
            if offset + size != len(data):
                raise ValueError(f"BGZF EOF marker is not final at block offset {offset}")
            break
        uncompressed_starts.append(sum(len(chunk) for chunk in chunks))
        chunks.append(inflated)
        blocks.append({
            "compressed_offset": offset,
            "compressed_size": size,
            "isize": isize,
            "first_record_ordinal": 0,
            "records_starting": 0,
            "first_record_offset": None,
        })
        offset += size

    has_eof = len(data) >= len(BGZF_EOF) and data.endswith(BGZF_EOF)
    payload = b"".join(chunks)
    records = parse_fastq(payload) if payload else []
    for ordinal, (start, _name, _seq, _qual) in enumerate(records):
        block_index = bisect.bisect_right(uncompressed_starts, start) - 1
        block = blocks[block_index]
        if block["records_starting"] == 0:
            block["first_record_ordinal"] = ordinal
            block["first_record_offset"] = start - uncompressed_starts[block_index]
        block["records_starting"] += 1
    cumulative = 0
    for block in blocks:
        if block["records_starting"] == 0:
            block["first_record_ordinal"] = cumulative
        cumulative += block["records_starting"]

    checksum = 0
    for _start, name, seq, qual in records:
        checksum = (checksum + fnv1a(name + b"\0" + seq + b"\0" + qual)) & MASK64
    return {
        "path": str(Path(path)), "bgzf": True, "has_eof": has_eof,
        "file_size": len(data), "blocks": blocks,
        "record_count": len(records), "checksum64": f"{checksum:016x}",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--no-crc", action="store_true")
    args = parser.parse_args()
    try:
        result = scan(args.input, not args.no_crc)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"ERROR: {exc}\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
