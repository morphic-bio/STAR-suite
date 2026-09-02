#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 [--block-bytes N] INPUT.fastq[.gz] OUTPUT.fastq.gz" >&2
}

block_bytes=32768
if [[ "${1:-}" == "--block-bytes" ]]; then
    [[ $# -ge 4 ]] || { usage; exit 2; }
    block_bytes="$2"
    shift 2
fi
[[ $# -eq 2 ]] || { usage; exit 2; }

python3 - "$1" "$2" "$block_bytes" <<'PY'
import gzip
import os
import struct
import sys
import tempfile
import zlib
from pathlib import Path

source = Path(sys.argv[1])
dest = Path(sys.argv[2])
block_bytes = int(sys.argv[3])
if not 1 <= block_bytes <= 65536:
    raise SystemExit("--block-bytes must be in [1, 65536]")

BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)

def member(payload: bytes) -> bytes:
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    body = compressor.compress(payload) + compressor.flush()
    total = 18 + len(body) + 8
    if total > 65536:
        raise ValueError(f"compressed BGZF member is too large: {total} bytes")
    header = (
        b"\x1f\x8b\x08\x04" + b"\0" * 4 + b"\0\xff" +
        struct.pack("<H", 6) + b"BC" + struct.pack("<H", 2) +
        struct.pack("<H", total - 1)
    )
    trailer = struct.pack("<II", zlib.crc32(payload), len(payload) & 0xffffffff)
    return header + body + trailer

opener = gzip.open if source.suffix == ".gz" else open
dest.parent.mkdir(parents=True, exist_ok=True)
fd, tmp_name = tempfile.mkstemp(prefix=dest.name + ".", dir=str(dest.parent))
try:
    with opener(source, "rb") as src, os.fdopen(fd, "wb") as out:
        while True:
            chunk = src.read(block_bytes)
            if not chunk:
                break
            while True:
                try:
                    out.write(member(chunk))
                    break
                except ValueError:
                    if len(chunk) == 1:
                        raise
                    split = len(chunk) // 2
                    out.write(member(chunk[:split]))
                    chunk = chunk[split:]
        out.write(BGZF_EOF)
    os.replace(tmp_name, dest)
except BaseException:
    try:
        os.unlink(tmp_name)
    except FileNotFoundError:
        pass
    raise
PY
