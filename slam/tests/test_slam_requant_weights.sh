#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

STAR_BIN="${STAR_BIN:-$ROOT_DIR/source/STAR}"
REQ_BIN="${REQ_BIN:-$ROOT_DIR/tools/slam_requant/slam_requant}"
FASTQ="${FASTQ:-$ROOT_DIR/test/fixtures/slam/raw/slam_100000_reads_SRR32576116.fastq.gz}"
STAR_INDEX="${STAR_INDEX:-$ROOT_DIR/test/fixtures/slam/ref/star_index}"
SNPS_BED="${SNPS_BED:-$ROOT_DIR/test/fixtures/slam/ref/snps.bed}"

WORK="${WORK:-$ROOT_DIR/test/tmp_slam_requant_weights}"
OUT_PREFIX="${OUT_PREFIX:-$WORK/fixture_}"
DUMP_PATH="${DUMP_PATH:-$WORK/fixture.dump.bin}"
WEIGHT_PATH="${WEIGHT_PATH:-$WORK/fixture.weights.bin}"
WEIGHT_SHUFFLED="${WEIGHT_SHUFFLED:-$WORK/fixture.weights.reversed.bin}"

COMPARE="${COMPARE:-$SCRIPT_DIR/slam/compare_star_outputs.py}"

mkdir -p "$WORK"

if [[ ! -x "$STAR_BIN" || ! -x "$REQ_BIN" ]]; then
  echo "Missing STAR or slam_requant binary"
  exit 1
fi

echo "=== STAR: build dump + weights ==="
"$STAR_BIN" \
  --runThreadN 4 \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$FASTQ" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUT_PREFIX" \
  --outSAMtype None \
  --clip3pAdapterSeq AGATCGGAAGAG \
  --clip3pAdapterMMp 0.1 \
  --slamQuantMode 1 \
  --slamSnpMaskIn "$SNPS_BED" \
  --slamDumpBinary "$DUMP_PATH" \
  --slamDumpWeights "$WEIGHT_PATH" \
  > "${OUT_PREFIX}slam.log" 2>&1

echo "=== Rewrite weights in reverse order ==="
python3 - "$WEIGHT_PATH" "$WEIGHT_SHUFFLED" <<'PY'
import struct
import sys

in_path = sys.argv[1]
out_path = sys.argv[2]

with open(in_path, "rb") as f:
    header = f.read(8 + 4 + 4 + 8 + 4)
    if len(header) != 28:
        raise SystemExit("bad header")
    magic = header[:8]
    if magic != b"SLAMWGT1":
        raise SystemExit("bad magic")
    version, flags, n_reads, weight_mode = struct.unpack("<IIQI", header[8:])
    recs = []
    rec_fmt = "<QQd"
    rec_size = struct.calcsize(rec_fmt)
    for _ in range(n_reads):
        data = f.read(rec_size)
        if len(data) != rec_size:
            raise SystemExit("short record")
        recs.append(data)

with open(out_path, "wb") as f:
    f.write(header)
    for rec in reversed(recs):
        f.write(rec)
PY

echo "=== slam_requant: replay dump (key match) ==="
"$REQ_BIN" \
  --dump "$DUMP_PATH" \
  --out "$WORK/requant_key_" \
  --slamSnpMaskIn "$SNPS_BED" \
  --slamWeightFile "$WEIGHT_SHUFFLED" \
  --slamWeightMatch key

echo "=== Compare STAR vs slam_requant (keyed weights) ==="
python3 "$COMPARE" \
  --reference "${OUT_PREFIX}SlamQuant.out" \
  --test "$WORK/requant_key_SlamQuant.out"
