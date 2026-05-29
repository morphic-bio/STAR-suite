#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -f "$SCRIPT_DIR/external_fixtures_env.sh" ]]; then
  # shellcheck source=/dev/null
  source "$SCRIPT_DIR/external_fixtures_env.sh"
fi

DEFAULT_REQ_BIN="$ROOT_DIR/slam/tools/slam_requant/slam_requant"
if [[ ! -x "$DEFAULT_REQ_BIN" && -x "$ROOT_DIR/tools/slam_requant/slam_requant" ]]; then
  DEFAULT_REQ_BIN="$ROOT_DIR/tools/slam_requant/slam_requant"
fi

STAR_BIN="${STAR_BIN:-$ROOT_DIR/core/legacy/source/STAR}"
REQ_BIN="${REQ_BIN:-$DEFAULT_REQ_BIN}"
FASTQ="${FASTQ:-${SLAM_FIXTURE_FASTQ:-$ROOT_DIR/test/fixtures/slam/raw/slam_100000_reads_SRR32576116.fastq.gz}}"
STAR_INDEX="${STAR_INDEX:-${SLAM_FIXTURE_STAR_INDEX:-$ROOT_DIR/test/fixtures/slam/ref/star_index}}"
SNPS_BED="${SNPS_BED:-${SLAM_FIXTURE_SNPS_BED:-$ROOT_DIR/test/fixtures/slam/ref/snps.bed}}"
EXACT_MIN_PEARSON="${EXACT_MIN_PEARSON:-0.999999}"
ALIGN_MIN_PEARSON="${ALIGN_MIN_PEARSON:-0.98}"
EXACT_MAX_ABS_DELTA="${EXACT_MAX_ABS_DELTA:-1e-6}"

WORK="${WORK:-$ROOT_DIR/test/tmp_slam_requant}"
OUT_PREFIX="${OUT_PREFIX:-$WORK/fixture_}"
DUMP_PATH="${DUMP_PATH:-$WORK/fixture.dump.bin}"
WEIGHT_PATH="${WEIGHT_PATH:-$WORK/fixture.weights.bin}"

COMPARE="${COMPARE:-$SCRIPT_DIR/slam/compare_star_outputs.py}"

mkdir -p "$WORK"

if [[ ! -x "$STAR_BIN" || ! -x "$REQ_BIN" ]]; then
  echo "Missing STAR or slam_requant binary"
  exit 1
fi

echo "=== STAR: build dump + SlamQuant ==="
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

echo "=== slam_requant: replay dump (dump weights) ==="
"$REQ_BIN" \
  --dump "$DUMP_PATH" \
  --out "$WORK/requant_dump_" \
  --slamSnpMaskIn "$SNPS_BED"

echo "=== slam_requant: replay dump (alignments weights) ==="
"$REQ_BIN" \
  --dump "$DUMP_PATH" \
  --out "$WORK/requant_align_" \
  --slamSnpMaskIn "$SNPS_BED" \
  --slamWeightMode alignments

echo "=== slam_requant: replay dump (weight sidecar) ==="
"$REQ_BIN" \
  --dump "$DUMP_PATH" \
  --out "$WORK/requant_weight_" \
  --slamSnpMaskIn "$SNPS_BED" \
  --slamWeightFile "$WEIGHT_PATH"

echo "=== Compare STAR vs slam_requant (dump weights) ==="
python3 "$COMPARE" \
  --reference "${OUT_PREFIX}SlamQuant.out" \
  --test "$WORK/requant_dump_SlamQuant.out" \
  --min-pearson "$EXACT_MIN_PEARSON" \
  --max-abs-delta "$EXACT_MAX_ABS_DELTA" \
  --require-same-genes

echo "=== Compare STAR vs slam_requant (alignments weights) ==="
python3 "$COMPARE" \
  --reference "${OUT_PREFIX}SlamQuant.out" \
  --test "$WORK/requant_align_SlamQuant.out" \
  --min-pearson "$ALIGN_MIN_PEARSON"

echo "=== Compare STAR vs slam_requant (weight sidecar) ==="
python3 "$COMPARE" \
  --reference "${OUT_PREFIX}SlamQuant.out" \
  --test "$WORK/requant_weight_SlamQuant.out" \
  --min-pearson "$EXACT_MIN_PEARSON" \
  --max-abs-delta "$EXACT_MAX_ABS_DELTA" \
  --require-same-genes
