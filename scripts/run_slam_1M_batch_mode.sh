#!/usr/bin/env bash
set -euo pipefail

# Batch SLAM-seq run on 1M downsampled FASTQs (blank-first)
# Uses batch mode to reuse the loaded genome and propagate blank error rate.
# Adds TranscriptVB quantification (Salmon-style).

OUT_BASE=${OUT_BASE:-/mnt/pikachu/slam_batch_1M_$(date +%Y%m%d_%H%M%S)}
THREADS=${THREADS:-24}

STAR_BIN=${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}
STAR_INDEX=${STAR_INDEX:-/storage/autoindex_110_44/bulk_index}
FASTQ_DIR=${FASTQ_DIR:-/mnt/pikachu/NW-5-21/SLAM-Seq-1M}
SNP_MASK=${SNP_MASK:-/mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz}

# Default: run all 1M R1 FASTQs, blank first
# If FASTQ_LIST_OVERRIDE is set, it must be a comma-separated list with blank first.
USE_ALL_FASTQS=${USE_ALL_FASTQS:-1}
FASTQ_LIST_OVERRIDE=${FASTQ_LIST_OVERRIDE:-}

# Optional explicit subset (used if USE_ALL_FASTQS=0)
FASTQ_BLANK=${FASTQ_BLANK:-$FASTQ_DIR/ARID1A-no4su_S50_R1_001.fastq.gz}
FASTQ_0H=${FASTQ_0H:-$FASTQ_DIR/ARID1A-0h-1_S40_R1_001.fastq.gz}
FASTQ_6H=${FASTQ_6H:-$FASTQ_DIR/ARID1A-6h-1_S43_R1_001.fastq.gz}
FASTQ_24H=${FASTQ_24H:-$FASTQ_DIR/ARID1A-24h-1_S46_R1_001.fastq.gz}

FASTQS=()
if [[ -n "$FASTQ_LIST_OVERRIDE" ]]; then
  IFS=, read -r -a FASTQS <<< "$FASTQ_LIST_OVERRIDE"
elif [[ "$USE_ALL_FASTQS" == "1" ]]; then
  # Find blank first (prefer FASTQ_BLANK if it exists)
  if [[ -f "$FASTQ_BLANK" ]]; then
    BLANK_PATH="$FASTQ_BLANK"
  else
    BLANK_PATH=$(find "$FASTQ_DIR" -maxdepth 1 -type f -name '*no4su*_R1_*.fastq.gz' | sort | head -n 1)
  fi
  if [[ -z "${BLANK_PATH:-}" || ! -f "$BLANK_PATH" ]]; then
    echo "ERROR: blank FASTQ not found (no4su). Set FASTQ_BLANK or place a no4su file in $FASTQ_DIR." >&2
    exit 1
  fi

  mapfile -t FASTQS_ALL < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name '*_R1_*.fastq.gz' | sort)
  if [[ ${#FASTQS_ALL[@]} -eq 0 ]]; then
    echo "ERROR: no R1 FASTQs found in $FASTQ_DIR" >&2
    exit 1
  fi

  # Remove blank from the list, then prepend it
  FASTQS=("$BLANK_PATH")
  for fq in "${FASTQS_ALL[@]}"; do
    if [[ "$fq" != "$BLANK_PATH" ]]; then
      FASTQS+=("$fq")
    fi
  done
else
  FASTQS=("$FASTQ_BLANK" "$FASTQ_0H" "$FASTQ_6H" "$FASTQ_24H")
fi

QUANT_MODE=${QUANT_MODE:-TranscriptVB}
# For single-end SLAM data, auto-detect (A) fails; default to U (unstranded SE)
QUANT_VB_LIBTYPE=${QUANT_VB_LIBTYPE:-U}
SLAM_DUMP_BINARY=${SLAM_DUMP_BINARY:-1}
SLAM_DUMP_WEIGHTS=${SLAM_DUMP_WEIGHTS:-1}

mkdir -p "$OUT_BASE"

# Basic checks
if [[ ! -x "$STAR_BIN" ]]; then
  echo "ERROR: STAR_BIN not executable: $STAR_BIN" >&2
  exit 1
fi
if [[ ! -d "$STAR_INDEX" ]]; then
  echo "ERROR: STAR_INDEX not found: $STAR_INDEX" >&2
  exit 1
fi
if [[ ! -f "$SNP_MASK" ]]; then
  echo "ERROR: SNP mask not found: $SNP_MASK" >&2
  exit 1
fi
for fq in "${FASTQS[@]}"; do
  if [[ ! -f "$fq" ]]; then
    echo "ERROR: FASTQ missing: $fq" >&2
    exit 1
  fi
 done

# Join FASTQs into a single comma-separated list (blank must be first)
FASTQ_LIST=$(IFS=,; echo "${FASTQS[*]}")

cat <<INFO
====================================================================
SLAM 1M Batch Mode
--------------------------------------------------------------------
Output base: $OUT_BASE
STAR:        $STAR_BIN
Index:       $STAR_INDEX
SNP mask:    $SNP_MASK
Threads:     $THREADS
QuantMode:   $QUANT_MODE
DumpBinary:  $SLAM_DUMP_BINARY
DumpWeights: $SLAM_DUMP_WEIGHTS
FASTQ count: ${#FASTQS[@]}
FASTQs:      $FASTQ_LIST
====================================================================
INFO

QB_ARGS=()
if [[ "$QUANT_MODE" == "TranscriptVB" ]]; then
  QB_ARGS=(--quantVBLibType "$QUANT_VB_LIBTYPE")
fi

"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --runMode alignReads \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$FASTQ_LIST" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUT_BASE/" \
  --outFileNamePrefixAuto 1 \
  --outSAMtype None \
  --quantMode "$QUANT_MODE" \
  "${QB_ARGS[@]}" \
  --emitNoYBAM yes \
  --emitYNoYFastq yes \
  --slamQuantMode 1 \
  --slamGrandSlamOut 1 \
  --slamSnpMaskIn "$SNP_MASK" \
  --slamStrandness Sense \
  --slamErrorRateFromBlank 1 \
  --slamDumpBinary "$SLAM_DUMP_BINARY" \
  --slamDumpWeights "$SLAM_DUMP_WEIGHTS" \
  --autoTrim - \
  --slamCompatTrim5p 0 \
  --slamCompatTrim3p 0 \
  --batchMode 1

echo "DONE. Outputs in: $OUT_BASE"
