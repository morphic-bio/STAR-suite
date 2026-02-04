#!/usr/bin/env bash
set -euo pipefail

# Batch SLAM-seq run on full production FASTQs (blank-first)
# Uses batch mode to reuse the loaded genome and propagate blank error rate.

OUT_BASE=${OUT_BASE:-/mnt/pikachu/NW-5-21/SLAM-Seq-processed}
THREADS=${THREADS:-24}

STAR_BIN=${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}
STAR_INDEX=${STAR_INDEX:-/storage/autoindex_110_44/bulk_index}
FASTQ_DIR=${FASTQ_DIR:-/mnt/pikachu/NW-5-21/SLAM-Seq}
SNP_MASK=${SNP_MASK:-/mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz}

# If FASTQ_LIST_OVERRIDE is set, it must be a comma-separated list with blank first.
FASTQ_LIST_OVERRIDE=${FASTQ_LIST_OVERRIDE:-}

# Optional explicit blank override
FASTQ_BLANK=${FASTQ_BLANK:-}

QUANT_MODE=${QUANT_MODE:-TranscriptVB}
# For single-end SLAM data, auto-detect (A) can fail; default to U (unstranded SE)
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
if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "ERROR: FASTQ_DIR not found: $FASTQ_DIR" >&2
  exit 1
fi

FASTQS=()
if [[ -n "$FASTQ_LIST_OVERRIDE" ]]; then
  IFS=, read -r -a FASTQS <<< "$FASTQ_LIST_OVERRIDE"
else
  # Determine blank (no4su) first
  if [[ -n "$FASTQ_BLANK" ]]; then
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
SLAM Production Batch Mode
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
  --outSAMtype BAM Unsorted \
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
