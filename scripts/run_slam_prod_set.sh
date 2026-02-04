#!/usr/bin/env bash
set -euo pipefail

# Run STAR-SLAM across the full production FASTQ set, using a no4su blank
# for background error-rate estimation and SNP masking.
#
# This script mirrors the batch-run prompt we used for autoslam tests,
# but loops over the full production set.

OUT_BASE=/mnt/pikachu/slam_blank_artifacts_20260201/prod_run_$(date +%Y%m%d)
FASTQ_DIR=/mnt/pikachu/NW-5-21/SLAM-Seq
STAR_BIN=/mnt/pikachu/STAR-suite/core/legacy/source/STAR
STAR_INDEX=/storage/autoindex_110_44/bulk_index
SNP_MASK=/mnt/pikachu/slam_blank_artifacts_20260201/mask/snps_from_vcf.bed.gz
THREADS=24
BLANK_FASTQ=/mnt/pikachu/NW-5-21/SLAM-Seq/ARID1A-no4su_S50_R1_001.fastq.gz

mkdir -p "$OUT_BASE/logs"

if [[ ! -x "$STAR_BIN" ]]; then
  echo "ERROR: STAR binary not found: $STAR_BIN"
  exit 1
fi

if [[ ! -d "$STAR_INDEX" ]]; then
  echo "ERROR: STAR index not found: $STAR_INDEX"
  exit 1
fi

if [[ ! -f "$SNP_MASK" ]]; then
  echo "ERROR: SNP mask not found: $SNP_MASK"
  exit 1
fi

if [[ ! -f "$BLANK_FASTQ" ]]; then
  echo "ERROR: Blank (no4su) FASTQ not found: $BLANK_FASTQ"
  exit 1
fi

sample_name_from_fastq() {
  local fq="$1"
  local base
  base="$(basename "$fq")"
  # strip extensions
  base="${base%.fastq.gz}"
  base="${base%.fq.gz}"
  base="${base%.fastq}"
  base="${base%.fq}"
  # strip read tokens
  base="${base%_R1_001}"
  base="${base%_R2_001}"
  base="${base%_R1}"
  base="${base%_R2}"
  base="${base%.R1}"
  base="${base%.R2}"
  echo "$base"
}

echo "========================================================================"
echo "STAR-SLAM Production Set"
echo "========================================================================"
echo "Date: $(date)"
echo "FASTQ dir: $FASTQ_DIR"
echo "Output base: $OUT_BASE"
echo "Blank (no4su): $BLANK_FASTQ"
echo "STAR: $STAR_BIN"
echo "STAR index: $STAR_INDEX"
echo "SNP mask: $SNP_MASK"
echo "Threads: $THREADS"
echo "========================================================================"

# Collect FASTQs (single-end only)
mapfile -t FASTQS_ALL < <(ls -1 "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | sort)

if [[ ${#FASTQS_ALL[@]} -eq 0 ]]; then
  echo "ERROR: No FASTQ files found in $FASTQ_DIR"
  exit 1
fi

# Ensure blank is first in the list (required for trimScope first + blank seeding)
READS_IN=("$BLANK_FASTQ")
for fq in "${FASTQS_ALL[@]}"; do
  if [[ "$fq" != "$BLANK_FASTQ" ]]; then
    READS_IN+=("$fq")
  fi
done

echo "Running single STAR command with ${#READS_IN[@]} FASTQs..."
LOG="$OUT_BASE/logs/prod_run.log"

"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "${READS_IN[@]}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$OUT_BASE/" \
  --outFileNamePrefixAuto 1 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM MD \
  --outSAMprimaryFlag OneBestScore \
  --emitNoYBAM yes \
  --keepBAM yes \
  --slamQuantMode 1 \
  --slamGrandSlamOut 1 \
  --slamSnpMaskIn "$SNP_MASK" \
  --slamStrandness Sense \
  --trimScope first \
  --trimSource "$BLANK_FASTQ" \
  --slamErrorRateFromBlank 1 \
  > "$LOG" 2>&1

echo ""
echo "DONE. Outputs under: $OUT_BASE"
