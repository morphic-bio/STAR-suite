#!/usr/bin/env bash
set -euo pipefail

# Stage 2: STAR-SLAM alignment on 1M FASTQ subset
# Uses SNP masking but no trimming (compat mode)

OUT_BASE=${OUT_BASE:-/mnt/pikachu/slam_blank_artifacts_20260201}
FASTQ_DIR=${FASTQ_DIR:-/mnt/pikachu/NW-5-21/SLAM-Seq-1M}
STAR_BIN=${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}
STAR_INDEX=${STAR_INDEX:-/storage/autoindex_110_44/bulk_index}
THREADS=${THREADS:-24}
STAR_EXTRA_ARGS=${STAR_EXTRA_ARGS:-}

# SNP mask from Stage 0a
SNP_MASK="$OUT_BASE/mask/snps_from_vcf.bed.gz"

# Output directory for 1M runs
RUN_DIR="$OUT_BASE/fastq_1M_runs"
mkdir -p "$RUN_DIR/logs"

echo "========================================================================"
echo "STAR-SLAM 1M FASTQ Alignment (Compat Mode - No Trimming)"
echo "========================================================================"
echo "Date: $(date)"
echo "Output: $RUN_DIR"
echo "FASTQ dir: $FASTQ_DIR"
echo "STAR: $STAR_BIN"
echo "SNP mask: $SNP_MASK"
echo "Threads: $THREADS"
echo "========================================================================"

# Verify inputs
if [[ ! -x "$STAR_BIN" ]]; then
  echo "ERROR: STAR binary not found: $STAR_BIN"
  exit 1
fi
if [[ ! -f "$SNP_MASK" ]]; then
  echo "ERROR: SNP mask not found: $SNP_MASK"
  exit 1
fi

# Define samples to process (0h-1, 24h-1, no4su for ARID1A group)
declare -A SAMPLES=(
  ["ARID1A_0h_1"]="$FASTQ_DIR/ARID1A-0h-1_S40_R1_001.fastq.gz"
  ["ARID1A_24h_1"]="$FASTQ_DIR/ARID1A-24h-1_S46_R1_001.fastq.gz"
  ["ARID1A_6h_1"]="$FASTQ_DIR/ARID1A-6h-1_S43_R1_001.fastq.gz"
  ["ARID1A_no4su"]="$FASTQ_DIR/ARID1A-no4su_S50_R1_001.fastq.gz"
)

run_star_slam() {
  local name=$1
  local fastq=$2
  local prefix="$RUN_DIR/${name}_"
  local log="$RUN_DIR/logs/${name}.log"
  
  if [[ -s "${prefix}SlamQuant.out" && -s "${prefix}Aligned.sortedByCoord.out.bam" ]]; then
    echo "✓ $name already processed, skipping"
    return 0
  fi
  
  echo "Processing $name..."
  "$STAR_BIN" \
    --runThreadN "$THREADS" \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$fastq" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$prefix" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM MD \
    --outSAMprimaryFlag OneBestScore \
    --slamQuantMode 1 \
    --slamGrandSlamOut 1 \
    --slamSnpMaskIn "$SNP_MASK" \
    --slamStrandness Sense \
    $STAR_EXTRA_ARGS \
    > "$log" 2>&1
  
  echo "✓ $name completed"
}

echo ""
echo "Processing samples..."
for name in "${!SAMPLES[@]}"; do
  fastq="${SAMPLES[$name]}"
  if [[ ! -f "$fastq" ]]; then
    echo "WARNING: FASTQ not found: $fastq"
    continue
  fi
  run_star_slam "$name" "$fastq"
done

echo ""
echo "========================================================================"
echo "STAR-SLAM outputs:"
echo "========================================================================"
ls -la "$RUN_DIR"/*.out 2>/dev/null || echo "No .out files yet"

echo ""
echo "========================================================================"
echo "Transitions files:"
echo "========================================================================"
for f in "$RUN_DIR"/*_SlamQuant.out.transitions.tsv; do
  if [[ -f "$f" ]]; then
    echo "=== $(basename "$f") ==="
    head -5 "$f"
    echo ""
  fi
done

echo ""
echo "DONE. Outputs stored in: $RUN_DIR"
