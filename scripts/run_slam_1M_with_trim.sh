#!/usr/bin/env bash
set -euo pipefail

# Stage 4: STAR-SLAM alignment on 1M FASTQ subset with trimming
# Uses blank-derived trim constants

OUT_BASE=${OUT_BASE:-/mnt/pikachu/slam_blank_artifacts_20260201}
FASTQ_DIR=${FASTQ_DIR:-/mnt/pikachu/NW-5-21/SLAM-Seq-1M}
STAR_BIN=${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}
STAR_INDEX=${STAR_INDEX:-/storage/autoindex_110_44/bulk_index}
THREADS=${THREADS:-24}

# SNP mask from Stage 0a
SNP_MASK="$OUT_BASE/mask/snps_from_vcf.bed.gz"

# Trims from QC JSON
QC_JSON="$OUT_BASE/qc/trim_6h.slam_qc.json"
TRIM5=$(python3 -c "import json; print(int(json.load(open('$QC_JSON'))['trim5p']))")
TRIM3=$(python3 -c "import json; print(int(json.load(open('$QC_JSON'))['trim3p']))")

# Output directory for trimmed runs
RUN_DIR="$OUT_BASE/fastq_1M_runs_trimmed"
mkdir -p "$RUN_DIR/logs"

echo "========================================================================"
echo "STAR-SLAM 1M FASTQ Alignment (With Trimming)"
echo "========================================================================"
echo "Date: $(date)"
echo "Output: $RUN_DIR"
echo "FASTQ dir: $FASTQ_DIR"
echo "STAR: $STAR_BIN"
echo "SNP mask: $SNP_MASK"
echo "Trim: 5'=$TRIM5, 3'=$TRIM3"
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
    --slamCompatTrim5p "$TRIM5" \
    --slamCompatTrim3p "$TRIM3" \
    --slamStrandness Sense \
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
echo "STAR-SLAM (trimmed) outputs:"
echo "========================================================================"
ls -la "$RUN_DIR"/*_SlamQuant.out 2>/dev/null || echo "No SlamQuant outputs yet"

echo ""
echo "========================================================================"
echo "Comparing trimmed vs compat mode T→C rates:"
echo "========================================================================"
echo "Sample          Compat(%)  Trimmed(%)  Delta(%)"
for name in "${!SAMPLES[@]}"; do
  compat_file="$OUT_BASE/fastq_1M_runs/${name}_SlamQuant.out.transitions.tsv"
  trimmed_file="$RUN_DIR/${name}_SlamQuant.out.transitions.tsv"
  
  if [[ -f "$compat_file" && -f "$trimmed_file" ]]; then
    compat=$(awk -F'\t' '$1=="T" && $2=="C" {print $4/$3*100}' "$compat_file")
    trimmed=$(awk -F'\t' '$1=="T" && $2=="C" {print $4/$3*100}' "$trimmed_file")
    delta=$(python3 - <<PY
compat = float("$compat")
trimmed = float("$trimmed")
print(f"{trimmed - compat:.4f}")
PY
)
    printf "%-15s %8.4f   %8.4f   %+8.4f\n" "$name" "$compat" "$trimmed" "$delta"
  fi
done

echo ""
echo "DONE. Outputs stored in: $RUN_DIR"
