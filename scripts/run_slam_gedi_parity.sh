#!/usr/bin/env bash
set -euo pipefail

# Stage 3: GEDI parity comparison on 1M subset
# Compares GEDI vs STAR-SLAM transition rates

OUT_BASE=${OUT_BASE:-/mnt/pikachu/slam_blank_artifacts_20260201}
RUN_DIR="$OUT_BASE/fastq_1M_runs"
GEDI_BIN=${GEDI_BIN:-/mnt/pikachu/STAR-Flex/gedi}
GEDI_GENOME=${GEDI_GENOME:-/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml}
SNP_MASK="$OUT_BASE/mask/snps_from_vcf.bed.gz"
COMPARE_SCRIPT="/mnt/pikachu/STAR-suite/tests/slam/compare_fixture.py"
BEDTOOLS_BIN=${BEDTOOLS_BIN:-bedtools}

# Extract trim and error values from QC JSON (override via env if needed)
QC_JSON="$OUT_BASE/qc/trim_6h.slam_qc.json"
TRIM5_DEFAULT=$(python3 -c "import json; print(int(json.load(open('$QC_JSON'))['trim5p']))")
TRIM3_DEFAULT=$(python3 -c "import json; print(int(json.load(open('$QC_JSON'))['trim3p']))")
SNP_ERR_DEFAULT=$(python3 -c "import json; print(f\"{json.load(open('$QC_JSON')).get('snp_err_used', 0.001):.6f}\")")

TRIM5=${TRIM5_OVERRIDE:-$TRIM5_DEFAULT}
TRIM3=${TRIM3_OVERRIDE:-$TRIM3_DEFAULT}
SNP_ERR=${SNP_ERR_OVERRIDE:-$SNP_ERR_DEFAULT}

GEDI_DIR="$RUN_DIR/gedi"
COMPARE_DIR="$RUN_DIR/compare"
MASKED_DIR="$GEDI_DIR/bam_masked"
mkdir -p "$GEDI_DIR" "$COMPARE_DIR" "$MASKED_DIR" "$RUN_DIR/logs"

echo "========================================================================"
echo "GEDI Parity Comparison (1M subset)"
echo "========================================================================"
echo "Date: $(date)"
echo "Run dir: $RUN_DIR"
echo "GEDI: $GEDI_BIN"
echo "Genome: $GEDI_GENOME"
echo "Trim: 5'=$TRIM5, 3'=$TRIM3"
echo "SNP error: $SNP_ERR"
echo "========================================================================"

# Verify GEDI
if [[ ! -x "$GEDI_BIN" ]]; then
  echo "WARNING: GEDI not found: $GEDI_BIN"
  echo "Skipping GEDI runs..."
  exit 0
fi
if [[ ! -f "$GEDI_GENOME" ]]; then
  echo "WARNING: GEDI genome not found: $GEDI_GENOME"
  echo "Skipping GEDI runs..."
  exit 0
fi
if ! command -v "$BEDTOOLS_BIN" >/dev/null 2>&1; then
  echo "ERROR: bedtools not found; required to SNP-mask BAMs for GEDI parity"
  exit 1
fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not found; required to index SNP-masked BAMs for GEDI"
  exit 1
fi

# Samples to process
SAMPLES=("ARID1A_0h_1" "ARID1A_24h_1" "ARID1A_6h_1" "ARID1A_no4su")

mask_bam_for_gedi() {
  local bam=$1
  local name=$2
  local masked="$MASKED_DIR/${name}.snpmask.bam"

  if [[ ! -s "$masked" ]]; then
    echo "SNP-masking BAM for GEDI: $bam -> $masked" >&2
    "$BEDTOOLS_BIN" intersect -v -abam "$bam" -b "$SNP_MASK" > "$masked"
  fi
  if [[ ! -f "${masked}.bai" ]]; then
    samtools index "$masked"
  fi
  echo "$masked"
}

# Run GEDI on each sample
echo ""
echo "Running GEDI on samples..."
for name in "${SAMPLES[@]}"; do
  bam="$RUN_DIR/${name}_Aligned.sortedByCoord.out.bam"
  
  if [[ ! -f "$bam" ]]; then
    echo "WARNING: BAM not found: $bam"
    continue
  fi
  
  if [[ -s "$GEDI_DIR/${name}.tsv.gz" ]]; then
    echo "✓ GEDI $name already processed, skipping"
    continue
  fi

  masked_bam=$(mask_bam_for_gedi "$bam" "$name")
  
  echo "Running GEDI on $name..."
  "$GEDI_BIN" -e Slam \
    -reads "$masked_bam" \
    -genomic "$GEDI_GENOME" \
    -prefix "$GEDI_DIR/${name}" \
    -trim5p "$TRIM5" \
    -trim3p "$TRIM3" \
    -strandness Sense \
    -err "$SNP_ERR" \
    -D \
    > "$RUN_DIR/logs/gedi_${name}.log" 2>&1
  
  echo "✓ GEDI $name completed"
done

# Compare GEDI vs STAR-SLAM
echo ""
echo "Comparing GEDI vs STAR-SLAM..."
for name in "${SAMPLES[@]}"; do
  ref="$GEDI_DIR/${name}.tsv.gz"
  test="$RUN_DIR/${name}_SlamQuant.out"
  output="$COMPARE_DIR/compare_${name}.txt"
  
  if [[ ! -f "$ref" || ! -f "$test" ]]; then
    echo "WARNING: Missing files for $name (ref=$ref, test=$test)"
    continue
  fi

  # Skip compare if GEDI output lacks MAP column (can happen for blank/no4su)
  if ! python3 - "$ref" <<'PY'
import gzip, sys
path = sys.argv[1]
opener = gzip.open if path.endswith(".gz") else open
with opener(path, "rt") as f:
    header = f.readline()
sys.exit(0 if "MAP" in header else 1)
PY
  then
    echo "WARNING: GEDI output missing MAP column for $name; skipping parity compare"
    continue
  fi
  
  echo "Comparing $name..."
  python3 "$COMPARE_SCRIPT" \
    --reference "$ref" \
    --test "$test" \
    --corr-min 0.9 \
    --thresholds 20,50,100 \
    | tee "$output"
  echo ""
done

# Summary
echo "========================================================================"
echo "Parity Summary"
echo "========================================================================"
for f in "$COMPARE_DIR"/compare_*.txt; do
  if [[ -f "$f" ]]; then
    name=$(basename "$f" .txt | sed 's/compare_//')
    result=$(grep -E "^(PASS|FAIL)" "$f" | head -1)
    pearson=$(grep "NTR correlation (Pearson)" "$f" | head -1 | sed 's/.*: //' | cut -d' ' -f1)
    echo "$name: $result (Pearson=$pearson)"
  fi
done

echo ""
echo "DONE. Outputs in: $GEDI_DIR, $COMPARE_DIR"
exit 0
