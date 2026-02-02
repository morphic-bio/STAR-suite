#!/usr/bin/env bash
set -euo pipefail

# Stage 0a: SLAM-seq blank + time-course artifact generation
# - Builds SNP mask from VCF (first run)
# - Runs STAR-SLAM (internal) on 0h/6h/24h/no4su
# - Runs external STAR-SLAM on 0h/6h/24h/no4su (if available)
# - Runs GEDI/GRAND-SLAM on internal BAMs (if available)
# - Compares STAR vs GEDI correlations (Pearson >= 0.9)

OUT_BASE=${OUT_BASE:-/mnt/pikachu/slam_blank_artifacts_$(date +%Y%m%d)}
THREADS=${THREADS:-16}

STAR_BIN=${STAR_BIN:-/mnt/pikachu/STAR-suite/core/legacy/source/STAR}
STAR_EXT_BIN=${STAR_EXT_BIN:-/mnt/pikachu/STAR-Flex/source/STAR}
STAR_INDEX=${STAR_INDEX:-/storage/autoindex_110_44/bulk_index}

VCF_IN=${VCF_IN:-/storage/HPSI0114i-kolf_2.wgs.gatk.haplotype_caller.20170425.genotypes.vcf.gz}
VCF_SAMPLE=${VCF_SAMPLE:-HPSI0114i-kolf_2}

GEDI_BIN=${GEDI_BIN:-/mnt/pikachu/STAR-Flex/gedi}
GEDI_GENOME=${GEDI_GENOME:-/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml}
BEDTOOLS_BIN=${BEDTOOLS_BIN:-bedtools}

FASTQ_0H=${FASTQ_0H:-/mnt/pikachu/NW-5-21/SLAM-Seq/ARID1A-0h-1_S40_R1_001.fastq.gz}
FASTQ_6H=${FASTQ_6H:-/mnt/pikachu/NW-5-21/SLAM-Seq/ARID1A-6h-1_S43_R1_001.fastq.gz}
FASTQ_24H=${FASTQ_24H:-/mnt/pikachu/NW-5-21/SLAM-Seq/ARID1A-24h-1_S46_R1_001.fastq.gz}
FASTQ_BLANK=${FASTQ_BLANK:-/mnt/pikachu/NW-5-21/SLAM-Seq/ARID1A-no4su_S50_R1_001.fastq.gz}

COMPARE_SCRIPT=${COMPARE_SCRIPT:-/mnt/pikachu/STAR-suite/tests/slam/compare_fixture.py}
CORR_MIN=${CORR_MIN:-0.9}

mkdir -p "$OUT_BASE"/{logs,mask,star_internal,star_external,gedi,compare,qc}
mkdir -p "$OUT_BASE/gedi/bam_masked"
OUT_BASE=$(realpath "$OUT_BASE")

LOG_FILE="$OUT_BASE/logs/run.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================================================"
echo "SLAM Blank/Time-course Artifact Generation"
echo "========================================================================"
echo "Date: $(date)"
echo "Output base: $OUT_BASE"
echo "STAR (internal): $STAR_BIN"
echo "STAR (external): $STAR_EXT_BIN"
echo "STAR index: $STAR_INDEX"
echo "VCF: $VCF_IN (sample: $VCF_SAMPLE)"
echo "FASTQ 0h: $FASTQ_0H"
echo "FASTQ 6h: $FASTQ_6H"
echo "FASTQ 24h: $FASTQ_24H"
echo "FASTQ blank(no4su): $FASTQ_BLANK"
echo "GEDI: $GEDI_BIN"
echo "GEDI genome: $GEDI_GENOME"
echo "========================================================================"

if [[ ! -x "$STAR_BIN" ]]; then
  echo "ERROR: STAR_BIN not executable: $STAR_BIN"
  exit 1
fi
if [[ ! -d "$STAR_INDEX" ]]; then
  echo "ERROR: STAR_INDEX not found: $STAR_INDEX"
  exit 1
fi
for fq in "$FASTQ_0H" "$FASTQ_6H" "$FASTQ_24H" "$FASTQ_BLANK"; do
  if [[ ! -f "$fq" ]]; then
    echo "ERROR: FASTQ missing: $fq"
    exit 1
  fi
done
if [[ ! -f "$VCF_IN" ]]; then
  echo "ERROR: VCF missing: $VCF_IN"
  exit 1
fi

MASK_BED="$OUT_BASE/mask/snps_from_vcf.bed.gz"
MASK_SUMMARY="$OUT_BASE/mask/snps_from_vcf_summary.tsv"

QC_PREFIX="$OUT_BASE/qc/trim_6h"

run_star_internal() {
  local label=$1
  local fastq=$2
  local prefix="$OUT_BASE/star_internal/${label}_"
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
    --slamSnpMaskIn "$MASK_BED" \
    --autoTrim - \
    --slamCompatTrim5p "$TRIM5" \
    --slamCompatTrim3p "$TRIM3" \
    --slamStrandness Sense \
    > "$OUT_BASE/logs/star_internal_${label}.log" 2>&1
}

run_star_external() {
  local label=$1
  local fastq=$2
  local prefix="$OUT_BASE/star_external/${label}_"
  "$STAR_EXT_BIN" \
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
    --slamSnpMaskIn "$MASK_BED" \
    --autoTrim - \
    --slamCompatTrim5p "$TRIM5" \
    --slamCompatTrim3p "$TRIM3" \
    --slamStrandness Sense \
    > "$OUT_BASE/logs/star_external_${label}.log" 2>&1
}

ensure_bam_index() {
  local bam="$1"
  if [[ -f "${bam}.bai" ]]; then
    return 0
  fi
  if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not found; required to index BAM for GEDI"
    exit 1
  fi
  samtools index "$bam"
}

echo "[1/6] Running STAR (internal) on 6h with VCF mask + autoTrim..."
if [[ -s "$OUT_BASE/star_internal/6h_SlamQuant.out" && -s "$OUT_BASE/star_internal/6h_Aligned.sortedByCoord.out.bam" \
      && -s "$QC_PREFIX.slam_qc.json" && -s "$MASK_BED" ]]; then
  echo "✓ STAR 6h + VCF mask already present, skipping step 1"
else
  "$STAR_BIN" \
    --runThreadN "$THREADS" \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$FASTQ_6H" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$OUT_BASE/star_internal/6h_" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM MD \
    --slamQuantMode 1 \
    --slamGrandSlamOut 1 \
    --slamSnpMaskVcfIn "$VCF_IN" \
    --slamSnpMaskVcfSample "$VCF_SAMPLE" \
    --slamSnpMaskBedOut "$MASK_BED" \
    --slamSnpMaskSummaryOut "$MASK_SUMMARY" \
    --autoTrim variance \
    --trimSource "$FASTQ_6H" \
    --autoTrimDetectionReads 100000 \
    --slamQcReport "$QC_PREFIX" \
    --slamStrandness Sense \
    > "$OUT_BASE/logs/star_internal_6h_trim.log" 2>&1
fi

if [[ ! -f "$MASK_BED" ]]; then
  echo "ERROR: SNP BED not created: $MASK_BED"
  exit 1
fi
if [[ ! -f "$QC_PREFIX.slam_qc.json" ]]; then
  echo "ERROR: QC JSON not created: $QC_PREFIX.slam_qc.json"
  exit 1
fi

echo "[2/6] Parsing trims + SNP error from QC JSON..."
TRIM5=$(python3 - "$QC_PREFIX.slam_qc.json" << 'PY'
import json, sys
with open(sys.argv[1]) as f:
    data = json.load(f)
print(int(data.get("trim5p", 0)))
PY
)
TRIM3=$(python3 - "$QC_PREFIX.slam_qc.json" << 'PY'
import json, sys
with open(sys.argv[1]) as f:
    data = json.load(f)
print(int(data.get("trim3p", 0)))
PY
)
SNP_ERR=$(python3 - "$QC_PREFIX.slam_qc.json" << 'PY'
import json, sys
with open(sys.argv[1]) as f:
    data = json.load(f)
print(f"{data.get('snp_err_used', 0.001):.6f}")
PY
)
echo "  trim5p=$TRIM5 trim3p=$TRIM3 snp_err=$SNP_ERR"

echo "[3/6] Running STAR (internal) on 0h/24h/no4su with fixed trims..."
for label in 0h 24h blank; do
  out_bam="$OUT_BASE/star_internal/${label}_Aligned.sortedByCoord.out.bam"
  out_slam="$OUT_BASE/star_internal/${label}_SlamQuant.out"
  if [[ -s "$out_bam" && -s "$out_slam" ]]; then
    echo "✓ STAR internal $label already present, skipping"
  else
    if [[ "$label" == "0h" ]]; then
      run_star_internal "0h" "$FASTQ_0H"
    elif [[ "$label" == "24h" ]]; then
      run_star_internal "24h" "$FASTQ_24H"
    else
      run_star_internal "blank" "$FASTQ_BLANK"
    fi
  fi
done

echo "[4/6] Running STAR (external) on 0h/6h/24h/no4su (if available)..."
HAS_SLAM_EXT=0
if [[ -x "$STAR_EXT_BIN" ]]; then
  if "$STAR_EXT_BIN" --help 2>/dev/null | rg -qi "\bslam"; then
    HAS_SLAM_EXT=1
  else
    echo "WARNING: external STAR lacks SLAM support; skipping: $STAR_EXT_BIN"
  fi
else
  echo "WARNING: external STAR not found/executable: $STAR_EXT_BIN"
fi

if [[ "$HAS_SLAM_EXT" == "1" ]]; then
  for label in 0h 6h 24h blank; do
    out_bam="$OUT_BASE/star_external/${label}_Aligned.sortedByCoord.out.bam"
    out_slam="$OUT_BASE/star_external/${label}_SlamQuant.out"
    if [[ -s "$out_bam" && -s "$out_slam" ]]; then
      echo "✓ STAR external $label already present, skipping"
    else
      if [[ "$label" == "0h" ]]; then
        run_star_external "0h" "$FASTQ_0H"
      elif [[ "$label" == "6h" ]]; then
        run_star_external "6h" "$FASTQ_6H"
      elif [[ "$label" == "24h" ]]; then
        run_star_external "24h" "$FASTQ_24H"
      else
        run_star_external "blank" "$FASTQ_BLANK"
      fi
    fi
  done
fi

echo "[5/6] Running GEDI/GRAND-SLAM (if available)..."
RUN_GEDI=1
if [[ ! -x "$GEDI_BIN" ]]; then
  echo "WARNING: GEDI binary not found: $GEDI_BIN"
  RUN_GEDI=0
fi
if [[ ! -f "$GEDI_GENOME" ]]; then
  echo "WARNING: GEDI genome not found: $GEDI_GENOME"
  RUN_GEDI=0
fi

if [[ "$RUN_GEDI" == "1" ]]; then
  if ! command -v "$BEDTOOLS_BIN" >/dev/null 2>&1; then
    echo "ERROR: bedtools not found; required to SNP-mask BAMs for GEDI parity"
    exit 1
  fi
  if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not found; required to index SNP-masked BAMs for GEDI"
    exit 1
  fi

  mask_bam_for_gedi() {
    local bam=$1
    local label=$2
    local masked="$OUT_BASE/gedi/bam_masked/${label}.snpmask.bam"

    if [[ ! -s "$masked" ]]; then
      echo "SNP-masking BAM for GEDI: $bam -> $masked" >&2
      "$BEDTOOLS_BIN" intersect -v -abam "$bam" -b "$MASK_BED" > "$masked"
    fi
    if [[ ! -f "${masked}.bai" ]]; then
      samtools index "$masked"
    fi
    echo "$masked"
  }

  for label in 0h 6h 24h blank; do
    bam="$OUT_BASE/star_internal/${label}_Aligned.sortedByCoord.out.bam"
    if [[ ! -f "$bam" ]]; then
      echo "ERROR: Missing BAM for GEDI: $bam"
      exit 1
    fi
    if [[ -s "$OUT_BASE/gedi/${label}.tsv.gz" ]]; then
      echo "✓ GEDI $label already present, skipping"
    else
      masked_bam=$(mask_bam_for_gedi "$bam" "$label")
      "$GEDI_BIN" -e Slam \
        -reads "$masked_bam" \
        -genomic "$GEDI_GENOME" \
        -prefix "$OUT_BASE/gedi/${label}" \
        -trim5p "$TRIM5" \
        -trim3p "$TRIM3" \
        -strandness Sense \
        -err "$SNP_ERR" \
        -D \
        > "$OUT_BASE/logs/gedi_${label}.log" 2>&1
    fi
  done
fi

echo "[6/6] Comparing STAR vs GEDI correlations (Pearson >= $CORR_MIN)..."
if [[ -f "$COMPARE_SCRIPT" && "$RUN_GEDI" == "1" ]]; then
  for label in 0h 6h 24h; do
    ref="$OUT_BASE/gedi/${label}.tsv.gz"
    test_int="$OUT_BASE/star_internal/${label}_SlamQuant.out"
    test_ext="$OUT_BASE/star_external/${label}_SlamQuant.out"
    if [[ -f "$ref" && -f "$test_int" ]]; then
      python3 "$COMPARE_SCRIPT" \
        --reference "$ref" \
        --test "$test_int" \
        --corr-min "$CORR_MIN" \
        --thresholds 20,50,100 \
        | tee "$OUT_BASE/compare/compare_internal_${label}.txt"
    fi
    if [[ -f "$ref" && -f "$test_ext" ]]; then
      python3 "$COMPARE_SCRIPT" \
        --reference "$ref" \
        --test "$test_ext" \
        --corr-min "$CORR_MIN" \
        --thresholds 20,50,100 \
        | tee "$OUT_BASE/compare/compare_external_${label}.txt"
    fi
  done
else
  echo "Skipping comparison (GEDI or compare script unavailable)."
fi

echo "DONE. Outputs stored in: $OUT_BASE"
