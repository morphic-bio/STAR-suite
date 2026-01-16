#!/usr/bin/env bash
# Build a small, repeatable dev fixture for STAR↔GEDI SNP parity iterations:
#  1) Downsample WT FASTQ (SE 50bp) to 10M reads (deterministic seed)
#  2) Align with STAR to produce BAM (SortedByCoordinate)
#  3) Run GEDI Slam to produce *.snpdata (GEDI defaults)
#  4) Convert GEDI *.snpdata -> deduplicated 3-col BED goal set (0-based)
#
# Outputs are written to a timestamped /storage directory and include a manifest.
#
# Usage:
#   ./run_downsample10M_alignment_and_gedi_goalset.sh [--threads 8] [--seed 1] [--reads 10000000] [--overwrite 0|1]
#
set -euo pipefail

THREADS=8
SEED=1
READS=10000000
OVERWRITE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads) THREADS="${2:?}"; shift 2 ;;
    --seed) SEED="${2:?}"; shift 2 ;;
    --reads) READS="${2:?}"; shift 2 ;;
    --overwrite) OVERWRITE="${2:?}"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [--threads 8] [--seed 1] [--reads 10000000] [--overwrite 0|1]"
      exit 0
      ;;
    *)
      echo "Unknown arg: $1" >&2
      exit 2
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Inputs
FASTQ_IN="/mnt/pikachu/NW-5-21/SLAM-Seq/ARID1A-no4su_S50_R1_001.fastq.gz"
STAR_BIN="/mnt/pikachu/STAR-Flex/source/STAR"
STAR_INDEX="/storage/autoindex_110_44/bulk_index"
GEDI_BIN="/mnt/pikachu/STAR-Flex/gedi"
GEDI_GENOME="/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml"

if [[ ! -f "$FASTQ_IN" ]]; then
  echo "ERROR: FASTQ not found: $FASTQ_IN" >&2
  exit 2
fi
if [[ ! -x "$STAR_BIN" ]]; then
  echo "ERROR: STAR binary not found/executable: $STAR_BIN" >&2
  exit 2
fi
if [[ ! -d "$STAR_INDEX" ]]; then
  echo "ERROR: STAR genomeDir not found: $STAR_INDEX" >&2
  exit 2
fi
if [[ ! -x "$GEDI_BIN" ]]; then
  echo "ERROR: GEDI wrapper not found/executable: $GEDI_BIN" >&2
  exit 2
fi
if [[ ! -f "$GEDI_GENOME" ]]; then
  echo "ERROR: GEDI genome OML not found: $GEDI_GENOME" >&2
  exit 2
fi

WORK_DIR="/storage/slam_dev_arid1a_ds${READS}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$WORK_DIR"/{fastq,star,gedi,goal,logs}

FASTQ_DS="$WORK_DIR/fastq/ARID1A-no4su_S50_R1_001.ds${READS}.seed${SEED}.fastq.gz"
STAR_PREFIX="$WORK_DIR/star/ds_"
BAM="$WORK_DIR/star/ds_Aligned.sortedByCoord.out.bam"
GEDI_PREFIX="$WORK_DIR/gedi/ds"
SNPDATA="$GEDI_PREFIX.snpdata"
GOAL_BED="$WORK_DIR/goal/gedi_default.cov6.ratio0.3.p0.001.uniq.bed"
MANIFEST="$WORK_DIR/manifest.tsv"

echo "=== Downsample+Align+GEDI Goalset (dev fixture) ==="
echo "WORK_DIR: $WORK_DIR"
echo "FASTQ_IN: $FASTQ_IN"
echo "READS: $READS  SEED: $SEED  THREADS: $THREADS"
echo

echo "== [1/4] Downsample FASTQ =="
if [[ -f "$FASTQ_DS" && "$OVERWRITE" -eq 0 ]]; then
  echo "Found existing downsampled FASTQ, skipping: $FASTQ_DS"
else
  "$SCRIPT_DIR/downsample_fastq_gz.sh" \
    --in "$FASTQ_IN" \
    --out "$FASTQ_DS" \
    --reads "$READS" \
    --seed "$SEED" \
    --mode seqtk \
    > "$WORK_DIR/logs/downsample.log" 2>&1
fi
ls -lh "$FASTQ_DS"
echo

echo "== [2/4] STAR alignment (produce BAM) =="
if [[ -f "$BAM" && -f "$BAM.bai" && "$OVERWRITE" -eq 0 ]]; then
  echo "Found existing BAM+bai, skipping: $BAM"
else
  "$STAR_BIN" \
    --runThreadN "$THREADS" \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$FASTQ_DS" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$STAR_PREFIX" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS nM MD \
    --slamQuantMode 1 \
    --outFilterMultimapNmax 999 \
    > "$WORK_DIR/logs/star_align.log" 2>&1

  # Index BAM
  samtools index "$BAM"
fi
ls -lh "$BAM" "$BAM.bai"
echo

echo "== [3/4] Run GEDI Slam (GEDI defaults) =="
# We explicitly set the documented defaults for SNP calling.
if [[ -f "$SNPDATA" && "$OVERWRITE" -eq 0 ]]; then
  echo "Found existing snpdata, skipping: $SNPDATA"
else
  "$GEDI_BIN" -e Slam \
    -reads "$BAM" \
    -genomic "$GEDI_GENOME" \
    -prefix "$GEDI_PREFIX" \
    -strandness Sense \
    -err 0.001 \
    -snppval 0.001 \
    -snpConv 0.3 \
    -keep \
    -D \
    > "$WORK_DIR/logs/gedi.log" 2>&1
fi
if [[ ! -f "$SNPDATA" ]]; then
  echo "ERROR: GEDI did not produce snpdata: $SNPDATA" >&2
  echo "Check log: $WORK_DIR/logs/gedi.log" >&2
  exit 1
fi
ls -lh "$SNPDATA"
echo

echo "== [4/4] Convert GEDI snpdata -> BED goal set (dedup, 0-based) =="
python3 "$SCRIPT_DIR/convert_gedi_snpdata_to_goal_bed.py" \
  --snpdata "$SNPDATA" \
  --out "$GOAL_BED" \
  --minCov 6 \
  --minRatio 0.3 \
  --maxPval 0.001 \
  2> "$WORK_DIR/logs/goalset_convert.log"

echo "Goal BED:"
wc -l "$GOAL_BED"
head -5 "$GOAL_BED"
echo

echo "== Manifest =="
{
  echo -e "key\tvalue"
  echo -e "work_dir\t$WORK_DIR"
  echo -e "fastq_in\t$FASTQ_IN"
  echo -e "fastq_downsampled\t$FASTQ_DS"
  echo -e "reads\t$READS"
  echo -e "seed\t$SEED"
  echo -e "star_bin\t$STAR_BIN"
  echo -e "star_index\t$STAR_INDEX"
  echo -e "bam\t$BAM"
  echo -e "gedi_bin\t$GEDI_BIN"
  echo -e "gedi_genome\t$GEDI_GENOME"
  echo -e "gedi_prefix\t$GEDI_PREFIX"
  echo -e "gedi_snpdata\t$SNPDATA"
  echo -e "goal_bed\t$GOAL_BED"
  echo -e "gedi_err\t0.001"
  echo -e "gedi_snppval\t0.001"
  echo -e "gedi_snpConv\t0.3"
  echo -e "gedi_strandness\tSense"
} | tee "$MANIFEST"

echo
echo "DONE."
echo "BAM: $BAM"
echo "Goal BED: $GOAL_BED"
echo "Logs: $WORK_DIR/logs/"

