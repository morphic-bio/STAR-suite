#!/usr/bin/env bash
#
# Run VB-weighted SLAM vs GEDI using a GEDI-derived common SNP mask.
#
# This script:
#  1) Runs STAR on 6h with TranscriptVB + SNP mask + auto-trim (builds dumps + vbGene weights)
#  2) Runs STAR on 0h with same trims + SNP mask (builds dumps + vbGene weights)
#  3) Re-quantifies both dumps via slam_requant using vbGene weights
#  4) Runs GEDI on both BAMs (same trims + err)
#  5) Compares GEDI vs VB-weighted slam_requant and computes delta correlation (6h-0h)
#
# Notes:
# - GEDI does not accept an external SNP mask. "Common mask" here means a GEDI-
#   derived mask applied to STAR; GEDI uses its own internal SNP mask logic.
# - Default dump cap is 1,000,000 reads. Increase via --dump-max for full parity.
#
# Usage:
#   bash tests/run_slam_end_to_end_vb_commonmask.sh [options]
#
# Options:
#   --work-dir DIR     Working directory (default: /storage/slam_e2e_vb_commonmask_YYYYMMDD_HHMMSS)
#   --threads N        Threads for STAR (default: 8)
#   --dump-max N       Max reads to dump (default: 1000000)
#   --common-mask PATH Path to GEDI-derived mask BED.gz
#   --overwrite        Rerun even if outputs exist
#   --fail-fast        Exit immediately on first failure
#   --help             Show this help message

set -o pipefail

WORK_DIR="${WORK_DIR:-/storage/slam_e2e_vb_commonmask_$(date +%Y%m%d_%H%M%S)}"
THREADS="${THREADS:-8}"
DUMP_MAX="${DUMP_MAX:-1000000}"
OVERWRITE="${OVERWRITE:-0}"
FAIL_FAST="${FAIL_FAST:-0}"

COMMON_MASK="${COMMON_MASK:-}"

export WORK_DIR

READY=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --dump-max) DUMP_MAX="$2"; shift 2 ;;
    --common-mask) COMMON_MASK="$2"; shift 2 ;;
    --overwrite) OVERWRITE=1; shift ;;
    --fail-fast) FAIL_FAST=1; shift ;;
    --help)
      grep "^# " "$0" | tail -n +4
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

STAR_BIN="${STAR_BIN:-$PROJECT_ROOT/core/legacy/source/STAR}"
REQ_BIN="${REQ_BIN:-$PROJECT_ROOT/tools/slam_requant/slam_requant}"
GEDI_BIN="${GEDI_BIN:-$PROJECT_ROOT/gedi}"
STAR_INDEX="${STAR_INDEX:-/storage/autoindex_110_44/bulk_index}"
GEDI_GENOME="${GEDI_GENOME:-/home/lhhung/.gedi/genomic/homo_sapiens_110_44.oml}"

WT_FASTQ="${WT_FASTQ:-/storage/SLAM-Seq-prod-compare-20260109/input/WDHD1-0h-3_S201_R1_001.fastq.gz}"
H6_FASTQ="${H6_FASTQ:-/storage/SLAM-Seq-prod-compare-20260109/input/ARID1A-6h-1_S43_R1_001.fastq.gz}"

QUANT_VB_LIBTYPE="${QUANT_VB_LIBTYPE:-U}"
SLAM_STRANDNESS="${SLAM_STRANDNESS:-Sense}"
REQUANT_STRANDNESS="${REQUANT_STRANDNESS:-sense}"

mkdir -p "$WORK_DIR"/{mask,star,requant,gedi,qc,report}
LOG_FILE="$WORK_DIR/report/run_vb_commonmask.log"
exec > >(tee -a "$LOG_FILE") 2>&1

FAILURES=()

log() { echo "[$(date +%H:%M:%S)] $*"; }

record_fail() {
  FAILURES+=("$1")
  if [[ "$FAIL_FAST" -eq 1 ]]; then
    log "FAIL_FAST enabled; exiting."
    exit 1
  fi
}

run_cmd_log() {
  local label="$1"
  local log_path="$2"
  shift 2
  log "$label"
  "$@" > "$log_path" 2>&1
  local rc=$?
  if [[ $rc -ne 0 ]]; then
    log "ERROR: $label failed (exit $rc). See $log_path"
    record_fail "$label"
  fi
  return $rc
}

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "$path" ]]; then
    log "ERROR: missing $label: $path"
    record_fail "missing:$label"
    READY=0
    return 1
  fi
  return 0
}

ensure_bam_index() {
  local bam="$1"
  if [[ -f "${bam}.bai" ]]; then
    return 0
  fi
  if ! command -v samtools >/dev/null 2>&1; then
    log "WARNING: samtools not found; cannot index $bam for GEDI"
    record_fail "samtools-missing"
    return 1
  fi
  log "Indexing BAM for GEDI: $bam"
  samtools index "$bam"
}

if [[ -z "$COMMON_MASK" ]]; then
  if [[ -f "/storage/gedi_common_mask_20260114_202352/gedi_mask.bed.gz" ]]; then
    COMMON_MASK="/storage/gedi_common_mask_20260114_202352/gedi_mask.bed.gz"
  fi
fi

log "========================================================================"
log "SLAM VB-Weighted End-to-End (Common GEDI Mask)"
log "========================================================================"
log "Work dir: $WORK_DIR"
log "STAR: $STAR_BIN"
log "slam_requant: $REQ_BIN"
log "GEDI: $GEDI_BIN"
log "STAR index: $STAR_INDEX"
log "0h FASTQ: $WT_FASTQ"
log "6h FASTQ: $H6_FASTQ"
log "Common mask: ${COMMON_MASK:-UNSET}"
log "Threads: $THREADS"
log "Dump max reads: $DUMP_MAX"
log "SLAM strandness (STAR): $SLAM_STRANDNESS"
log "Strandness (slam_requant): $REQUANT_STRANDNESS"
log "Fail fast: $FAIL_FAST"
log "========================================================================"

require_file "$STAR_BIN" "STAR binary"
require_file "$REQ_BIN" "slam_requant binary"
require_file "$STAR_INDEX/Genome" "STAR index"
require_file "$WT_FASTQ" "0h FASTQ"
require_file "$H6_FASTQ" "6h FASTQ"

if [[ -n "$COMMON_MASK" ]]; then
  require_file "$COMMON_MASK" "common SNP mask"
else
  log "ERROR: common SNP mask not set. Use --common-mask or export COMMON_MASK."
  record_fail "missing:common-mask"
  READY=0
fi

if [[ "$READY" -eq 0 ]]; then
  log "Missing required inputs; stopping before run."
  log "Check log: $LOG_FILE"
  exit 0
fi

RUN_GEDI=1
if [[ ! -x "$GEDI_BIN" ]]; then
  log "WARNING: GEDI binary not found; skipping GEDI runs."
  RUN_GEDI=0
elif [[ ! -f "$GEDI_GENOME" ]]; then
  log "WARNING: GEDI genome not found; skipping GEDI runs."
  RUN_GEDI=0
fi

STAR_6H_OUT="$WORK_DIR/star/6h_SlamQuant.out"
STAR_0H_OUT="$WORK_DIR/star/0h_SlamQuant.out"
QC_PREFIX="$WORK_DIR/qc/trim_6h"
TRIM_JSON="$QC_PREFIX.slam_qc.json"
TRIM_HTML="$QC_PREFIX.slam_qc.html"

log "[1/7] STAR 6h (TranscriptVB + mask + auto-trim + dump)"
if [[ -f "$STAR_6H_OUT" && -f "$WORK_DIR/star/6h_Aligned.sortedByCoord.out.bam" \
      && -f "$TRIM_JSON" && "$OVERWRITE" -eq 0 ]]; then
  log "✓ STAR 6h outputs exist, skipping (use --overwrite to rerun)"
else
  run_cmd_log "STAR 6h" "$WORK_DIR/report/star_6h_vb.log" \
    "$STAR_BIN" \
      --runThreadN "$THREADS" \
      --genomeDir "$STAR_INDEX" \
      --readFilesIn "$H6_FASTQ" \
      --readFilesCommand zcat \
      --outFileNamePrefix "$WORK_DIR/star/6h_" \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS nM MD \
      --slamQuantMode 1 \
      --slamStrandness "$SLAM_STRANDNESS" \
      --slamSnpMaskIn "$COMMON_MASK" \
      --autoTrim variance \
      --trimSource "$H6_FASTQ" \
      --autoTrimDetectionReads 100000 \
      --slamQcReport "$QC_PREFIX" \
      --quantMode TranscriptVB \
      --quantVBLibType "$QUANT_VB_LIBTYPE" \
      --slamDumpBinary "$WORK_DIR/requant/6h.dump.bin" \
      --slamDumpWeights "$WORK_DIR/requant/6h.weights.bin" \
      --slamDumpWeightsMode vbGene \
      --slamDumpMaxReads "$DUMP_MAX"
fi

if [[ ! -f "$STAR_6H_OUT" ]]; then
  log "ERROR: STAR 6h output missing: $STAR_6H_OUT"
  record_fail "star-6h-missing"
fi

if [[ ! -f "$TRIM_JSON" ]]; then
  log "ERROR: QC JSON missing: $TRIM_JSON"
  record_fail "trim-json-missing"
fi

TRIM5=0
TRIM3=0
SNP_ERR="0.001"
if [[ -f "$TRIM_JSON" ]]; then
  TRIM5=$(python3 - "$TRIM_JSON" << 'PY'
import json, sys
data = json.load(open(sys.argv[1]))
print(int(data.get("trim5p", 0)))
PY
)
  TRIM3=$(python3 - "$TRIM_JSON" << 'PY'
import json, sys
data = json.load(open(sys.argv[1]))
print(int(data.get("trim3p", 0)))
PY
)
  SNP_ERR=$(python3 - "$TRIM_JSON" << 'PY'
import json, sys
data = json.load(open(sys.argv[1]))
print(f"{data.get('snp_err_used', 0.001):.6f}")
PY
)
fi

log "Trim values: trim5p=$TRIM5 trim3p=$TRIM3"
log "SNP error rate: snp_err_used=$SNP_ERR"

log "[2/7] STAR 0h (TranscriptVB + mask + fixed trims + dump)"
if [[ -f "$STAR_0H_OUT" && -f "$WORK_DIR/star/0h_Aligned.sortedByCoord.out.bam" \
      && "$OVERWRITE" -eq 0 ]]; then
  log "✓ STAR 0h outputs exist, skipping (use --overwrite to rerun)"
else
  run_cmd_log "STAR 0h" "$WORK_DIR/report/star_0h_vb.log" \
    "$STAR_BIN" \
      --runThreadN "$THREADS" \
      --genomeDir "$STAR_INDEX" \
      --readFilesIn "$WT_FASTQ" \
      --readFilesCommand zcat \
      --outFileNamePrefix "$WORK_DIR/star/0h_" \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS nM MD \
      --slamQuantMode 1 \
      --slamStrandness "$SLAM_STRANDNESS" \
      --slamSnpMaskIn "$COMMON_MASK" \
      --autoTrim - \
      --slamCompatTrim5p "$TRIM5" \
      --slamCompatTrim3p "$TRIM3" \
      --quantMode TranscriptVB \
      --quantVBLibType "$QUANT_VB_LIBTYPE" \
      --slamDumpBinary "$WORK_DIR/requant/0h.dump.bin" \
      --slamDumpWeights "$WORK_DIR/requant/0h.weights.bin" \
      --slamDumpWeightsMode vbGene \
      --slamDumpMaxReads "$DUMP_MAX"
fi

if [[ ! -f "$STAR_0H_OUT" ]]; then
  log "ERROR: STAR 0h output missing: $STAR_0H_OUT"
  record_fail "star-0h-missing"
fi

log "[3/7] slam_requant with vbGene weights"
if [[ -f "$WORK_DIR/requant/6h_SlamQuant.out" && "$OVERWRITE" -eq 0 ]]; then
  log "✓ requant 6h outputs exist, skipping (use --overwrite to rerun)"
else
  run_cmd_log "slam_requant 6h" "$WORK_DIR/report/requant_6h.log" \
    "$REQ_BIN" \
      --dump "$WORK_DIR/requant/6h.dump.bin" \
      --out "$WORK_DIR/requant/6h_" \
      --slamSnpMaskIn "$COMMON_MASK" \
      --trim5p "$TRIM5" \
      --trim3p "$TRIM3" \
      --strandness "$REQUANT_STRANDNESS" \
      --slamWeightFile "$WORK_DIR/requant/6h.weights.bin" \
      --slamWeightMatch key
fi

if [[ -f "$WORK_DIR/requant/0h_SlamQuant.out" && "$OVERWRITE" -eq 0 ]]; then
  log "✓ requant 0h outputs exist, skipping (use --overwrite to rerun)"
else
  run_cmd_log "slam_requant 0h" "$WORK_DIR/report/requant_0h.log" \
    "$REQ_BIN" \
      --dump "$WORK_DIR/requant/0h.dump.bin" \
      --out "$WORK_DIR/requant/0h_" \
      --slamSnpMaskIn "$COMMON_MASK" \
      --trim5p "$TRIM5" \
      --trim3p "$TRIM3" \
      --strandness "$REQUANT_STRANDNESS" \
      --slamWeightFile "$WORK_DIR/requant/0h.weights.bin" \
      --slamWeightMatch key
fi

if [[ -f "$WORK_DIR/requant/6h_normal_SlamQuant.out" && "$OVERWRITE" -eq 0 ]]; then
  log "✓ requant normal 6h outputs exist, skipping (use --overwrite to rerun)"
else
  run_cmd_log "slam_requant 6h (normal weights)" "$WORK_DIR/report/requant_6h_normal.log" \
    "$REQ_BIN" \
      --dump "$WORK_DIR/requant/6h.dump.bin" \
      --out "$WORK_DIR/requant/6h_normal_" \
      --slamSnpMaskIn "$COMMON_MASK" \
      --trim5p "$TRIM5" \
      --trim3p "$TRIM3" \
      --strandness "$REQUANT_STRANDNESS" \
      --slamWeightMode dump
fi

if [[ -f "$WORK_DIR/requant/0h_normal_SlamQuant.out" && "$OVERWRITE" -eq 0 ]]; then
  log "✓ requant normal 0h outputs exist, skipping (use --overwrite to rerun)"
else
  run_cmd_log "slam_requant 0h (normal weights)" "$WORK_DIR/report/requant_0h_normal.log" \
    "$REQ_BIN" \
      --dump "$WORK_DIR/requant/0h.dump.bin" \
      --out "$WORK_DIR/requant/0h_normal_" \
      --slamSnpMaskIn "$COMMON_MASK" \
      --trim5p "$TRIM5" \
      --trim3p "$TRIM3" \
      --strandness "$REQUANT_STRANDNESS" \
      --slamWeightMode dump
fi

log "[4/7] GEDI runs (optional)"
if [[ "$RUN_GEDI" -eq 1 ]]; then
  ensure_bam_index "$WORK_DIR/star/0h_Aligned.sortedByCoord.out.bam"
  ensure_bam_index "$WORK_DIR/star/6h_Aligned.sortedByCoord.out.bam"

  if [[ -f "$WORK_DIR/gedi/0h.tsv.gz" && "$OVERWRITE" -eq 0 ]]; then
    log "✓ GEDI 0h output exists, skipping"
  else
    run_cmd_log "GEDI 0h" "$WORK_DIR/report/gedi_0h.log" \
      "$GEDI_BIN" -e Slam \
        -reads "$WORK_DIR/star/0h_Aligned.sortedByCoord.out.bam" \
        -genomic "$GEDI_GENOME" \
        -prefix "$WORK_DIR/gedi/0h" \
        -trim5p "$TRIM5" \
        -trim3p "$TRIM3" \
        -strandness Sense \
        -err "$SNP_ERR" \
        -D
  fi

  if [[ -f "$WORK_DIR/gedi/6h.tsv.gz" && "$OVERWRITE" -eq 0 ]]; then
    log "✓ GEDI 6h output exists, skipping"
  else
    run_cmd_log "GEDI 6h" "$WORK_DIR/report/gedi_6h.log" \
      "$GEDI_BIN" -e Slam \
        -reads "$WORK_DIR/star/6h_Aligned.sortedByCoord.out.bam" \
        -genomic "$GEDI_GENOME" \
        -prefix "$WORK_DIR/gedi/6h" \
        -trim5p "$TRIM5" \
        -trim3p "$TRIM3" \
        -strandness Sense \
        -err "$SNP_ERR" \
        -D
  fi
else
  log "Skipping GEDI (binary or genome missing)."
fi

log "[5/7] STAR vs GEDI comparison (VB-weighted slam_requant)"
COMPARE_SCRIPT="$SCRIPT_DIR/slam/compare_fixture.py"
if [[ "$RUN_GEDI" -eq 1 && -f "$COMPARE_SCRIPT" ]]; then
  run_cmd_log "Compare 0h" "$WORK_DIR/report/compare_vb_0h.txt" \
    python3 "$COMPARE_SCRIPT" \
      --reference "$WORK_DIR/gedi/0h.tsv.gz" \
      --test "$WORK_DIR/requant/0h_SlamQuant.out" \
      --thresholds 20,50,100

  run_cmd_log "Compare 6h" "$WORK_DIR/report/compare_vb_6h.txt" \
    python3 "$COMPARE_SCRIPT" \
      --reference "$WORK_DIR/gedi/6h.tsv.gz" \
      --test "$WORK_DIR/requant/6h_SlamQuant.out" \
      --thresholds 20,50,100

  run_cmd_log "Compare 0h (normal)" "$WORK_DIR/report/compare_normal_0h.txt" \
    python3 "$COMPARE_SCRIPT" \
      --reference "$WORK_DIR/gedi/0h.tsv.gz" \
      --test "$WORK_DIR/requant/0h_normal_SlamQuant.out" \
      --thresholds 20,50,100

  run_cmd_log "Compare 6h (normal)" "$WORK_DIR/report/compare_normal_6h.txt" \
    python3 "$COMPARE_SCRIPT" \
      --reference "$WORK_DIR/gedi/6h.tsv.gz" \
      --test "$WORK_DIR/requant/6h_normal_SlamQuant.out" \
      --thresholds 20,50,100
else
  log "Skipping comparison (GEDI missing or compare script not found)."
fi

log "[6/7] Delta correlation (6h - 0h)"
if [[ "$RUN_GEDI" -eq 1 && -f "$WORK_DIR/gedi/0h.tsv.gz" && -f "$WORK_DIR/gedi/6h.tsv.gz" \
      && -f "$WORK_DIR/requant/0h_SlamQuant.out" && -f "$WORK_DIR/requant/6h_SlamQuant.out" ]]; then
  python3 - << 'PY' | tee "$WORK_DIR/report/delta_corr_vb.txt"
import csv
import gzip
import math
import os

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def find_col(header, names):
    lower = {name.lower(): i for i, name in enumerate(header)}
    for name in names:
        idx = lower.get(name.lower())
        if idx is not None:
            return idx
    return None

def find_suffix(header, suffixes):
    for suffix in suffixes:
        for i, name in enumerate(header):
            if name.endswith(" " + suffix) or name.endswith(suffix):
                return i
    return None

def parse_ref(path):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        gene_idx = find_col(header, ["Gene", "GeneID"])
        read_idx = find_suffix(header, ["Readcount"])
        map_idx = find_suffix(header, ["MAP"])
        data = {}
        for row in reader:
            if not row or len(row) <= map_idx:
                continue
            gene = row[gene_idx]
            if not gene:
                continue
            data[gene] = {
                "readcount": float(row[read_idx]),
                "ntr": float(row[map_idx]),
            }
        return data

def parse_test(path):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        gene_idx = find_col(header, ["Gene", "GeneID"])
        read_idx = find_col(header, ["ReadCount", "Readcount"])
        ntr_idx = find_col(header, ["NTR", "NTR_MAP", "MAP"])
        data = {}
        for row in reader:
            if not row or len(row) <= ntr_idx:
                continue
            gene = row[gene_idx]
            if not gene:
                continue
            data[gene] = {
                "readcount": float(row[read_idx]),
                "ntr": float(row[ntr_idx]),
            }
        return data

def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = sum((x - mean_x) ** 2 for x in xs)
    den_y = sum((y - mean_y) ** 2 for y in ys)
    den = math.sqrt(den_x * den_y)
    return float("nan") if den == 0.0 else num / den

def spearman(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    rank_x = sorted(range(n), key=lambda i: xs[i])
    rank_y = sorted(range(n), key=lambda i: ys[i])
    rank_x_dict = {rank_x[i]: i + 1 for i in range(n)}
    rank_y_dict = {rank_y[i]: i + 1 for i in range(n)}
    d_sq = sum((rank_x_dict[i] - rank_y_dict[i]) ** 2 for i in range(n))
    return 1.0 - (6.0 * d_sq) / (n * (n * n - 1))

work_dir = os.environ["WORK_DIR"]
gedi_0h = parse_ref(f"{work_dir}/gedi/0h.tsv.gz")
gedi_6h = parse_ref(f"{work_dir}/gedi/6h.tsv.gz")
star_0h = parse_test(f"{work_dir}/requant/0h_SlamQuant.out")
star_6h = parse_test(f"{work_dir}/requant/6h_SlamQuant.out")

shared = sorted(set(gedi_0h) & set(gedi_6h) & set(star_0h) & set(star_6h))
thresholds = [0, 20, 50, 100]

print("=== Delta Correlation (GEDI vs STAR VB-weighted, 6h - 0h) ===")
print("Filter: min(gedi_0h, gedi_6h, star_0h, star_6h) readcount >= threshold")
print(f"{'Threshold':<10} {'N Genes':<10} {'Delta Pearson':<14} {'Delta Spearman':<14}")
print("-" * 60)

for thr in thresholds:
    xs, ys = [], []
    for gene in shared:
        if min(gedi_0h[gene]["readcount"], gedi_6h[gene]["readcount"],
               star_0h[gene]["readcount"], star_6h[gene]["readcount"]) < thr:
            continue
        delta_ref = gedi_6h[gene]["ntr"] - gedi_0h[gene]["ntr"]
        delta_test = star_6h[gene]["ntr"] - star_0h[gene]["ntr"]
        xs.append(delta_ref)
        ys.append(delta_test)
    pear = pearson(xs, ys)
    spear = spearman(xs, ys)
    print(f">={thr:<8} {len(xs):<10} {pear:>13.6f} {spear:>13.6f}")
PY
else
  log "Skipping delta correlation (inputs missing)."
fi

if [[ "$RUN_GEDI" -eq 1 && -f "$WORK_DIR/gedi/0h.tsv.gz" && -f "$WORK_DIR/gedi/6h.tsv.gz" \
      && -f "$WORK_DIR/requant/0h_normal_SlamQuant.out" && -f "$WORK_DIR/requant/6h_normal_SlamQuant.out" ]]; then
  python3 - << 'PY' | tee "$WORK_DIR/report/delta_corr_normal.txt"
import csv
import gzip
import math
import os

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def find_col(header, names):
    lower = {name.lower(): i for i, name in enumerate(header)}
    for name in names:
        idx = lower.get(name.lower())
        if idx is not None:
            return idx
    return None

def find_suffix(header, suffixes):
    for suffix in suffixes:
        for i, name in enumerate(header):
            if name.endswith(" " + suffix) or name.endswith(suffix):
                return i
    return None

def parse_ref(path):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        gene_idx = find_col(header, ["Gene", "GeneID"])
        read_idx = find_suffix(header, ["Readcount"])
        map_idx = find_suffix(header, ["MAP"])
        data = {}
        for row in reader:
            if not row or len(row) <= map_idx:
                continue
            gene = row[gene_idx]
            if not gene:
                continue
            data[gene] = {
                "readcount": float(row[read_idx]),
                "ntr": float(row[map_idx]),
            }
        return data

def parse_test(path):
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        gene_idx = find_col(header, ["Gene", "GeneID"])
        read_idx = find_col(header, ["ReadCount", "Readcount"])
        ntr_idx = find_col(header, ["NTR", "NTR_MAP", "MAP"])
        data = {}
        for row in reader:
            if not row or len(row) <= ntr_idx:
                continue
            gene = row[gene_idx]
            if not gene:
                continue
            data[gene] = {
                "readcount": float(row[read_idx]),
                "ntr": float(row[ntr_idx]),
            }
        return data

def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den_x = sum((x - mean_x) ** 2 for x in xs)
    den_y = sum((y - mean_y) ** 2 for y in ys)
    den = math.sqrt(den_x * den_y)
    return float("nan") if den == 0.0 else num / den

def spearman(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    rank_x = sorted(range(n), key=lambda i: xs[i])
    rank_y = sorted(range(n), key=lambda i: ys[i])
    rank_x_dict = {rank_x[i]: i + 1 for i in range(n)}
    rank_y_dict = {rank_y[i]: i + 1 for i in range(n)}
    d_sq = sum((rank_x_dict[i] - rank_y_dict[i]) ** 2 for i in range(n))
    return 1.0 - (6.0 * d_sq) / (n * (n * n - 1))

work_dir = os.environ["WORK_DIR"]
gedi_0h = parse_ref(f"{work_dir}/gedi/0h.tsv.gz")
gedi_6h = parse_ref(f"{work_dir}/gedi/6h.tsv.gz")
star_0h = parse_test(f"{work_dir}/requant/0h_normal_SlamQuant.out")
star_6h = parse_test(f"{work_dir}/requant/6h_normal_SlamQuant.out")

shared = sorted(set(gedi_0h) & set(gedi_6h) & set(star_0h) & set(star_6h))
thresholds = [0, 20, 50, 100]

print("=== Delta Correlation (GEDI vs STAR normal weights, 6h - 0h) ===")
print("Filter: min(gedi_0h, gedi_6h, star_0h, star_6h) readcount >= threshold")
print(f"{'Threshold':<10} {'N Genes':<10} {'Delta Pearson':<14} {'Delta Spearman':<14}")
print("-" * 60)

for thr in thresholds:
    xs, ys = [], []
    for gene in shared:
        if min(gedi_0h[gene]["readcount"], gedi_6h[gene]["readcount"],
               star_0h[gene]["readcount"], star_6h[gene]["readcount"]) < thr:
            continue
        delta_ref = gedi_6h[gene]["ntr"] - gedi_0h[gene]["ntr"]
        delta_test = star_6h[gene]["ntr"] - star_0h[gene]["ntr"]
        xs.append(delta_ref)
        ys.append(delta_test)
    pear = pearson(xs, ys)
    spear = spearman(xs, ys)
    print(f">={thr:<8} {len(xs):<10} {pear:>13.6f} {spear:>13.6f}")
PY
else
  log "Skipping delta correlation (normal weights inputs missing)."
fi

log "[7/7] Done."

if [[ ${#FAILURES[@]} -gt 0 ]]; then
  log "Completed with warnings/errors:"
  for f in "${FAILURES[@]}"; do
    log "  - $f"
  done
  log "Check log: $LOG_FILE"
  if [[ "$FAIL_FAST" -eq 1 ]]; then
    exit 1
  fi
else
  log "All steps completed successfully."
fi
