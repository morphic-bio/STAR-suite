#!/usr/bin/env bash
# Canonical GEX-only STAR invocation for Velocyto exactness debugging.
# This intentionally removes pfMultiConfig / feature calling and runs only the
# GEX library so Velocyto correctness can be debugged on the smallest useful surface.
#
# Usage:
#   run_star_velocyto_gexonly_canonical.sh --gex-dir /path/to/GEX --threads N --out-prefix /path/to/run/ [--prepare-mex]
#
# Optional env:
#   STAR_BIN
#   UCSF_VELOCYTO_GEX_GENOME_DIR
#   UCSF_VELOCYTO_GEX_CB_WHITELIST
#   UCSF_VELOCYTO_GEX_STRAND               (default: Forward)
#   UCSF_VELOCYTO_GEX_MULTI_MAPPERS        (default: Unique)
#   UCSF_VELOCYTO_GEX_READFILES_COMMAND    (default: zcat)
#   STAR_VELOCYTO_DETERMINISTIC_REPLAY
#   STAR_VELOCYTO_INTEGRATED_HASH
#   UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1      (only if intentionally reusing an outdir)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
GENOME_DIR="${UCSF_VELOCYTO_GEX_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
WHITELIST="${UCSF_VELOCYTO_GEX_CB_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt}"
SOLO_STRAND="${UCSF_VELOCYTO_GEX_STRAND:-Forward}"
SOLO_MULTI_MAPPERS="${UCSF_VELOCYTO_GEX_MULTI_MAPPERS:-Unique}"
READFILES_COMMAND="${UCSF_VELOCYTO_GEX_READFILES_COMMAND:-zcat}"

GEX_DIR=""
THREADS=""
OUT_PREFIX=""
PREPARE_MEX=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gex-dir)
      GEX_DIR="${2:?}"
      shift 2
      ;;
    --threads)
      THREADS="${2:?}"
      shift 2
      ;;
    --out-prefix)
      OUT_PREFIX="${2:?}"
      shift 2
      ;;
    --prepare-mex)
      PREPARE_MEX=1
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

die() {
  echo "ERROR: $*" >&2
  exit 1
}

join_by_comma() {
  local IFS=,
  echo "$*"
}

[[ -n "${GEX_DIR}" ]] || die "--gex-dir is required"
[[ -n "${THREADS}" ]] || die "--threads is required"
[[ -n "${OUT_PREFIX}" ]] || die "--out-prefix is required"
[[ -x "${STAR_BIN}" ]] || die "STAR not executable: ${STAR_BIN}"
[[ -d "${GEX_DIR}" ]] || die "GEX dir missing: ${GEX_DIR}"
[[ -f "${WHITELIST}" ]] || die "Whitelist missing: ${WHITELIST}"
[[ -d "${GENOME_DIR}" ]] || die "Genome dir missing: ${GENOME_DIR}"

# -L follows the corrected EBs2_2 symlinked FASTQ layout.
mapfile -t R2_FILES < <(find -L "${GEX_DIR}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' | sort)
mapfile -t R1_FILES < <(find -L "${GEX_DIR}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort)
[[ "${#R2_FILES[@]}" -gt 0 ]] || die "No R2 FASTQs under ${GEX_DIR}"
[[ "${#R2_FILES[@]}" -eq "${#R1_FILES[@]}" ]] || die "R1/R2 FASTQ count mismatch under ${GEX_DIR}"

R2_LIST="$(join_by_comma "${R2_FILES[@]}")"
R1_LIST="$(join_by_comma "${R1_FILES[@]}")"
OUT_PREFIX="${OUT_PREFIX%/}"

if [[ "${UCSF_VELOCYTO_REUSE_STAR_OUTDIR:-0}" != "1" ]]; then
  if [[ -e "${OUT_PREFIX}/Solo.out" || -e "${OUT_PREFIX}/Log.out" || -e "${OUT_PREFIX}/Aligned.out.bam" || -e "${OUT_PREFIX}/Aligned.out.sam" ]]; then
    die "Refusing non-fresh --out-prefix ${OUT_PREFIX} (Solo.out, Log.out, or aligned output present). Use a new directory or UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1"
  fi
  if [[ -d "${OUT_PREFIX}/outs/raw_velocyto_feature_bc_matrix" || -d "${OUT_PREFIX}/outs/filtered_velocyto_feature_bc_matrix" ]]; then
    die "Refusing non-fresh --out-prefix ${OUT_PREFIX} (outs/*velocyto* present). Use a new directory or UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1"
  fi
fi

mkdir -p "${OUT_PREFIX}"

echo "=== run_star_velocyto_gexonly_canonical ==="
echo "gex_dir=${GEX_DIR}"
echo "threads=${THREADS}"
echo "out=${OUT_PREFIX}"
echo "strand=${SOLO_STRAND} multimappers=${SOLO_MULTI_MAPPERS} readFilesCommand=${READFILES_COMMAND}"
if [[ -n "${STAR_VELOCYTO_DETERMINISTIC_REPLAY:-}" ]]; then
  echo "STAR_VELOCYTO_DETERMINISTIC_REPLAY=${STAR_VELOCYTO_DETERMINISTIC_REPLAY}"
fi
if [[ -n "${STAR_VELOCYTO_INTEGRATED_HASH:-}" ]]; then
  echo "STAR_VELOCYTO_INTEGRATED_HASH=${STAR_VELOCYTO_INTEGRATED_HASH}"
fi

"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${R2_LIST}" "${R1_LIST}" \
  --readFilesCommand "${READFILES_COMMAND}" \
  --outFileNamePrefix "${OUT_PREFIX}/" \
  --outSAMtype None \
  --clipAdapterType CellRanger4 \
  --clip3pPolyG yes \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${WHITELIST}" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers "${SOLO_MULTI_MAPPERS}" \
  --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no \
  --soloStrand "${SOLO_STRAND}" \
  --soloFeatures Gene GeneFull Velocyto

if [[ "${PREPARE_MEX}" -eq 1 ]]; then
  [[ -f "${OUT_PREFIX}/outs/velocyto_feature_bc_matrix_manifest.json" ]] || die "Missing native Velocyto MEX manifest in ${OUT_PREFIX}/outs"
  [[ -f "${OUT_PREFIX}/outs/raw_velocyto_feature_bc_matrix/matrix.mtx.gz" ]] || die "Missing native raw Velocyto MEX in ${OUT_PREFIX}/outs"
  [[ -f "${OUT_PREFIX}/outs/filtered_velocyto_feature_bc_matrix/matrix.mtx.gz" ]] || die "Missing native filtered Velocyto MEX in ${OUT_PREFIX}/outs"
fi

echo "OK: ${OUT_PREFIX}"
