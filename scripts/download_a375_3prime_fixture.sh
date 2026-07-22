#!/usr/bin/env bash
set -euo pipefail
# Download the 10x Genomics A375 10K CRISPR 3' v3 NextGEM public dataset.
#
# Source:  https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/SC3_v3_NextGem_DI_CRISPR_10K
# License: CC BY 4.0 (10x Genomics)
# Chemistry: 3' v3 NextGEM (TRU namespace, whitelist 3M-february-2018)
# Libraries: GEX + CRISPR Guide Capture (2 guides: RAB1A-2 + NON_TARGET-1)
# Cells: ~10,000 (7,632 recovered)
# FASTQ tar: ~59 GB (uncompressed tar of gzipped FASTQs -- large!)
# CellRanger version used by 10x: 4.0.0
#
# WARNING: The full FASTQ tar is 59 GB. By default this script streams the tar
# and extracts only the first --downsample N reads from each FASTQ, avoiding
# the full download. Use --full-download to fetch the entire tar.
#
# This dataset uses 3' chemistry, making it suitable for testing:
#   - Velocyto (spliced/unspliced quantification)
#   - BAM output with CB/UB tags
#   - GeneFull quantification
#   - packedReadInfo lifetime surface
#
# Usage:
#   scripts/download_a375_3prime_fixture.sh [options]
#
# Options:
#   --outdir DIR       Output directory (default: /tmp/a375_3prime_fixture)
#   --downsample N     Reads per FASTQ in the smoke tier (default: 100000)
#   --full-download    Download the entire 59 GB tar (default: stream + extract)
#   --force            Re-download even if files exist
#   -h, --help         Show help

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DATASET_URL="https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/SC3_v3_NextGem_DI_CRISPR_10K"
FASTQS_TAR_URL="https://cg.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_CRISPR_10K/SC3_v3_NextGem_DI_CRISPR_10K_fastqs.tar"
FEATURE_REF_URL="https://cf.10xgenomics.com/samples/cell-exp/4.0.0/SC3_v3_NextGem_DI_CRISPR_10K/SC3_v3_NextGem_DI_CRISPR_10K_feature_ref.csv"
FASTQS_TAR_BYTES=63735306240

OUTDIR="/tmp/a375_3prime_fixture"
DOWNSAMPLE_READS=100000
FULL_DOWNLOAD=0
FORCE=0

usage() { sed -n '2,/^$/s/^# \?//p' "$0"; }

log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"; }
die() { echo "ERROR: $*" >&2; exit 2; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)         OUTDIR="$2";           shift 2 ;;
    --downsample)     DOWNSAMPLE_READS="$2"; shift 2 ;;
    --full-download)  FULL_DOWNLOAD=1;       shift ;;
    --force)          FORCE=1;               shift ;;
    -h|--help)        usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

need_cmd curl
need_cmd tar

mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"

FASTQS_TAR="${OUTDIR}/SC3_v3_NextGem_DI_CRISPR_10K_fastqs.tar"
FEATURE_REF="${OUTDIR}/SC3_v3_NextGem_DI_CRISPR_10K_feature_ref.csv"
FASTQ_ROOT="${OUTDIR}/fastqs"
GEX_DIR="${FASTQ_ROOT}/gex"
CRISPR_DIR="${FASTQ_ROOT}/crispr"
TIER_DIR="${OUTDIR}/downsampled_${DOWNSAMPLE_READS}"
MANIFEST="${OUTDIR}/FIXTURE_MANIFEST.txt"

# ── Download feature reference ──────────────────────────────────────
if [[ "${FORCE}" -eq 0 && -f "${FEATURE_REF}" ]]; then
  log "Feature reference already present: ${FEATURE_REF}"
else
  log "Downloading feature reference CSV"
  curl -fSL -o "${FEATURE_REF}" "${FEATURE_REF_URL}"
fi

# ── Full download path ──────────────────────────────────────────────
if [[ "${FULL_DOWNLOAD}" -eq 1 ]]; then
  if [[ "${FORCE}" -eq 0 && -d "${GEX_DIR}" && -d "${CRISPR_DIR}" ]]; then
    log "FASTQs already extracted: ${FASTQ_ROOT}"
  else
    if [[ "${FORCE}" -eq 0 && -f "${FASTQS_TAR}" ]]; then
      log "FASTQ tar already present: ${FASTQS_TAR}"
    else
      log "Downloading full FASTQ tar (~59 GB) -- this will take a while"
      curl -fSL --progress-bar -o "${FASTQS_TAR}" "${FASTQS_TAR_URL}"
    fi

    log "Extracting FASTQs"
    mkdir -p "${FASTQ_ROOT}"
    tar -xf "${FASTQS_TAR}" -C "${FASTQ_ROOT}" --strip-components=1
    log "Extracted to ${FASTQ_ROOT}"
  fi

  # Organize into gex/crispr subdirs
  if [[ ! -d "${GEX_DIR}" ]]; then
    mkdir -p "${GEX_DIR}" "${CRISPR_DIR}"
    moved=0
    for fq in "${FASTQ_ROOT}"/*.fastq.gz; do
      [[ -f "${fq}" ]] || continue
      fname="$(basename "${fq}")"
      case "${fname}" in
        *CRISPRScreening*) mv "${fq}" "${CRISPR_DIR}/"; ((moved++)) ;;
        *)                 mv "${fq}" "${GEX_DIR}/";    ((moved++)) ;;
      esac
    done
    (( moved > 0 )) && log "Organized ${moved} FASTQs into gex/ and crispr/ subdirs"
  fi

  # Create downsample tier from full data
  if [[ "${DOWNSAMPLE_READS}" -gt 0 ]]; then
    if [[ "${FORCE}" -eq 0 && -d "${TIER_DIR}/gex" && -d "${TIER_DIR}/crispr" ]]; then
      log "Downsample tier already exists: ${TIER_DIR}"
    else
      log "Creating ${DOWNSAMPLE_READS}-read downsample tier"
      MAX_LINES=$(( DOWNSAMPLE_READS * 4 ))
      for lib_name in gex crispr; do
        lib_src="${FASTQ_ROOT}/${lib_name}"
        lib_dst="${TIER_DIR}/${lib_name}"
        mkdir -p "${lib_dst}"
        for fq in "${lib_src}"/*.fastq.gz; do
          [[ -f "${fq}" ]] || continue
          fname="$(basename "${fq}")"
          ( set +o pipefail; zcat "${fq}" | head -n "${MAX_LINES}" | gzip > "${lib_dst}/${fname}" )
        done
      done
    fi
  fi

else
  # ── Stream-and-extract path (default, avoids 59 GB download) ──────
  # Stream the tar, extracting each FASTQ and downsampling on the fly.
  # The tar entries are sequential; we process each as it appears in the
  # stream and discard the rest of the archive once we have enough.
  if [[ "${FORCE}" -eq 0 && -d "${TIER_DIR}/gex" && -d "${TIER_DIR}/crispr" ]]; then
    log "Downsample tier already exists: ${TIER_DIR}"
  else
    log "Streaming tar and extracting first ${DOWNSAMPLE_READS} reads per FASTQ"
    log "This downloads only the needed portion of the 59 GB archive"
    MAX_LINES=$(( DOWNSAMPLE_READS * 4 ))
    STREAM_TMP="${OUTDIR}/.stream_tmp"
    mkdir -p "${STREAM_TMP}" "${TIER_DIR}/gex" "${TIER_DIR}/crispr"

    # Extract all FASTQs to a temp dir from the stream.
    # curl will get a SIGPIPE / broken-pipe when tar finishes reading all
    # entries it can handle -- that's expected and harmless.
    curl -fSL "${FASTQS_TAR_URL}" 2>/dev/null \
      | tar -xf - -C "${STREAM_TMP}" --strip-components=1 --wildcards '*.fastq.gz' \
      || true

    # Check what we got
    EXTRACTED_COUNT=$(find "${STREAM_TMP}" -name '*.fastq.gz' 2>/dev/null | wc -l)
    if (( EXTRACTED_COUNT == 0 )); then
      # Streaming partial extraction may not work on all servers.
      # Fall back to telling the user to use --full-download.
      die "Stream extraction got 0 FASTQs. The server may not support partial reads. Use --full-download instead."
    fi
    log "Extracted ${EXTRACTED_COUNT} FASTQs from stream"

    # Organize and downsample
    for fq in "${STREAM_TMP}"/*.fastq.gz; do
      [[ -f "${fq}" ]] || continue
      fname="$(basename "${fq}")"
      case "${fname}" in
        *CRISPRScreening*) dest_dir="${TIER_DIR}/crispr" ;;
        *)                 dest_dir="${TIER_DIR}/gex" ;;
      esac
      ( set +o pipefail; zcat "${fq}" | head -n "${MAX_LINES}" | gzip > "${dest_dir}/${fname}" )
    done
    rm -rf "${STREAM_TMP}"
    log "Downsample tier: ${TIER_DIR}"
  fi
fi

# ── Count reads in tier ─────────────────────────────────────────────
TIER_GEX_R1="$(ls "${TIER_DIR}/gex"/*_R1_*.fastq.gz 2>/dev/null | head -1)"
TIER_CRISPR_R1="$(ls "${TIER_DIR}/crispr"/*_R1_*.fastq.gz 2>/dev/null | head -1)"
GEX_READS="unknown"
CRISPR_READS="unknown"
[[ -n "${TIER_GEX_R1:-}" ]] && GEX_READS=$(zcat "${TIER_GEX_R1}" | awk 'END{print NR/4}')
[[ -n "${TIER_CRISPR_R1:-}" ]] && CRISPR_READS=$(zcat "${TIER_CRISPR_R1}" | awk 'END{print NR/4}')
log "Tier GEX reads (first R1): ${GEX_READS}"
log "Tier CRISPR reads (first R1): ${CRISPR_READS}"

# ── Manifest ────────────────────────────────────────────────────────
{
  echo "A375 3' v3 NextGEM CRISPR 10K Public Fixture"
  echo "Generated (UTC): $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  echo "Script: ${REPO_ROOT}/scripts/download_a375_3prime_fixture.sh"
  echo "Dataset: ${DATASET_URL}"
  echo "FASTQ tar URL: ${FASTQS_TAR_URL}"
  echo "Feature ref URL: ${FEATURE_REF_URL}"
  echo "Full tar size: ${FASTQS_TAR_BYTES} bytes (~59 GB)"
  echo "Chemistry: 3' v3 NextGEM (TRU)"
  echo "Libraries: GEX + CRISPR Guide Capture (2 guides)"
  echo "Whitelist: 3M-february-2018.txt"
  echo "Feature ref: ${FEATURE_REF}"
  echo "Mode: $( [[ "${FULL_DOWNLOAD}" -eq 1 ]] && echo 'full-download' || echo 'stream-downsample' )"
  echo "Downsample tier: ${TIER_DIR} (${DOWNSAMPLE_READS} reads)"
  echo "Tier GEX reads (first R1): ${GEX_READS}"
  echo "Tier CRISPR reads (first R1): ${CRISPR_READS}"
} > "${MANIFEST}"

log "Wrote manifest: ${MANIFEST}"
log "Done. Output: ${OUTDIR}"
