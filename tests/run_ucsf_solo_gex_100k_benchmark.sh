#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

SAMPLE_ROOT="${UCSF_SAMPLE_ROOT:-/mnt/pikachu/ucsf-perturb-seq-corrected/iPSC2_1}"
GEX_DIR="${UCSF_GEX_DIR:-${SAMPLE_ROOT}/GEX}"
TOTAL_READS="${UCSF_GEX_TOTAL_READS:-100000}"
THREADS="${UCSF_GEX_THREADS:-16}"
GENOME_DIR="${UCSF_GEX_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
WHITELIST="${UCSF_GEX_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt}"
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
ARTIFACT_ROOT="${UCSF_GEX_ARTIFACT_ROOT:-/storage/100K/ucsf_solo_optimization_20260324/iPSC2_1_GEX_100k_total}"
FASTQ_DIR="${UCSF_GEX_FASTQ_DIR:-${ARTIFACT_ROOT}/fastq_downsampled}"
OUTPREFIX="${UCSF_GEX_OUTPREFIX:-${ARTIFACT_ROOT}/run/}"
TMPDIR="${UCSF_GEX_TMPDIR:-${ARTIFACT_ROOT}/tmp}"
TIME_LOG="${UCSF_GEX_TIME_LOG:-${ARTIFACT_ROOT}/time.txt}"
EXTRA_STAR_ARGS_STR="${UCSF_GEX_EXTRA_STAR_ARGS:-}"
SOLO_UMI_FILTERING="${UCSF_GEX_SOLO_UMI_FILTERING:-MultiGeneUMI_CR}"
SOLO_UMI_DEDUP="${UCSF_GEX_SOLO_UMI_DEDUP:-1MM_CR}"
SOLO_MULTI_MAPPERS="${UCSF_GEX_SOLO_MULTI_MAPPERS:-Rescue}"
SOLO_CELL_FILTER="${UCSF_GEX_SOLO_CELL_FILTER:-EmptyDrops_CR}"
SOLO_INLINE_HASH_MODE="${UCSF_GEX_SOLO_INLINE_HASH_MODE:-}"
EXTRA_STAR_ARGS=()

if [[ -n "${EXTRA_STAR_ARGS_STR}" ]]; then
  # shellcheck disable=SC2206
  EXTRA_STAR_ARGS=( ${EXTRA_STAR_ARGS_STR} )
fi

if [[ "${OUTPREFIX}" != */ ]]; then
  OUTPREFIX="${OUTPREFIX}/"
fi

if [[ ! -d "${GEX_DIR}" ]]; then
  echo "Missing GEX directory: ${GEX_DIR}" >&2
  exit 1
fi

if [[ ! -d "${GENOME_DIR}" ]]; then
  echo "Missing genome index: ${GENOME_DIR}" >&2
  exit 1
fi

if [[ ! -f "${WHITELIST}" ]]; then
  echo "Missing whitelist: ${WHITELIST}" >&2
  exit 1
fi

if [[ ! -x "${STAR_BIN}" ]]; then
  echo "STAR binary not found at ${STAR_BIN}. Build STAR first." >&2
  exit 1
fi

mkdir -p "${ARTIFACT_ROOT}" "${FASTQ_DIR}" "${OUTPREFIX}" "$(dirname "${TMPDIR}")"

lane_count=$(find -L "${GEX_DIR}" -maxdepth 1 -type f -name "*_R1_*.fastq.gz" | wc -l | tr -d ' ')
if [[ "${lane_count}" -eq 0 ]]; then
  echo "No R1 FASTQs found under ${GEX_DIR}" >&2
  exit 1
fi

reads_per_fastq=$(( (TOTAL_READS + lane_count - 1) / lane_count ))
lines_to_keep=$(( reads_per_fastq * 4 ))

while IFS= read -r fastq; do
  base="$(basename "${fastq}")"
  dest="${FASTQ_DIR}/${base}"
  if [[ -s "${dest}" ]]; then
    continue
  fi
  (gzip -dc -- "${fastq}" || true) | head -n "${lines_to_keep}" | gzip -c > "${dest}"
done < <(find -L "${GEX_DIR}" -maxdepth 1 -type f -name "*.fastq.gz" | sort)

R1_FILES=$(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "*_R1_*.fastq.gz" | sort | paste -sd, -)
R2_FILES=$(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "*_R2_*.fastq.gz" | sort | paste -sd, -)

if [[ -z "${R1_FILES}" || -z "${R2_FILES}" ]]; then
  echo "Failed to build downsampled FASTQ list under ${FASTQ_DIR}" >&2
  exit 1
fi

rm -rf "${OUTPREFIX}" "${TMPDIR}"
mkdir -p "${OUTPREFIX}"

/usr/bin/time -v -o "${TIME_LOG}" \
  "${STAR_BIN}" \
    --runThreadN "${THREADS}" \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${R2_FILES}" "${R1_FILES}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${OUTPREFIX}" \
    --outTmpDir "${TMPDIR}" \
    --outSAMtype None \
    --clipAdapterType CellRanger4 \
    --alignEndsType Local \
    --chimSegmentMin 1000000 \
    --soloType CB_UMI_Simple \
    --soloCBstart 1 \
    --soloCBlen 16 \
    --soloUMIstart 17 \
    --soloUMIlen 12 \
    --soloBarcodeReadLength 0 \
    --soloCBwhitelist "${WHITELIST}" \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --soloUMIfiltering "${SOLO_UMI_FILTERING}" \
    --soloUMIdedup "${SOLO_UMI_DEDUP}" \
    --soloMultiMappers "${SOLO_MULTI_MAPPERS}" \
    --soloCbUbRequireTogether no \
    --soloCellFilter "${SOLO_CELL_FILTER}" \
    --soloStrand Unstranded \
    --soloFeatures GeneFull \
    ${SOLO_INLINE_HASH_MODE:+--soloInlineHashMode "${SOLO_INLINE_HASH_MODE}"} \
    "${EXTRA_STAR_ARGS[@]}"

echo "Artifact root: ${ARTIFACT_ROOT}"
echo "STAR output: ${OUTPREFIX}"
echo "Wall-clock log: ${TIME_LOG}"
echo "Downsample target: ${TOTAL_READS} total reads across ${lane_count} lanes (~${reads_per_fastq} reads per FASTQ)"
echo "Solo timings:"
grep "Solo timing:" "${OUTPREFIX}/Log.out" || true
