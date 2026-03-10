#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXTERNAL_ENV="${REPO_ROOT}/tests/external_fixtures_env.sh"
if [[ -f "${EXTERNAL_ENV}" ]]; then
  # shellcheck disable=SC1090
  source "${EXTERNAL_ENV}"
fi

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
REMOVE_Y_READS_BIN="${REMOVE_Y_READS_BIN:-${REPO_ROOT}/core/features/yremove_fastq/tools/remove_y_reads/remove_y_reads}"
AGGREGATE_HELPER="${SCRIPT_DIR}/aggregate_remove_y_reads_outputs.py"
COMPARE_HELPER="${SCRIPT_DIR}/compare_fastq_records.py"

OUT_DIR="${SOLO_YREMOVE_SMOKE_OUTDIR:-${SCRIPT_DIR}/solo_yremove_smoke_output_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${OUT_DIR}/run"
REF_DIR="${OUT_DIR}/reference_split"
SUMMARY_FILE="${OUT_DIR}/SUMMARY.txt"

GEX_DIR="${SOLO_YREMOVE_GEX_DIR:-${UCSF_100K_GEX_DIR:-/storage/ucsf-2M/fixtures/ucsf2m_iPSC2_AALG2_100k_pfconfig/GEX/iPSC2_1_AALG2}}"
WHITELIST="${SOLO_YREMOVE_CB_WHITELIST:-${UCSF_100K_CB_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}}"
GENOME_DIR="${SOLO_YREMOVE_GENOME_DIR:-${UCSF_100K_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}}"
THREADS="${SOLO_YREMOVE_THREADS:-1}"

die() {
  echo "ERROR: $*" >&2
  exit 1
}

join_by_comma() {
  local IFS=,
  echo "$*"
}

find_one() {
  local root="$1"
  local pattern="$2"
  mapfile -t matches < <(find "$root" -maxdepth 2 -type f -name "$pattern" | sort)
  [[ "${#matches[@]}" -eq 1 ]] || die "Expected one match for ${pattern} under ${root}, found ${#matches[@]}"
  echo "${matches[0]}"
}

[[ -x "${STAR_BIN}" ]] || die "STAR binary not found: ${STAR_BIN}"
if [[ ! -x "${REMOVE_Y_READS_BIN}" ]]; then
  make -C "${REPO_ROOT}/core/features/yremove_fastq/tools/remove_y_reads" -j8
fi
[[ -x "${REMOVE_Y_READS_BIN}" ]] || die "remove_y_reads binary not found: ${REMOVE_Y_READS_BIN}"

mapfile -t R2_FILES < <(find "${GEX_DIR}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' | sort)
mapfile -t R1_FILES < <(find "${GEX_DIR}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort)

[[ "${#R2_FILES[@]}" -gt 0 ]] || die "No vanilla Solo cDNA FASTQs found under ${GEX_DIR}"
[[ "${#R1_FILES[@]}" -eq "${#R2_FILES[@]}" ]] || die "Mismatched R1/R2 file counts under ${GEX_DIR}"

rm -rf "${OUT_DIR}"
mkdir -p "${RUN_DIR}" "${REF_DIR}"

R2_LIST="$(join_by_comma "${R2_FILES[@]}")"
R1_LIST="$(join_by_comma "${R1_FILES[@]}")"

echo "=== Solo Y-Removal Smoke ==="
echo "Output: ${OUT_DIR}"
echo "Threads: ${THREADS}"
echo "GEX lanes: ${#R2_FILES[@]}"

"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${R2_LIST}" "${R1_LIST}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${RUN_DIR}/" \
  --outSAMtype BAM Unsorted \
  --emitNoYBAM yes \
  --emitYNoYFastq yes \
  --emitYNoYFastqCompression gz \
  --clipAdapterType CellRanger4 \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${WHITELIST}" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique \
  --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no \
  --soloStrand Unstranded \
  --soloFeatures GeneFull

Y_COUNT="$(samtools view -c "${RUN_DIR}/Aligned.out_Y.bam")"
[[ "${Y_COUNT}" -gt 0 ]] || die "Expected at least one Y read in vanilla Solo Y BAM"
NOY_CHRY="$(samtools view "${RUN_DIR}/Aligned.out_noY.bam" | awk '$3=="chrY"{c++} END{print c+0}')"
[[ "${NOY_CHRY}" -eq 0 ]] || die "Vanilla Solo noY BAM still contains chrY reads"

"${REMOVE_Y_READS_BIN}" \
  -y "${RUN_DIR}/Aligned.out_Y.bam" \
  -o "${REF_DIR}" \
  "${R2_FILES[@]}" "${R1_FILES[@]}"

REF_Y_R2="${REF_DIR}/reference_Y_R2_001.fastq.gz"
REF_NOY_R2="${REF_DIR}/reference_noY_R2_001.fastq.gz"
python3 "${AGGREGATE_HELPER}" \
  --split-dir "${REF_DIR}" \
  --output-y "${REF_Y_R2}" \
  --output-noy "${REF_NOY_R2}" \
  "${R2_FILES[@]}"

REF_Y_R1="${REF_DIR}/reference_Y_R1_001.fastq.gz"
REF_NOY_R1="${REF_DIR}/reference_noY_R1_001.fastq.gz"
python3 "${AGGREGATE_HELPER}" \
  --split-dir "${REF_DIR}" \
  --output-y "${REF_Y_R1}" \
  --output-noy "${REF_NOY_R1}" \
  "${R1_FILES[@]}"

ACTUAL_Y_R2="$(find_one "${RUN_DIR}/y_separated" '*_Y_R2_001.fastq.gz')"
ACTUAL_NOY_R2="$(find_one "${RUN_DIR}/y_separated" '*_noY_R2_001.fastq.gz')"
ACTUAL_Y_R1="$(find_one "${RUN_DIR}/y_separated" '*_Y_R1_001.fastq.gz')"
ACTUAL_NOY_R1="$(find_one "${RUN_DIR}/y_separated" '*_noY_R1_001.fastq.gz')"

python3 "${COMPARE_HELPER}" --actual "${ACTUAL_Y_R2}" --reference "${REF_Y_R2}" --ignore-header-comments --label "Solo Y R2"
python3 "${COMPARE_HELPER}" --actual "${ACTUAL_NOY_R2}" --reference "${REF_NOY_R2}" --ignore-header-comments --label "Solo noY R2"
python3 "${COMPARE_HELPER}" --actual "${ACTUAL_Y_R1}" --reference "${REF_Y_R1}" --ignore-header-comments --label "Solo Y R1"
python3 "${COMPARE_HELPER}" --actual "${ACTUAL_NOY_R1}" --reference "${REF_NOY_R1}" --ignore-header-comments --label "Solo noY R1"

[[ -f "${RUN_DIR}/Solo.out/GeneFull/raw/matrix.mtx" ]] || die "Missing vanilla Solo GeneFull raw matrix"
[[ -f "${RUN_DIR}/Solo.out/GeneFull/filtered/matrix.mtx" ]] || die "Missing vanilla Solo GeneFull filtered matrix"

{
  echo "Solo Y-removal smoke: PASS"
  echo "Output: ${OUT_DIR}"
  echo "Threads: ${THREADS}"
  echo "GEX lanes: ${#R2_FILES[@]}"
  echo "Y BAM reads: ${Y_COUNT}"
  echo "Actual Y R2 FASTQ: ${ACTUAL_Y_R2}"
  echo "Actual noY R2 FASTQ: ${ACTUAL_NOY_R2}"
  echo "Actual Y R1 FASTQ: ${ACTUAL_Y_R1}"
  echo "Actual noY R1 FASTQ: ${ACTUAL_NOY_R1}"
} > "${SUMMARY_FILE}"

echo "PASS: Solo Y-removal smoke"
echo "Summary: ${SUMMARY_FILE}"
