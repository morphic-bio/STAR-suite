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

OUT_DIR="${FLEX_YREMOVE_SMOKE_OUTDIR:-${SCRIPT_DIR}/flex_yremove_smoke_output_$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${OUT_DIR}/run"
REF_DIR="${OUT_DIR}/reference_split"
PER_SAMPLE_DIR="${RUN_DIR}/per_sample"
SUMMARY_FILE="${OUT_DIR}/SUMMARY.txt"

GENOME_DIR="${FLEX_INDEX:-/storage/flex_filtered_reference/star_index}"
FASTQ_BASE="${YCHROM_FLEX_FASTQ_BASE:-/storage/downsampled/SC2300771}"
SAMPLE_NAME="${FLEX_YREMOVE_SAMPLE_NAME:-SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5}"
THREADS="${FLEX_YREMOVE_THREADS:-1}"

WHITELIST="${FLEX_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SAMPLE_WHITELIST="${FLEX_SAMPLE_WHITELIST:-/storage/SC2300771_filtered_2M/sample_whitelist.tsv}"
PROBE_LIST="${FLEX_PROBE_LIST:-/storage/flex_filtered_reference/filtered_reference/probe_list.txt}"
SAMPLE_PROBES="${FLEX_SAMPLE_PROBES:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
ALLOWED_TAGS="${FLEX_ALLOWED_TAGS:-/storage/SC2300771_filtered_2M/sample_whitelist.tsv}"

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

mapfile -t R2_FILES < <(find "${FASTQ_BASE}" -maxdepth 1 -type f -name "${SAMPLE_NAME}_L00*_R2_001.fastq.gz" | sort)
mapfile -t R1_FILES < <(find "${FASTQ_BASE}" -maxdepth 1 -type f -name "${SAMPLE_NAME}_L00*_R1_001.fastq.gz" | sort)
[[ "${#R2_FILES[@]}" -gt 0 ]] || die "No Flex cDNA FASTQs found under ${FASTQ_BASE}"
[[ "${#R1_FILES[@]}" -eq "${#R2_FILES[@]}" ]] || die "Flex barcode/cDNA lane count mismatch"

rm -rf "${OUT_DIR}"
mkdir -p "${RUN_DIR}" "${REF_DIR}" "${PER_SAMPLE_DIR}"

R2_LIST="$(join_by_comma "${R2_FILES[@]}")"
R1_LIST="$(join_by_comma "${R1_FILES[@]}")"

echo "=== Flex Y-Removal Smoke ==="
echo "Output: ${OUT_DIR}"
echo "Threads: ${THREADS}"
echo "cDNA lanes: ${#R2_FILES[@]}"

"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --soloType CB_UMI_Simple \
  --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${WHITELIST}" \
  --flex yes \
  --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist "${SAMPLE_WHITELIST}" \
  --soloProbeList "${PROBE_LIST}" \
  --soloSampleProbes "${SAMPLE_PROBES}" \
  --soloSampleProbeOffset 68 \
  --soloFlexAllowedTags "${ALLOWED_TAGS}" \
  --soloFlexOutputPrefix "${PER_SAMPLE_DIR}" \
  --limitIObufferSize 50000000 50000000 \
  --outSJtype None \
  --outBAMcompression 6 \
  --soloMultiMappers Rescue \
  --alignIntronMax 500000 \
  --outFilterMismatchNmax 6 \
  --outFilterMismatchNoverReadLmax 1.0 \
  --outFilterMatchNmin 25 \
  --outSAMunmapped None \
  --outFilterMatchNminOverLread 0 \
  --outFilterMultimapNmax 10000 \
  --outFilterMultimapScoreRange 4 \
  --outSAMmultNmax 10000 \
  --winAnchorMultimapNmax 200 \
  --outSAMprimaryFlag AllBestScore \
  --outFilterScoreMin 0 \
  --outFilterScoreMinOverLread 0 \
  --outSAMattributes NH HI AS nM NM GX GN \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCellFilter None \
  --clipAdapterType CellRanger4 \
  --soloFeatures Gene \
  --alignEndsType Local \
  --soloStrand Unstranded \
  --chimSegmentMin 1000000 \
  --soloKeysCompat cr \
  --soloSampleSearchNearby no \
  --readFilesIn "${R2_LIST}" "${R1_LIST}" \
  --outSAMtype BAM Unsorted \
  --emitNoYBAM yes \
  --emitYNoYFastq yes \
  --emitYNoYFastqCompression gz \
  --outFileNamePrefix "${RUN_DIR}/"

Y_COUNT="$(samtools view -c "${RUN_DIR}/Aligned.out_Y.bam")"
[[ "${Y_COUNT}" -gt 0 ]] || die "Expected at least one Y read in Flex Y BAM"
NOY_CHRY="$(samtools view "${RUN_DIR}/Aligned.out_noY.bam" | awk '$3=="chrY"{c++} END{print c+0}')"
[[ "${NOY_CHRY}" -eq 0 ]] || die "Flex noY BAM still contains chrY reads"

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

python3 "${COMPARE_HELPER}" --actual "${ACTUAL_Y_R2}" --reference "${REF_Y_R2}" --ignore-header-comments --label "Flex Y R2"
python3 "${COMPARE_HELPER}" --actual "${ACTUAL_NOY_R2}" --reference "${REF_NOY_R2}" --ignore-header-comments --label "Flex noY R2"
python3 "${COMPARE_HELPER}" --actual "${ACTUAL_Y_R1}" --reference "${REF_Y_R1}" --ignore-header-comments --label "Flex Y R1"
python3 "${COMPARE_HELPER}" --actual "${ACTUAL_NOY_R1}" --reference "${REF_NOY_R1}" --ignore-header-comments --label "Flex noY R1"

[[ -f "${PER_SAMPLE_DIR}/flexfilter_summary.tsv" ]] || die "Missing Flex per-sample summary"

{
  echo "Flex Y-removal smoke: PASS"
  echo "Output: ${OUT_DIR}"
  echo "Threads: ${THREADS}"
  echo "cDNA lanes: ${#R2_FILES[@]}"
  echo "Y BAM reads: ${Y_COUNT}"
  echo "Actual Y R2 FASTQ: ${ACTUAL_Y_R2}"
  echo "Actual noY R2 FASTQ: ${ACTUAL_NOY_R2}"
  echo "Actual Y R1 FASTQ: ${ACTUAL_Y_R1}"
  echo "Actual noY R1 FASTQ: ${ACTUAL_NOY_R1}"
} > "${SUMMARY_FILE}"

echo "PASS: Flex Y-removal smoke"
echo "Summary: ${SUMMARY_FILE}"
