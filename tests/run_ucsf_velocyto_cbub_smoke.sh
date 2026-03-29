#!/usr/bin/env bash
# UCSF 100K host-local smoke for the shared packedReadInfo lifetime:
# - production-style corrected perturb surface
# - GeneFull + Velocyto
# - unsorted BAM with CB/UB tag injection requested
# - Y/noY BAM emission enabled
#
# This catches regressions where packedReadInfo is freed before either:
# - Velocyto consumes it during later Solo feature processing
# - BAM CB/UB tag injection consumes it after Solo processing
#
# IMPORTANT:
# - this smoke depends on protected UCSF inputs / staged outputs
# - generated artifacts must remain host-local and must not be redistributed

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
EXTERNAL_ENV="${REPO_ROOT}/tests/external_fixtures_env.sh"
if [[ -f "${EXTERNAL_ENV}" ]]; then
  # shellcheck disable=SC1090
  source "${EXTERNAL_ENV}"
fi

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR.release}"
BATCH_RUNNER="${REPO_ROOT}/scripts/run_ucsf_perturb_yremove_batch.sh"
PREPARE_VELOCYTO_MEX="${REPO_ROOT}/scripts/prepare_velocyto_mex.py"

SAMPLE="${UCSF_VELOCYTO_CBUB_SAMPLE:-EBs2_2}"
THREADS="${UCSF_VELOCYTO_CBUB_THREADS:-24}"
READS="${UCSF_VELOCYTO_CBUB_READS:-100000}"
DOWNSAMPLE_SEED="${UCSF_VELOCYTO_CBUB_DOWNSAMPLE_SEED:-1}"

DATASET_ROOT="${UCSF_VELOCYTO_CBUB_DATASET_ROOT:-/mnt/pikachu/ucsf-perturb-seq-corrected}"
FEATURE_REF="${UCSF_VELOCYTO_CBUB_FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
GENOME_DIR="${UCSF_VELOCYTO_CBUB_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
SOLO_CB_WHITELIST="${UCSF_VELOCYTO_CBUB_SOLO_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
CR_WHITELIST="${UCSF_VELOCYTO_CBUB_CR_WHITELIST:-${SOLO_CB_WHITELIST}}"

SOURCE_SAMPLE_ROOT="${UCSF_VELOCYTO_CBUB_SOURCE_SAMPLE_ROOT:-}"
GENERATED_STAGE_OUT_ROOT="${UCSF_VELOCYTO_CBUB_STAGE_OUT_ROOT:-/tmp/ucsf_velocyto_cbub_stage_${SAMPLE}_${READS}}"

STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${UCSF_VELOCYTO_CBUB_OUTDIR:-/tmp/ucsf_velocyto_cbub_smoke_${STAMP}}"
RUN_DIR="${OUT_DIR}/run"
SUMMARY_FILE="${OUT_DIR}/SUMMARY.txt"

skip() {
  echo "SKIP: $*"
  exit 0
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

join_by_comma() {
  local IFS=,
  echo "$*"
}

require_generated_stage() {
  [[ -x "${BATCH_RUNNER}" ]] || skip "Missing batch runner: ${BATCH_RUNNER}"
  [[ -d "${DATASET_ROOT}/${SAMPLE}/GEX" ]] || skip "Missing corrected GEX dir: ${DATASET_ROOT}/${SAMPLE}/GEX"
  [[ -d "${DATASET_ROOT}/${SAMPLE}/guides" ]] || skip "Missing corrected guide dir: ${DATASET_ROOT}/${SAMPLE}/guides"
  [[ -f "${FEATURE_REF}" ]] || skip "Missing feature reference: ${FEATURE_REF}"
  [[ -d "${GENOME_DIR}" ]] || skip "Missing genome dir: ${GENOME_DIR}"
  [[ -f "${SOLO_CB_WHITELIST}" ]] || skip "Missing Solo CB whitelist: ${SOLO_CB_WHITELIST}"
  [[ -f "${CR_WHITELIST}" ]] || skip "Missing CR whitelist: ${CR_WHITELIST}"
  [[ -x "${STAR_BIN}" ]] || skip "Missing STAR binary: ${STAR_BIN}"
  command -v seqtk >/dev/null 2>&1 || skip "seqtk not found in PATH"
  command -v pigz >/dev/null 2>&1 || skip "pigz not found in PATH"

  rm -rf "${GENERATED_STAGE_OUT_ROOT}"
  mkdir -p "${GENERATED_STAGE_OUT_ROOT}"

  echo "Preparing staged 100K UCSF corrected fixture via batch-runner dry-run..."
  (
    export DOWNSAMPLE_SEED
    "${BATCH_RUNNER}" \
      --out-root "${GENERATED_STAGE_OUT_ROOT}" \
      --dataset-root "${DATASET_ROOT}" \
      --samples "${SAMPLE}" \
      --feature-ref "${FEATURE_REF}" \
      --genome-dir "${GENOME_DIR}" \
      --solo-cb-whitelist "${SOLO_CB_WHITELIST}" \
      --cr-whitelist "${CR_WHITELIST}" \
      --threads "${THREADS}" \
      --star-bin "${STAR_BIN}" \
      --downsample-reads "${READS}" \
      --dry-run
  )
  SOURCE_SAMPLE_ROOT="${GENERATED_STAGE_OUT_ROOT}/samples/${SAMPLE}"
}

[[ -x "${STAR_BIN}" ]] || skip "STAR binary not found: ${STAR_BIN}"
[[ -f "${PREPARE_VELOCYTO_MEX}" ]] || skip "Missing helper: ${PREPARE_VELOCYTO_MEX}"
command -v samtools >/dev/null 2>&1 || skip "samtools not found in PATH"
command -v python3 >/dev/null 2>&1 || skip "python3 not found in PATH"

if [[ -z "${SOURCE_SAMPLE_ROOT}" || ! -d "${SOURCE_SAMPLE_ROOT}/staged_input" || ! -f "${SOURCE_SAMPLE_ROOT}/pf_multi_config.csv" ]]; then
  require_generated_stage
fi

STAGED_ROOT="${SOURCE_SAMPLE_ROOT}/staged_input"
PF_CONFIG="${SOURCE_SAMPLE_ROOT}/pf_multi_config.csv"
STAGED_GEX_DIR="${STAGED_ROOT}/GEX/${SAMPLE}"
STAGED_GUIDE_DIR="${STAGED_ROOT}/guides/${SAMPLE}"

[[ -d "${STAGED_GEX_DIR}" ]] || die "Missing staged GEX dir: ${STAGED_GEX_DIR}"
[[ -d "${STAGED_GUIDE_DIR}" ]] || die "Missing staged guide dir: ${STAGED_GUIDE_DIR}"
[[ -f "${PF_CONFIG}" ]] || die "Missing staged pf_multi_config.csv: ${PF_CONFIG}"
[[ -f "${FEATURE_REF}" ]] || die "Missing feature reference: ${FEATURE_REF}"
[[ -d "${GENOME_DIR}" ]] || die "Missing genome dir: ${GENOME_DIR}"
[[ -f "${SOLO_CB_WHITELIST}" ]] || die "Missing Solo CB whitelist: ${SOLO_CB_WHITELIST}"
[[ -f "${CR_WHITELIST}" ]] || die "Missing CR whitelist: ${CR_WHITELIST}"

mapfile -t GEX_R2_FILES < <(find -L "${STAGED_GEX_DIR}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' | sort)
mapfile -t GEX_R1_FILES < <(find -L "${STAGED_GEX_DIR}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort)
mapfile -t GUIDE_R2_FILES < <(find -L "${STAGED_GUIDE_DIR}" -maxdepth 1 -type f -name '*_R2_001.fastq.gz' | sort)
mapfile -t GUIDE_R1_FILES < <(find -L "${STAGED_GUIDE_DIR}" -maxdepth 1 -type f -name '*_R1_001.fastq.gz' | sort)

[[ "${#GEX_R2_FILES[@]}" -gt 0 ]] || die "No staged GEX R2 FASTQs under ${STAGED_GEX_DIR}"
[[ "${#GUIDE_R2_FILES[@]}" -gt 0 ]] || die "No staged guide R2 FASTQs under ${STAGED_GUIDE_DIR}"
[[ "${#GEX_R2_FILES[@]}" -eq "${#GEX_R1_FILES[@]}" ]] || die "Staged GEX R1/R2 mismatch"
[[ "${#GUIDE_R2_FILES[@]}" -eq "${#GUIDE_R1_FILES[@]}" ]] || die "Staged guide R1/R2 mismatch"

ALL_R2=("${GEX_R2_FILES[@]}" "${GUIDE_R2_FILES[@]}")
ALL_R1=("${GEX_R1_FILES[@]}" "${GUIDE_R1_FILES[@]}")
R2_LIST="$(join_by_comma "${ALL_R2[@]}")"
R1_LIST="$(join_by_comma "${ALL_R1[@]}")"

rm -rf "${OUT_DIR}"
mkdir -p "${RUN_DIR}"

echo "=== UCSF Velocyto + CB/UB smoke ==="
echo "Sample: ${SAMPLE}"
echo "Source sample root: ${SOURCE_SAMPLE_ROOT}"
echo "Output: ${OUT_DIR}"
echo "Threads: ${THREADS}"

"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${R2_LIST}" "${R1_LIST}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${RUN_DIR}/" \
  --outTmpDir "${OUT_DIR}/tmp" \
  --outSAMtype BAM Unsorted \
  --outSAMattributes NH HI AS nM NM GX GN CB UB \
  --emitNoYBAM yes \
  --emitYNoYFastq yes \
  --emitYNoYFastqCompression gz \
  --clipAdapterType CellRanger4 \
  --clip3pPolyG yes \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 \
  --soloCBlen 16 \
  --soloUMIstart 17 \
  --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${SOLO_CB_WHITELIST}" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique \
  --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no \
  --soloStrand Forward \
  --soloFeatures GeneFull Velocyto \
  --soloCrGexFeature genefull \
  --soloCrMultimapRescue yes \
  --pfMultiConfig "${PF_CONFIG}" \
  --crFeatureRef "${FEATURE_REF}" \
  --crChemistry auto \
  --crOutputChemistry TRU \
  --crWhitelist "${CR_WHITELIST}" \
  --crMinUmi 3 \
  --crAssignMaxHamming 1 \
  --crAssignFeatureOffset -1 \
  --crAssignLimitSearch -1 \
  --crAssignMinCounts 0 \
  --crAssignMaxBarcodeMismatches 5 \
  --crAssignFeatureN 0 \
  --crAssignBarcodeN 1 \
  --crAssignConsumerThreads -1 \
  --crAssignSearchThreads 1 \
  --crAssignSkipQcOutputs 1 \
  --dynamicThreadInterface 1 \
  --dynamicThreadConstMapPermits "${THREADS}" \
  --dynamicThreadTelemetry 1

python3 "${PREPARE_VELOCYTO_MEX}" --run-dir "${RUN_DIR}"

for path in \
  "${RUN_DIR}/Aligned.out.bam" \
  "${RUN_DIR}/Aligned.out_Y.bam" \
  "${RUN_DIR}/Aligned.out_noY.bam" \
  "${RUN_DIR}/Solo.out/Velocyto/raw/spliced.mtx" \
  "${RUN_DIR}/Solo.out/Velocyto/raw/unspliced.mtx" \
  "${RUN_DIR}/Solo.out/Velocyto/raw/ambiguous.mtx" \
  "${RUN_DIR}/outs/velocyto_feature_bc_matrix_manifest.json" \
  "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix/matrix.mtx.gz" \
  "${RUN_DIR}/outs/filtered_velocyto_feature_bc_matrix/matrix.mtx.gz"; do
  [[ -f "${path}" ]] || die "Missing expected output: ${path}"
done

RAW_UNSPLICED_NNZ="$(awk 'NR==3 {print $3}' "${RUN_DIR}/Solo.out/Velocyto/raw/unspliced.mtx")"
[[ -n "${RAW_UNSPLICED_NNZ}" ]] || die "Failed to read unspliced raw nnz"
if [[ "${RAW_UNSPLICED_NNZ}" -le 0 ]]; then
  die "Velocyto raw unspliced nnz is zero"
fi

TAG_BAM="${RUN_DIR}/Aligned.out_noY.bam"
TOTAL_RECORDS="$(samtools view -c "${TAG_BAM}")"
CB_RECORDS="$(samtools view "${TAG_BAM}" | grep -c 'CB:Z:' || true)"
UB_RECORDS="$(samtools view "${TAG_BAM}" | grep -c 'UB:Z:' || true)"
CBUB_RECORDS="$(samtools view "${TAG_BAM}" | grep 'CB:Z:' | grep -c 'UB:Z:' || true)"

[[ "${TOTAL_RECORDS}" -gt 0 ]] || die "No records found in ${TAG_BAM}"
[[ "${CB_RECORDS}" -gt 0 ]] || die "No CB tags found in ${TAG_BAM}"
[[ "${UB_RECORDS}" -gt 0 ]] || die "No UB tags found in ${TAG_BAM}"

readarray -t METRICS < <(python3 - <<'PY' "${RUN_DIR}/outs/velocyto_feature_bc_matrix_manifest.json"
import json
import sys
from pathlib import Path

manifest = json.loads(Path(sys.argv[1]).read_text())
print(f"RAW_BARCODES={manifest['raw']['barcodes']}")
print(f"FILTERED_BARCODES={manifest['filtered']['barcodes']}")
print(f"RAW_NNZ={manifest['raw']['nnz_total']}")
print(f"FILTERED_NNZ={manifest['filtered']['nnz_total']}")
PY
)

declare -A STATS=()
for line in "${METRICS[@]}"; do
  key="${line%%=*}"
  value="${line#*=}"
  STATS["${key}"]="${value}"
done

FIRST_TAGGED_RECORD="$(samtools view "${TAG_BAM}" | grep 'CB:Z:.*UB:Z:' | head -n 1 || true)"

{
  echo "UCSF Velocyto + CB/UB smoke: PASS"
  echo "Sample: ${SAMPLE}"
  echo "Source sample root: ${SOURCE_SAMPLE_ROOT}"
  echo "Staged GEX dir: ${STAGED_GEX_DIR}"
  echo "Staged guide dir: ${STAGED_GUIDE_DIR}"
  echo "pfMultiConfig: ${PF_CONFIG}"
  echo "STAR binary: ${STAR_BIN}"
  echo "Threads: ${THREADS}"
  echo "Output: ${OUT_DIR}"
  echo "Tag BAM: ${TAG_BAM}"
  echo "Total BAM records: ${TOTAL_RECORDS}"
  echo "CB-tagged BAM records: ${CB_RECORDS}"
  echo "UB-tagged BAM records: ${UB_RECORDS}"
  echo "CB+UB-tagged BAM records: ${CBUB_RECORDS}"
  echo "Raw velocyto barcodes: ${STATS[RAW_BARCODES]}"
  echo "Filtered velocyto barcodes: ${STATS[FILTERED_BARCODES]}"
  echo "Raw velocyto nnz: ${STATS[RAW_NNZ]}"
  echo "Filtered velocyto nnz: ${STATS[FILTERED_NNZ]}"
  echo "Raw unspliced nnz: ${RAW_UNSPLICED_NNZ}"
  if [[ -n "${FIRST_TAGGED_RECORD}" ]]; then
    echo "Example tagged record:"
    echo "${FIRST_TAGGED_RECORD}"
  fi
} > "${SUMMARY_FILE}"

cat "${SUMMARY_FILE}"
