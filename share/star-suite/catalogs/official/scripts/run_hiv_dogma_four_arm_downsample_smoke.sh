#!/usr/bin/env bash
set -euo pipefail

STAR_SUITE_DIR="${STAR_SUITE_DIR:-}"
STAR_BIN_DEFAULTED="1"
ATAC_PEAK_MEX_BIN_DEFAULTED="1"
if [[ "${STAR_BIN+x}" == "x" ]]; then
  STAR_BIN_DEFAULTED="0"
fi
if [[ "${ATAC_PEAK_MEX_BIN+x}" == "x" ]]; then
  ATAC_PEAK_MEX_BIN_DEFAULTED="0"
fi
STAR_BIN="${STAR_BIN:-}"
ATAC_PEAK_MEX_BIN="${ATAC_PEAK_MEX_BIN:-}"
FULL_ROOT="${HIV_DOGMA_FULL_ROOT:-}"
DOWNSAMPLE_ROOT="${HIV_DOGMA_DOWNSAMPLE_ROOT:-}"
OUTPUT_ROOT="${HIV_DOGMA_E2E_OUTPUT_ROOT:-}"
RUN_ID="${HIV_DOGMA_E2E_RUN_ID:-star_four_arm_downsample_e2e_100k_$(date -u +%Y%m%dT%H%M%SZ)}"
THREADS="${HIV_DOGMA_E2E_THREADS:-8}"
MAX_READS="${HIV_DOGMA_E2E_MAX_READS:-100000}"
MACS3_QVALUE="${HIV_DOGMA_E2E_MACS3_QVALUE:-}"
RUN_SIGNAC_PROFILE="1"
DRY_RUN="0"

ATAC2GEX="${HIV_DOGMA_E2E_ATAC2GEX:-}"
PROTEIN_FEATURE_REF_DEFAULTED="1"
STATE_FEATURE_REF_DEFAULTED="1"
HIV_STATE_COUNTS_DEFAULTED="1"
if [[ "${HIV_DOGMA_PROTEIN_FEATURE_REF+x}" == "x" ]]; then
  PROTEIN_FEATURE_REF_DEFAULTED="0"
fi
if [[ "${HIV_DOGMA_STATE_FEATURE_REF+x}" == "x" ]]; then
  STATE_FEATURE_REF_DEFAULTED="0"
fi
if [[ "${HIV_DOGMA_HIV_STATE_COUNTS+x}" == "x" ]]; then
  HIV_STATE_COUNTS_DEFAULTED="0"
fi
PROTEIN_FEATURE_REF="${HIV_DOGMA_PROTEIN_FEATURE_REF:-${FULL_ROOT}/star_run_20260615_001327/cr_assign/Protein/adt_yw8/adt/feature_reference.csv}"
STATE_FEATURE_REF="${HIV_DOGMA_STATE_FEATURE_REF:-${FULL_ROOT}/star_run_20260615_001327/cr_assign/Custom/hiv_state_yw8/feature_reference.csv}"
HIV_STATE_COUNTS="${HIV_DOGMA_HIV_STATE_COUNTS:-${FULL_ROOT}/hiv_state_counts.tsv}"

usage() {
  cat <<'EOF'
Usage:
  run_hiv_dogma_four_arm_downsample_smoke.sh [options]

Runs the DOGMA-HIV YW8 four-arm downsample end-to-end smoke through STAR Suite:
GEX + Chromap ATAC + ADT + table-backed HIV state in one STAR process. The
recipe creates or reuses matched first-N FASTQs for every arm, then runs the
STAR-suite DOGMA table-backed smoke and an optional standalone Signac/MACS
BED-profile peak-MEX pass on the generated ATAC sidecar.

Options:
  --star-suite-dir PATH       STAR-suite checkout
  --star-bin PATH             STAR binary
  --atac-peak-mex-bin PATH    standalone star_multiome_atac_peak_mex helper
  --atac-to-gex PATH          ATAC-to-GEX barcode translation table
  --full-root PATH            full DOGMA local STAR input root
  --downsample-root PATH      compact FASTQ fixture root
  --output-root PATH          output root outside this repo
  --run-id ID                 output directory name under output-root
  --threads N                 STAR/Chromap/helper threads
  --max-reads N               first-N read count for every FASTQ
  --macs3-qvalue Q            optional inline FRAG q-value threshold
  --protein-feature-ref PATH  ADT feature reference
  --state-feature-ref PATH    HIV state feature reference
  --hiv-state-counts PATH     table-backed HIV counts TSV
  --skip-signac-profile       skip standalone Signac/MACS BED-profile pass
  --dry-run                   print resolved commands without running
  --help                      show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --star-suite-dir) STAR_SUITE_DIR="$2"; shift 2 ;;
    --star-bin) STAR_BIN="$2"; STAR_BIN_DEFAULTED="0"; shift 2 ;;
    --atac-peak-mex-bin) ATAC_PEAK_MEX_BIN="$2"; ATAC_PEAK_MEX_BIN_DEFAULTED="0"; shift 2 ;;
    --atac-to-gex) ATAC2GEX="$2"; shift 2 ;;
    --full-root) FULL_ROOT="$2"; shift 2 ;;
    --downsample-root) DOWNSAMPLE_ROOT="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --run-id) RUN_ID="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --max-reads) MAX_READS="$2"; shift 2 ;;
    --macs3-qvalue) MACS3_QVALUE="$2"; shift 2 ;;
    --protein-feature-ref) PROTEIN_FEATURE_REF="$2"; PROTEIN_FEATURE_REF_DEFAULTED="0"; shift 2 ;;
    --state-feature-ref) STATE_FEATURE_REF="$2"; STATE_FEATURE_REF_DEFAULTED="0"; shift 2 ;;
    --hiv-state-counts) HIV_STATE_COUNTS="$2"; HIV_STATE_COUNTS_DEFAULTED="0"; shift 2 ;;
    --skip-signac-profile) RUN_SIGNAC_PROFILE="0"; shift ;;
    --dry-run) DRY_RUN="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "${STAR_SUITE_DIR}" ]] || { echo "ERROR: --star-suite-dir is required" >&2; exit 1; }
[[ -n "${FULL_ROOT}" ]] || { echo "ERROR: --full-root is required" >&2; exit 1; }
[[ -n "${DOWNSAMPLE_ROOT}" ]] || { echo "ERROR: --downsample-root is required" >&2; exit 1; }
[[ -n "${OUTPUT_ROOT}" ]] || { echo "ERROR: --output-root is required" >&2; exit 1; }
[[ -n "${ATAC2GEX}" ]] || { echo "ERROR: --atac-to-gex is required" >&2; exit 1; }
if [[ "${STAR_BIN_DEFAULTED}" == "1" ]]; then
  STAR_BIN="${STAR_SUITE_DIR}/core/legacy/source/STAR"
fi
if [[ "${ATAC_PEAK_MEX_BIN_DEFAULTED}" == "1" ]]; then
  ATAC_PEAK_MEX_BIN="${STAR_SUITE_DIR}/core/features/libchromap_contract/star_multiome_atac_peak_mex"
fi
if [[ "${PROTEIN_FEATURE_REF_DEFAULTED}" == "1" ]]; then
  PROTEIN_FEATURE_REF="${FULL_ROOT}/star_run_20260615_001327/cr_assign/Protein/adt_yw8/adt/feature_reference.csv"
fi
if [[ "${STATE_FEATURE_REF_DEFAULTED}" == "1" ]]; then
  STATE_FEATURE_REF="${FULL_ROOT}/star_run_20260615_001327/cr_assign/Custom/hiv_state_yw8/feature_reference.csv"
fi
if [[ "${HIV_STATE_COUNTS_DEFAULTED}" == "1" ]]; then
  HIV_STATE_COUNTS="${FULL_ROOT}/hiv_state_counts.tsv"
fi

STAR_SUITE_DIR="$(realpath -m "${STAR_SUITE_DIR}")"
STAR_BIN="$(realpath -m "${STAR_BIN}")"
ATAC_PEAK_MEX_BIN="$(realpath -m "${ATAC_PEAK_MEX_BIN}")"
FULL_ROOT="$(realpath -m "${FULL_ROOT}")"
DOWNSAMPLE_ROOT="$(realpath -m "${DOWNSAMPLE_ROOT}")"
OUTPUT_ROOT="$(realpath -m "${OUTPUT_ROOT}")"
PROTEIN_FEATURE_REF="$(realpath -m "${PROTEIN_FEATURE_REF}")"
STATE_FEATURE_REF="$(realpath -m "${STATE_FEATURE_REF}")"
HIV_STATE_COUNTS="$(realpath -m "${HIV_STATE_COUNTS}")"

RUN_ROOT="${OUTPUT_ROOT}/${RUN_ID}"
LOG_DIR="${RUN_ROOT}/logs"
RUN_RECORD="${RUN_ROOT}/RUN_RECORD.txt"
STAR_SMOKE="${STAR_SUITE_DIR}/tests/multi_feature/test_hiv_dogma_four_arm_table_smoke.sh"

if [[ ! "${MAX_READS}" =~ ^[0-9]+$ ]] || [[ "${MAX_READS}" -le 0 ]]; then
  echo "ERROR: --max-reads must be a positive integer" >&2
  exit 1
fi

print_cmd() {
  printf '  '
  printf '%q ' "$@"
  printf '\n'
}

SMOKE_CMD=(bash "${STAR_SMOKE}")

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "DOGMA-HIV four-arm downsample smoke dry run"
  echo "RUN_ROOT=${RUN_ROOT}"
  echo "DOWNSAMPLE_ROOT=${DOWNSAMPLE_ROOT}"
  echo "STAR_BIN=${STAR_BIN}"
  echo "MAX_READS=${MAX_READS}"
  echo "THREADS=${THREADS}"
  echo "MACS3_QVALUE=${MACS3_QVALUE}"
  printf 'Downsample fixture: first %s reads per FASTQ from %s/fastq\n' "${MAX_READS}" "${FULL_ROOT}"
  printf 'STAR smoke command:\n'
  print_cmd env \
    HIV_DOGMA_ROOT="${FULL_ROOT}" \
    HIV_DOGMA_FASTQ_ROOT="${DOWNSAMPLE_ROOT}/fastq" \
    HIV_DOGMA_GEX_FASTQ_DIR="${DOWNSAMPLE_ROOT}/gex" \
    HIV_DOGMA_ADT_FASTQ_DIR="${DOWNSAMPLE_ROOT}/adt" \
    HIV_FOUR_ARM_ROOT="${RUN_ROOT}" \
    HIV_DOGMA_PROTEIN_FEATURE_REF="${PROTEIN_FEATURE_REF}" \
    HIV_DOGMA_STATE_FEATURE_REF="${STATE_FEATURE_REF}" \
    HIV_DOGMA_HIV_STATE_COUNTS="${HIV_STATE_COUNTS}" \
    HIV_DOGMA_MATERIALIZE_HIV_TABLE=0 \
    HIV_FOUR_ARM_THREADS="${THREADS}" \
    HIV_FOUR_ARM_MAX_READS="${MAX_READS}" \
    HIV_FOUR_ARM_MACS3_QVALUE="${MACS3_QVALUE}" \
    STAR_BIN="${STAR_BIN}" \
    "${SMOKE_CMD[@]}"
  exit 0
fi

for required in \
  "${STAR_SMOKE}" \
  "${STAR_BIN}" \
  "${FULL_ROOT}/fastq/gex_R1.fastq.gz" \
  "${FULL_ROOT}/fastq/gex_R2.fastq.gz" \
  "${FULL_ROOT}/fastq/adt_R1.fastq.gz" \
  "${FULL_ROOT}/fastq/adt_R2.fastq.gz" \
  "${FULL_ROOT}/fastq/atac_R1.fastq.gz" \
  "${FULL_ROOT}/fastq/atac_R3.fastq.gz" \
  "${FULL_ROOT}/fastq/atac_barcode.fastq.gz" \
  "${PROTEIN_FEATURE_REF}" \
  "${STATE_FEATURE_REF}" \
  "${HIV_STATE_COUNTS}" \
  "${ATAC2GEX}"
do
  [[ -e "${required}" ]] || { echo "ERROR: missing required input ${required}" >&2; exit 1; }
done
if [[ "${RUN_SIGNAC_PROFILE}" == "1" ]]; then
  [[ -x "${ATAC_PEAK_MEX_BIN}" ]] || { echo "ERROR: missing executable ${ATAC_PEAK_MEX_BIN}" >&2; exit 1; }
fi

make_first_n() {
  local src="$1" dst="$2" nreads="$3"
  local lines=$((nreads * 4))
  if [[ -s "${dst}" ]]; then
    return 0
  fi
  echo "writing ${dst}"
  set +o pipefail
  zcat "${src}" | head -n "${lines}" | pigz -p 4 > "${dst}"
  local rc=$?
  set -o pipefail
  return "${rc}"
}

ensure_downsample_fixture() {
  local manifest="${DOWNSAMPLE_ROOT}/MANIFEST.txt"
  if [[ -s "${manifest}" ]] && grep -qx "read_count=${MAX_READS}" "${manifest}"; then
    echo "Using cached DOGMA downsample fixture: ${DOWNSAMPLE_ROOT}"
    return 0
  fi
  rm -rf "${DOWNSAMPLE_ROOT}"
  mkdir -p "${DOWNSAMPLE_ROOT}/fastq" "${DOWNSAMPLE_ROOT}/gex" "${DOWNSAMPLE_ROOT}/adt"
  make_first_n "${FULL_ROOT}/fastq/gex_R1.fastq.gz" "${DOWNSAMPLE_ROOT}/fastq/gex_R1.fastq.gz" "${MAX_READS}"
  make_first_n "${FULL_ROOT}/fastq/gex_R2.fastq.gz" "${DOWNSAMPLE_ROOT}/fastq/gex_R2.fastq.gz" "${MAX_READS}"
  make_first_n "${FULL_ROOT}/fastq/adt_R1.fastq.gz" "${DOWNSAMPLE_ROOT}/fastq/adt_R1.fastq.gz" "${MAX_READS}"
  make_first_n "${FULL_ROOT}/fastq/adt_R2.fastq.gz" "${DOWNSAMPLE_ROOT}/fastq/adt_R2.fastq.gz" "${MAX_READS}"
  make_first_n "${FULL_ROOT}/fastq/atac_R1.fastq.gz" "${DOWNSAMPLE_ROOT}/fastq/atac_R1.fastq.gz" "${MAX_READS}"
  make_first_n "${FULL_ROOT}/fastq/atac_R3.fastq.gz" "${DOWNSAMPLE_ROOT}/fastq/atac_R3.fastq.gz" "${MAX_READS}"
  make_first_n "${FULL_ROOT}/fastq/atac_barcode.fastq.gz" "${DOWNSAMPLE_ROOT}/fastq/atac_barcode.fastq.gz" "${MAX_READS}"
  ln -sfn ../fastq/gex_R1.fastq.gz "${DOWNSAMPLE_ROOT}/gex/YW8_GEX_R1_001.fastq.gz"
  ln -sfn ../fastq/gex_R2.fastq.gz "${DOWNSAMPLE_ROOT}/gex/YW8_GEX_R2_001.fastq.gz"
  ln -sfn ../fastq/adt_R1.fastq.gz "${DOWNSAMPLE_ROOT}/adt/YW8_ADT_R1_001.fastq.gz"
  ln -sfn ../fastq/adt_R2.fastq.gz "${DOWNSAMPLE_ROOT}/adt/YW8_ADT_R2_001.fastq.gz"
  {
    printf 'read_count=%s\n' "${MAX_READS}"
    printf 'source_root=%s\n' "${FULL_ROOT}"
  } > "${manifest}"
}

latest_star_run() {
  find "${RUN_ROOT}" -maxdepth 1 -type d -name 'star_run_*' -printf '%T@ %p\n' \
    | sort -n \
    | tail -1 \
    | cut -d' ' -f2-
}

ensure_downsample_fixture
mkdir -p "${LOG_DIR}"
{
  printf 'date_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'recipe=recipes/run_hiv_dogma_four_arm_downsample_smoke.sh\n'
  printf 'star_suite_dir=%s\n' "${STAR_SUITE_DIR}"
  printf 'star_bin=%s\n' "${STAR_BIN}"
  printf 'atac_peak_mex_bin=%s\n' "${ATAC_PEAK_MEX_BIN}"
  printf 'full_root=%s\n' "${FULL_ROOT}"
  printf 'downsample_root=%s\n' "${DOWNSAMPLE_ROOT}"
  printf 'run_root=%s\n' "${RUN_ROOT}"
  printf 'max_reads=%s\n' "${MAX_READS}"
  printf 'threads=%s\n' "${THREADS}"
  printf 'macs3_qvalue=%s\n' "${MACS3_QVALUE}"
  printf 'protein_feature_ref=%s\n' "${PROTEIN_FEATURE_REF}"
  printf 'state_feature_ref=%s\n' "${STATE_FEATURE_REF}"
  printf 'hiv_state_counts=%s\n' "${HIV_STATE_COUNTS}"
  git -C "${STAR_SUITE_DIR}" rev-parse HEAD 2>/dev/null | sed 's/^/star_suite_commit=/' || true
  git -C "${STAR_SUITE_DIR}" status --short 2>/dev/null | sed 's/^/star_suite_status=/' || true
} > "${RUN_RECORD}"

echo "DOGMA-HIV E2E smoke output root: ${RUN_ROOT}"
env \
  HIV_DOGMA_ROOT="${FULL_ROOT}" \
  HIV_DOGMA_FASTQ_ROOT="${DOWNSAMPLE_ROOT}/fastq" \
  HIV_DOGMA_GEX_FASTQ_DIR="${DOWNSAMPLE_ROOT}/gex" \
  HIV_DOGMA_ADT_FASTQ_DIR="${DOWNSAMPLE_ROOT}/adt" \
  HIV_FOUR_ARM_ROOT="${RUN_ROOT}" \
  HIV_DOGMA_PROTEIN_FEATURE_REF="${PROTEIN_FEATURE_REF}" \
  HIV_DOGMA_STATE_FEATURE_REF="${STATE_FEATURE_REF}" \
  HIV_DOGMA_HIV_STATE_COUNTS="${HIV_STATE_COUNTS}" \
  HIV_DOGMA_MATERIALIZE_HIV_TABLE=0 \
  HIV_FOUR_ARM_THREADS="${THREADS}" \
  HIV_FOUR_ARM_MAX_READS="${MAX_READS}" \
  HIV_FOUR_ARM_MACS3_QVALUE="${MACS3_QVALUE}" \
  STAR_BIN="${STAR_BIN}" \
  "${SMOKE_CMD[@]}" 2>&1 | tee "${LOG_DIR}/hiv_dogma_four_arm_smoke.log"

STAR_RUN="$(latest_star_run)"
if [[ -z "${STAR_RUN}" || ! -d "${STAR_RUN}" ]]; then
  echo "ERROR: could not locate DOGMA star_run_* under ${RUN_ROOT}" >&2
  exit 1
fi

if [[ "${RUN_SIGNAC_PROFILE}" == "1" ]]; then
  SIG="${STAR_RUN}/signac_atac"
  mkdir -p "${SIG}"
  SIGNAC_CMD=(
    "${ATAC_PEAK_MEX_BIN}"
    --sidecar "${STAR_RUN}/atac_fragments.bin"
    --barcode-translate "${ATAC2GEX}"
    --barcode-translate-from-first
    --call-peaks-from-sidecar
    --peak-call-mode macs-bed
    --macs-profile signac-atac
    --peaks "${SIG}/atac_peaks.narrowPeak"
    --summits-out "${SIG}/atac_summits.bed"
    --out-dir "${SIG}/peak_mex"
    --metrics-tsv "${SIG}/atac_metrics.tsv"
    --temp-dir "${SIG}/tmp"
    --keep-intermediates-dir "${SIG}/keep"
    --threads "${THREADS}"
    --force
  )
  printf 'Signac-profile peak-MEX command:\n' | tee "${LOG_DIR}/signac_atac_peak_mex.log"
  print_cmd "${SIGNAC_CMD[@]}" | tee -a "${LOG_DIR}/signac_atac_peak_mex.log"
  "${SIGNAC_CMD[@]}" \
    >"${STAR_RUN}/signac_atac_peak_mex.stdout" \
    2>"${STAR_RUN}/signac_atac_peak_mex.stderr"
  cat "${STAR_RUN}/signac_atac_peak_mex.stderr" | tee -a "${LOG_DIR}/signac_atac_peak_mex.log"
fi

echo "${STAR_RUN}" > "${OUTPUT_ROOT}/LAST_HIV_DOGMA_FOUR_ARM_E2E_STAR_RUN.txt"
echo "DOGMA-HIV E2E smoke completed: ${STAR_RUN}"
