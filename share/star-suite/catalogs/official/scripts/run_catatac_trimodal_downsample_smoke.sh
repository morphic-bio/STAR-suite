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
OUTPUT_ROOT="${CATATAC_E2E_OUTPUT_ROOT:-}"
RUN_ID="${CATATAC_E2E_RUN_ID:-catatac_trimodal_e2e_100k_$(date -u +%Y%m%dT%H%M%SZ)}"
THREADS="${CATATAC_E2E_THREADS:-8}"
MAX_READS="${CATATAC_E2E_MAX_READS:-100000}"
MACS3_QVALUE="${CATATAC_E2E_MACS3_QVALUE:-}"
RUN_SIGNAC_PROFILE="1"
DRY_RUN="0"

ATAC2GEX="${CATATAC_E2E_ATAC2GEX:-}"

usage() {
  cat <<'EOF'
Usage:
  run_catatac_trimodal_downsample_smoke.sh [options]

Runs the CAT-ATAC DMSO1 downsample end-to-end smoke through STAR Suite:
GEX + Chromap ATAC + CRISPR guide capture in one STAR process. The recipe also
runs a standalone Signac/MACS BED-profile peak-MEX pass on the generated ATAC
sidecar unless --skip-signac-profile is set.

Options:
  --star-suite-dir PATH      STAR-suite checkout
  --star-bin PATH            STAR binary
  --atac-peak-mex-bin PATH   standalone star_multiome_atac_peak_mex helper
  --atac-to-gex PATH         ATAC-to-GEX barcode translation table
  --output-root PATH         output root outside this repo
  --run-id ID                output directory name under output-root
  --threads N                STAR/Chromap/helper threads
  --max-reads N              downsample read count for GEX and staged ATAC
  --macs3-qvalue Q           optional inline FRAG q-value threshold
  --skip-signac-profile      skip standalone Signac/MACS BED-profile pass
  --dry-run                  print resolved commands without running
  --help                     show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --star-suite-dir) STAR_SUITE_DIR="$2"; shift 2 ;;
    --star-bin) STAR_BIN="$2"; STAR_BIN_DEFAULTED="0"; shift 2 ;;
    --atac-peak-mex-bin) ATAC_PEAK_MEX_BIN="$2"; ATAC_PEAK_MEX_BIN_DEFAULTED="0"; shift 2 ;;
    --atac-to-gex) ATAC2GEX="$2"; shift 2 ;;
    --output-root) OUTPUT_ROOT="$2"; shift 2 ;;
    --run-id) RUN_ID="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --max-reads) MAX_READS="$2"; shift 2 ;;
    --macs3-qvalue) MACS3_QVALUE="$2"; shift 2 ;;
    --skip-signac-profile) RUN_SIGNAC_PROFILE="0"; shift ;;
    --dry-run) DRY_RUN="1"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "${STAR_SUITE_DIR}" ]] || { echo "ERROR: --star-suite-dir is required" >&2; exit 1; }
[[ -n "${OUTPUT_ROOT}" ]] || { echo "ERROR: --output-root is required" >&2; exit 1; }
[[ -n "${ATAC2GEX}" ]] || { echo "ERROR: --atac-to-gex is required" >&2; exit 1; }
if [[ "${STAR_BIN_DEFAULTED}" == "1" ]]; then
  STAR_BIN="${STAR_SUITE_DIR}/core/legacy/source/STAR"
fi
if [[ "${ATAC_PEAK_MEX_BIN_DEFAULTED}" == "1" ]]; then
  ATAC_PEAK_MEX_BIN="${STAR_SUITE_DIR}/core/features/libchromap_contract/star_multiome_atac_peak_mex"
fi

STAR_SUITE_DIR="$(realpath -m "${STAR_SUITE_DIR}")"
STAR_BIN="$(realpath -m "${STAR_BIN}")"
ATAC_PEAK_MEX_BIN="$(realpath -m "${ATAC_PEAK_MEX_BIN}")"
OUTPUT_ROOT="$(realpath -m "${OUTPUT_ROOT}")"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
STAR_SMOKE="${STAR_SUITE_DIR}/tests/test_catatac_trimodal_downsample_smoke.sh"
STAR_RUN="${RUN_DIR}/star_run"
LOG_DIR="${RUN_DIR}/logs"
RUN_RECORD="${RUN_DIR}/RUN_RECORD.txt"

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
SIGNAC_CMD=(
  "${ATAC_PEAK_MEX_BIN}"
  --sidecar "${STAR_RUN}/atac_fragments.bin"
  --barcode-translate "${ATAC2GEX}"
  --barcode-translate-from-first
  --call-peaks-from-sidecar
  --peak-call-mode macs-bed
  --macs-profile signac-atac
  --peaks "${STAR_RUN}/signac_atac/atac_peaks.narrowPeak"
  --summits-out "${STAR_RUN}/signac_atac/atac_summits.bed"
  --out-dir "${STAR_RUN}/signac_atac/peak_mex"
  --metrics-tsv "${STAR_RUN}/signac_atac/atac_metrics.tsv"
  --temp-dir "${STAR_RUN}/signac_atac/tmp"
  --keep-intermediates-dir "${STAR_RUN}/signac_atac/keep"
  --threads "${THREADS}"
  --force
)

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "CAT-ATAC trimodal downsample smoke dry run"
  echo "RUN_DIR=${RUN_DIR}"
  echo "STAR_BIN=${STAR_BIN}"
  echo "MAX_READS=${MAX_READS}"
  echo "THREADS=${THREADS}"
  echo "MACS3_QVALUE=${MACS3_QVALUE}"
  printf 'STAR smoke command:\n'
  print_cmd env \
    CATATAC_TRIMODAL_SMOKE_OUT="${RUN_DIR}" \
    CATATAC_TRIMODAL_INLINE_ATAC_PEAK_MEX=yes \
    CATATAC_TRIMODAL_THREADS="${THREADS}" \
    CATATAC_TRIMODAL_MAX_READS="${MAX_READS}" \
    CATATAC_TRIMODAL_MACS3_QVALUE="${MACS3_QVALUE}" \
    STAR_BIN="${STAR_BIN}" \
    "${SMOKE_CMD[@]}"
  if [[ "${RUN_SIGNAC_PROFILE}" == "1" ]]; then
    printf 'Signac-profile peak-MEX command:\n'
    print_cmd "${SIGNAC_CMD[@]}"
  fi
  exit 0
fi

for required in "${STAR_SMOKE}" "${STAR_BIN}" "${ATAC2GEX}"; do
  [[ -e "${required}" ]] || { echo "ERROR: missing required input ${required}" >&2; exit 1; }
done
if [[ "${RUN_SIGNAC_PROFILE}" == "1" ]]; then
  [[ -x "${ATAC_PEAK_MEX_BIN}" ]] || { echo "ERROR: missing executable ${ATAC_PEAK_MEX_BIN}" >&2; exit 1; }
fi

echo "CAT-ATAC E2E smoke output: ${RUN_DIR}"
env \
  CATATAC_TRIMODAL_SMOKE_OUT="${RUN_DIR}" \
  CATATAC_TRIMODAL_INLINE_ATAC_PEAK_MEX=yes \
  CATATAC_TRIMODAL_THREADS="${THREADS}" \
  CATATAC_TRIMODAL_MAX_READS="${MAX_READS}" \
  CATATAC_TRIMODAL_MACS3_QVALUE="${MACS3_QVALUE}" \
  STAR_BIN="${STAR_BIN}" \
  "${SMOKE_CMD[@]}"

mkdir -p "${LOG_DIR}"
{
  printf 'date_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'recipe=recipes/run_catatac_trimodal_downsample_smoke.sh\n'
  printf 'star_suite_dir=%s\n' "${STAR_SUITE_DIR}"
  printf 'star_bin=%s\n' "${STAR_BIN}"
  printf 'atac_peak_mex_bin=%s\n' "${ATAC_PEAK_MEX_BIN}"
  printf 'run_dir=%s\n' "${RUN_DIR}"
  printf 'star_run=%s\n' "${STAR_RUN}"
  printf 'max_reads=%s\n' "${MAX_READS}"
  printf 'threads=%s\n' "${THREADS}"
  printf 'macs3_qvalue=%s\n' "${MACS3_QVALUE}"
  git -C "${STAR_SUITE_DIR}" rev-parse HEAD 2>/dev/null | sed 's/^/star_suite_commit=/' || true
  git -C "${STAR_SUITE_DIR}" status --short 2>/dev/null | sed 's/^/star_suite_status=/' || true
} > "${RUN_RECORD}"

if [[ "${RUN_SIGNAC_PROFILE}" == "1" ]]; then
  mkdir -p "${STAR_RUN}/signac_atac"
  printf 'Signac-profile peak-MEX command:\n' | tee "${LOG_DIR}/signac_atac_peak_mex.log"
  print_cmd "${SIGNAC_CMD[@]}" | tee -a "${LOG_DIR}/signac_atac_peak_mex.log"
  "${SIGNAC_CMD[@]}" \
    >"${STAR_RUN}/signac_atac_peak_mex.stdout" \
    2>"${STAR_RUN}/signac_atac_peak_mex.stderr"
  cat "${STAR_RUN}/signac_atac_peak_mex.stderr" | tee -a "${LOG_DIR}/signac_atac_peak_mex.log"
fi

echo "${STAR_RUN}" > "${OUTPUT_ROOT}/LAST_CATATAC_TRIMODAL_E2E_STAR_RUN.txt"
echo "CAT-ATAC E2E smoke completed: ${STAR_RUN}"
