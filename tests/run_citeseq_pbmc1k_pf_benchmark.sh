#!/usr/bin/env bash
set -euo pipefail
# process_features-only benchmark for the public 10x PBMC 1k TotalSeq-B fixture.
#
# This starts from the official Cell Ranger filtered barcode set, runs only the
# ADT/Antibody Capture feature assignment path, and compares STAR-suite ADT MEX
# against the official Cell Ranger filtered feature-barcode MEX.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PF_DIR="${REPO_ROOT}/core/features/process_features"
ASSIGN="${PF_DIR}/assignBarcodes"
DOWNLOADER="${REPO_ROOT}/scripts/download_10x_citeseq_pbmc1k_fixture.sh"
COMPARE="${REPO_ROOT}/tests/compare_citeseq_mex.py"

FIXTURE="${CITESEQ_PBMC1K_FIXTURE:-/mnt/pikachu/citeseq_pbmc1k_totalseqb_v31}"
OUTDIR="${CITESEQ_PBMC1K_BENCH_OUT:-/mnt/pikachu/citeseq_pbmc1k_pf_benchmark_$(date +%Y%m%d_%H%M%S)}"
WHITELIST="${CITESEQ_NXT_WHITELIST:-/storage/scRNAseq_output/whitelists/3M-february-2018_NXT.txt}"
SOURCE_NAMESPACE="${CITESEQ_PBMC1K_SOURCE_NAMESPACE:-TRU}"
TARGET_NAMESPACE="${CITESEQ_PBMC1K_TARGET_NAMESPACE:-NXT}"
AUTO_DOWNLOAD=0
FORCE_DOWNLOAD=0
MIN_FEATURE_PEARSON="${CITESEQ_MIN_FEATURE_PEARSON:-0.95}"
MIN_BARCODE_PEARSON="${CITESEQ_MIN_BARCODE_PEARSON:-0.90}"
MAX_TOTAL_REL_DELTA="${CITESEQ_MAX_TOTAL_REL_DELTA:-0.25}"

usage() {
  sed -n '2,/^$/s/^# \?//p' "$0"
  cat <<USAGE

Usage:
  tests/run_citeseq_pbmc1k_pf_benchmark.sh [options]

Options:
  --fixture DIR              Staged fixture root. Default: ${FIXTURE}
  --outdir DIR               Benchmark output root. Default: timestamped /mnt/pikachu path
  --whitelist FILE           10x 3' v3 assignment whitelist. Default: ${WHITELIST}
  --source-namespace NS      Namespace of Cell Ranger filtered barcodes. Default: ${SOURCE_NAMESPACE}
  --target-namespace NS      Namespace used for ADT assignment output. Default: ${TARGET_NAMESPACE}
  --download                 Download/stage the fixture if missing (~6.3 GB FASTQ tar)
  --force-download           Re-download/re-stage fixture assets
  --min-feature-pearson X    Default: ${MIN_FEATURE_PEARSON}
  --min-barcode-pearson X    Default: ${MIN_BARCODE_PEARSON}
  --max-total-rel-delta X    Default: ${MAX_TOTAL_REL_DELTA}
  -h, --help                 Show help
USAGE
}

log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"; }
die() { echo "FAIL: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fixture)              FIXTURE="$2"; shift 2 ;;
    --outdir)               OUTDIR="$2"; shift 2 ;;
    --whitelist)            WHITELIST="$2"; shift 2 ;;
    --source-namespace)     SOURCE_NAMESPACE="$2"; shift 2 ;;
    --target-namespace)     TARGET_NAMESPACE="$2"; shift 2 ;;
    --download)             AUTO_DOWNLOAD=1; shift ;;
    --force-download)       AUTO_DOWNLOAD=1; FORCE_DOWNLOAD=1; shift ;;
    --min-feature-pearson)  MIN_FEATURE_PEARSON="$2"; shift 2 ;;
    --min-barcode-pearson)  MIN_BARCODE_PEARSON="$2"; shift 2 ;;
    --max-total-rel-delta)  MAX_TOTAL_REL_DELTA="$2"; shift 2 ;;
    -h|--help)              usage; exit 0 ;;
    *)                      die "Unknown argument: $1" ;;
  esac
done

if [[ "${AUTO_DOWNLOAD}" -eq 1 ]]; then
  dl_args=(--outdir "${FIXTURE}")
  [[ "${FORCE_DOWNLOAD}" -eq 1 ]] && dl_args+=(--force)
  bash "${DOWNLOADER}" "${dl_args[@]}"
fi

ADT_FASTQS="${FIXTURE}/fastqs/adt"
FEATURE_REF="${FIXTURE}/SC3_v3_NextGem_SI_PBMC_CSP_1K_feature_ref.csv"
CR_FILTERED_MEX="${FIXTURE}/cellranger/filtered_feature_bc_matrix"

[[ -d "${ADT_FASTQS}" ]] || die "Missing ADT FASTQ dir: ${ADT_FASTQS}. Re-run with --download."
[[ -s "${FEATURE_REF}" ]] || die "Missing feature ref: ${FEATURE_REF}. Re-run with --download."
[[ -d "${CR_FILTERED_MEX}" ]] || die "Missing Cell Ranger filtered MEX: ${CR_FILTERED_MEX}. Re-run with --download."
if [[ -s "${CR_FILTERED_MEX}/barcodes.tsv.gz" ]]; then
  CR_FILTERED_BARCODES="${CR_FILTERED_MEX}/barcodes.tsv.gz"
elif [[ -s "${CR_FILTERED_MEX}/barcodes.tsv" ]]; then
  CR_FILTERED_BARCODES="${CR_FILTERED_MEX}/barcodes.tsv"
else
  die "Cell Ranger MEX missing barcodes.tsv(.gz)"
fi
[[ -s "${WHITELIST}" ]] || die "Missing whitelist: ${WHITELIST}"
[[ "${SOURCE_NAMESPACE}" == "TRU" || "${SOURCE_NAMESPACE}" == "NXT" ]] || die "--source-namespace must be TRU or NXT"
[[ "${TARGET_NAMESPACE}" == "TRU" || "${TARGET_NAMESPACE}" == "NXT" ]] || die "--target-namespace must be TRU or NXT"

if [[ ! -x "${ASSIGN}" ]]; then
  log "assignBarcodes not found; building process_features tool"
  make -C "${PF_DIR}" assignBarcodes
fi
[[ -x "${ASSIGN}" ]] || die "assignBarcodes build failed: ${ASSIGN}"

mkdir -p "${OUTDIR}"
PF_OUT="${OUTDIR}/process_features_adt"
rm -rf "${PF_OUT}"
CR_FILTERED_BARCODES_FOR_ASSIGN="${CR_FILTERED_BARCODES}"
if [[ "${CR_FILTERED_BARCODES}" == *.gz ]]; then
  CR_FILTERED_BARCODES_FOR_ASSIGN="${OUTDIR}/cellranger_filtered_barcodes.tsv"
  gzip -cd "${CR_FILTERED_BARCODES}" \
    | awk '{ sub(/-[0-9]+$/, "", $1); print $1 }' \
    > "${CR_FILTERED_BARCODES_FOR_ASSIGN}"
else
  CR_FILTERED_BARCODES_FOR_ASSIGN="${OUTDIR}/cellranger_filtered_barcodes.tsv"
  awk '{ sub(/-[0-9]+$/, "", $1); print $1 }' "${CR_FILTERED_BARCODES}" \
    > "${CR_FILTERED_BARCODES_FOR_ASSIGN}"
fi

log "=== CITE-seq PBMC 1k process_features benchmark ==="
log "Fixture:       ${FIXTURE}"
log "ADT FASTQs:    ${ADT_FASTQS}"
log "Feature ref:   ${FEATURE_REF}"
log "CR MEX:        ${CR_FILTERED_MEX}"
log "Whitelist:     ${WHITELIST}"
log "Namespace:     Cell Ranger ${SOURCE_NAMESPACE} -> STAR assignment ${TARGET_NAMESPACE}"
log "Output:        ${OUTDIR}"

"${ASSIGN}" \
  --whitelist "${WHITELIST}" \
  --featurelist "${FEATURE_REF}" \
  --directory "${PF_OUT}" \
  --output-mode adt_mex \
  --filtered_barcodes "${CR_FILTERED_BARCODES_FOR_ASSIGN}" \
  --source_namespace "${SOURCE_NAMESPACE}" \
  --target_namespace "${TARGET_NAMESPACE}" \
  --skip_empty_drops \
  --skip_qc_outputs \
  "${ADT_FASTQS}" \
  -b 16 -u 12 \
  > "${OUTDIR}/assignBarcodes.log" 2>&1

STAR_MEX="$(find "${PF_OUT}" -type f -name protein_quant_summary.json -printf '%h\n' | head -n1)"
[[ -n "${STAR_MEX}" ]] || die "Could not locate STAR-suite ADT MEX under ${PF_OUT}; see ${OUTDIR}/assignBarcodes.log"

compare_args=(
  --cr-mex "${CR_FILTERED_MEX}"
  --star-mex "${STAR_MEX}"
  --feature-type "Antibody Capture"
  --min-feature-pearson "${MIN_FEATURE_PEARSON}"
  --min-barcode-pearson "${MIN_BARCODE_PEARSON}"
  --max-total-rel-delta "${MAX_TOTAL_REL_DELTA}"
  --report-json "${OUTDIR}/citeseq_pbmc1k_compare.json"
  --feature-totals-tsv "${OUTDIR}/citeseq_pbmc1k_feature_totals.tsv"
)
if [[ "${SOURCE_NAMESPACE}" != "${TARGET_NAMESPACE}" ]]; then
  compare_args+=(--translate-star-nxt-tru)
fi

python3 "${COMPARE}" \
  "${compare_args[@]}" \
  | tee "${OUTDIR}/compare.log"

{
  echo "CITE-seq PBMC 1k process_features benchmark"
  echo "Generated (UTC): $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  echo "Fixture: ${FIXTURE}"
  echo "Output: ${OUTDIR}"
  echo "STAR MEX: ${STAR_MEX}"
  echo "Cell Ranger MEX: ${CR_FILTERED_MEX}"
  echo "Whitelist: ${WHITELIST}"
  echo "Namespace: Cell Ranger ${SOURCE_NAMESPACE} -> STAR assignment ${TARGET_NAMESPACE}"
  echo "Compare JSON: ${OUTDIR}/citeseq_pbmc1k_compare.json"
} > "${OUTDIR}/BENCHMARK_SUMMARY.txt"

log "PASS: CITE-seq PBMC 1k process_features benchmark completed"
log "Summary: ${OUTDIR}/BENCHMARK_SUMMARY.txt"
