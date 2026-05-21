#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_multiome_cell_call_external_gex_from_arc.sh [options]

Required:
  --out-dir DIR                    Output directory

ARC inputs:
  --arc-run-dir DIR                Cell Ranger ARC run dir; defaults raw MEX to outs/raw_feature_bc_matrix
  --arc-per-barcode-metrics PATH   Path to outs/per_barcode_metrics.csv when not using --arc-run-dir
  --arc-raw-mex-dir DIR            Path to raw_feature_bc_matrix when not using --arc-run-dir

Optional GEX extraction/calling:
  --feature-type TYPE              Feature type to extract from raw MEX (default: Gene Expression)
  --simpleed-mode MODE             simple|full (default: full)
  --umi-min N                      Minimum UMI threshold override for scrna_simpleed
  --sim-n N                        Monte Carlo simulations for scrna_simpleed full mode
  --ed-retain-count N              Retained cell cap for scrna_simpleed
  --lower-testing-bound N          EmptyDrops lower testing bound override
  --ambient-umi-max N              Ambient profile max UMI override
  --use-bootstrap                  Enable CR9-style bootstrap in scrna_simpleed
  --use-raw-pvalue                 Use raw p-value gate instead of FDR
  --fdr X                          FDR threshold override
  --raw-pvalue X                   Raw p-value threshold override
  --include-zero-umis              Pass through to scrna_simpleed
  --gex-rescue-mode MODE           candidate_non_call|raw_p_non_call|fdr_non_call|none

Optional ATAC thresholds:
  --min-peak-region-cutsites N
  --min-peak-region-fragments N
  --min-atac-fragments N
  --min-peak-fraction X

Optional combine settings:
  --combine-mode MODE              gex_atac_rescue|gex_only|union|intersection

Other:
  --scripts-dir DIR                Override STAR-suite scripts dir
  -h, --help                       Show this help
EOF
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

OUT_DIR=""
ARC_RUN_DIR=""
ARC_PER_BARCODE_METRICS=""
ARC_RAW_MEX_DIR=""
FEATURE_TYPE="Gene Expression"
SIMPLEED_MODE="full"
UMI_MIN=""
SIM_N=""
ED_RETAIN_COUNT=""
LOWER_TESTING_BOUND=""
AMBIENT_UMI_MAX=""
USE_BOOTSTRAP=0
USE_RAW_PVALUE=0
FDR=""
RAW_PVALUE=""
INCLUDE_ZERO_UMIS=0
GEX_RESCUE_MODE="candidate_non_call"
MIN_PEAK_REGION_CUTSITES="1"
MIN_PEAK_REGION_FRAGMENTS="0"
MIN_ATAC_FRAGMENTS="0"
MIN_PEAK_FRACTION="0.0"
COMBINE_MODE="gex_atac_rescue"
SCRIPTS_DIR="${SCRIPT_DIR}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2;;
    --arc-run-dir) ARC_RUN_DIR="$2"; shift 2;;
    --arc-per-barcode-metrics) ARC_PER_BARCODE_METRICS="$2"; shift 2;;
    --arc-raw-mex-dir) ARC_RAW_MEX_DIR="$2"; shift 2;;
    --feature-type) FEATURE_TYPE="$2"; shift 2;;
    --simpleed-mode) SIMPLEED_MODE="$2"; shift 2;;
    --umi-min) UMI_MIN="$2"; shift 2;;
    --sim-n) SIM_N="$2"; shift 2;;
    --ed-retain-count) ED_RETAIN_COUNT="$2"; shift 2;;
    --lower-testing-bound) LOWER_TESTING_BOUND="$2"; shift 2;;
    --ambient-umi-max) AMBIENT_UMI_MAX="$2"; shift 2;;
    --use-bootstrap) USE_BOOTSTRAP=1; shift 1;;
    --use-raw-pvalue) USE_RAW_PVALUE=1; shift 1;;
    --fdr) FDR="$2"; shift 2;;
    --raw-pvalue) RAW_PVALUE="$2"; shift 2;;
    --include-zero-umis) INCLUDE_ZERO_UMIS=1; shift 1;;
    --gex-rescue-mode) GEX_RESCUE_MODE="$2"; shift 2;;
    --min-peak-region-cutsites) MIN_PEAK_REGION_CUTSITES="$2"; shift 2;;
    --min-peak-region-fragments) MIN_PEAK_REGION_FRAGMENTS="$2"; shift 2;;
    --min-atac-fragments) MIN_ATAC_FRAGMENTS="$2"; shift 2;;
    --min-peak-fraction) MIN_PEAK_FRACTION="$2"; shift 2;;
    --combine-mode) COMBINE_MODE="$2"; shift 2;;
    --scripts-dir) SCRIPTS_DIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage; exit 2;;
  esac
done

if [[ -z "${OUT_DIR}" ]]; then
  usage
  exit 2
fi

if [[ -n "${ARC_RUN_DIR}" ]]; then
  if [[ -n "${ARC_PER_BARCODE_METRICS}" || -n "${ARC_RAW_MEX_DIR}" ]]; then
    echo "When --arc-run-dir is provided, do not also set --arc-per-barcode-metrics or --arc-raw-mex-dir" >&2
    exit 2
  fi
  ARC_PER_BARCODE_METRICS="${ARC_RUN_DIR}/outs/per_barcode_metrics.csv"
  ARC_RAW_MEX_DIR="${ARC_RUN_DIR}/outs/raw_feature_bc_matrix"
fi

if [[ -z "${ARC_PER_BARCODE_METRICS}" || -z "${ARC_RAW_MEX_DIR}" ]]; then
  echo "Need ARC inputs: either --arc-run-dir or both --arc-per-barcode-metrics and --arc-raw-mex-dir" >&2
  exit 2
fi

if [[ ! -f "${ARC_PER_BARCODE_METRICS}" ]]; then
  echo "Missing ARC per_barcode_metrics.csv: ${ARC_PER_BARCODE_METRICS}" >&2
  exit 2
fi
if [[ ! -d "${ARC_RAW_MEX_DIR}" ]]; then
  echo "Missing ARC raw_feature_bc_matrix directory: ${ARC_RAW_MEX_DIR}" >&2
  exit 2
fi

EXTRACT_SCRIPT="${SCRIPTS_DIR}/extract_cr_feature_type_mex.py"
SIMPLEED_RUNNER="${SCRIPTS_DIR}/run_simpleed_fallback.sh"
MULTIOME_DRIVER="${SCRIPTS_DIR}/run_multiome_cell_call_from_arc.sh"

for req in "${EXTRACT_SCRIPT}" "${SIMPLEED_RUNNER}" "${MULTIOME_DRIVER}"; do
  if [[ ! -f "${req}" ]]; then
    echo "Missing required script: ${req}" >&2
    exit 2
  fi
done

mkdir -p "${OUT_DIR}"

GEX_RAW_MEX_DIR="${OUT_DIR}/gex_raw_mex"
python3 "${EXTRACT_SCRIPT}" \
  --input-mex-dir "${ARC_RAW_MEX_DIR}" \
  --feature-type "${FEATURE_TYPE}" \
  --out-dir "${GEX_RAW_MEX_DIR}" | tee "${OUT_DIR}/extract_gex_mex.log"

GEX_EMPTYDROPS_DIR="${OUT_DIR}/gex_emptydrops"
FILTERED_BARCODES="${GEX_EMPTYDROPS_DIR}/filtered_barcodes.txt"
simpleed_args=(
  --raw-mex "${GEX_RAW_MEX_DIR}"
  --filtered-barcodes "${FILTERED_BARCODES}"
  --force
  --mode "${SIMPLEED_MODE}"
  --out-dir "${GEX_EMPTYDROPS_DIR}"
)
if [[ -n "${UMI_MIN}" ]]; then
  simpleed_args+=(--umi-min "${UMI_MIN}")
fi
if [[ -n "${SIM_N}" ]]; then
  simpleed_args+=(--sim-n "${SIM_N}")
fi
if [[ -n "${ED_RETAIN_COUNT}" ]]; then
  simpleed_args+=(--ed-retain-count "${ED_RETAIN_COUNT}")
fi
if [[ -n "${LOWER_TESTING_BOUND}" ]]; then
  simpleed_args+=(--lower-testing-bound "${LOWER_TESTING_BOUND}")
fi
if [[ -n "${AMBIENT_UMI_MAX}" ]]; then
  simpleed_args+=(--ambient-umi-max "${AMBIENT_UMI_MAX}")
fi
if [[ "${USE_BOOTSTRAP}" -eq 1 ]]; then
  simpleed_args+=(--use-bootstrap)
fi
if [[ "${USE_RAW_PVALUE}" -eq 1 ]]; then
  simpleed_args+=(--use-raw-pvalue)
else
  simpleed_args+=(--use-fdr-gate)
fi
if [[ -n "${FDR}" ]]; then
  simpleed_args+=(--fdr "${FDR}")
fi
if [[ -n "${RAW_PVALUE}" ]]; then
  simpleed_args+=(--raw-pvalue "${RAW_PVALUE}")
fi
if [[ "${INCLUDE_ZERO_UMIS}" -eq 1 ]]; then
  simpleed_args+=(--include-zero-umis)
fi

"${SIMPLEED_RUNNER}" "${simpleed_args[@]}"

EMPTYDROPS_RESULTS="${GEX_EMPTYDROPS_DIR}/EmptyDrops/emptydrops_results.tsv"
driver_args=(
  --arc-per-barcode-metrics "${ARC_PER_BARCODE_METRICS}"
  --filtered-barcodes "${FILTERED_BARCODES}"
  --out-dir "${OUT_DIR}"
  --gex-rescue-mode "${GEX_RESCUE_MODE}"
  --combine-mode "${COMBINE_MODE}"
  --min-peak-region-cutsites "${MIN_PEAK_REGION_CUTSITES}"
  --min-peak-region-fragments "${MIN_PEAK_REGION_FRAGMENTS}"
  --min-atac-fragments "${MIN_ATAC_FRAGMENTS}"
  --min-peak-fraction "${MIN_PEAK_FRACTION}"
)
if [[ -f "${EMPTYDROPS_RESULTS}" ]]; then
  driver_args+=(--emptydrops-results "${EMPTYDROPS_RESULTS}")
fi

"${MULTIOME_DRIVER}" "${driver_args[@]}"

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "arc_run_dir=${ARC_RUN_DIR}"
  echo "arc_per_barcode_metrics=${ARC_PER_BARCODE_METRICS}"
  echo "arc_raw_mex_dir=${ARC_RAW_MEX_DIR}"
  echo "feature_type=${FEATURE_TYPE}"
  echo "gex_raw_mex_dir=${GEX_RAW_MEX_DIR}"
  echo "gex_emptydrops_dir=${GEX_EMPTYDROPS_DIR}"
  echo "filtered_barcodes=${FILTERED_BARCODES}"
  echo "emptydrops_results=${EMPTYDROPS_RESULTS}"
  echo "simpleed_mode=${SIMPLEED_MODE}"
  echo "umi_min=${UMI_MIN}"
  echo "sim_n=${SIM_N}"
  echo "ed_retain_count=${ED_RETAIN_COUNT}"
  echo "lower_testing_bound=${LOWER_TESTING_BOUND}"
  echo "ambient_umi_max=${AMBIENT_UMI_MAX}"
  echo "use_bootstrap=${USE_BOOTSTRAP}"
  echo "use_raw_pvalue=${USE_RAW_PVALUE}"
  echo "fdr=${FDR}"
  echo "raw_pvalue=${RAW_PVALUE}"
  echo "include_zero_umis=${INCLUDE_ZERO_UMIS}"
  echo "gex_rescue_mode=${GEX_RESCUE_MODE}"
  echo "combine_mode=${COMBINE_MODE}"
  echo "min_peak_region_cutsites=${MIN_PEAK_REGION_CUTSITES}"
  echo "min_peak_region_fragments=${MIN_PEAK_REGION_FRAGMENTS}"
  echo "min_atac_fragments=${MIN_ATAC_FRAGMENTS}"
  echo "min_peak_fraction=${MIN_PEAK_FRACTION}"
} > "${OUT_DIR}/external_gex_run_manifest.txt"

echo "External-GEX multiome run complete." >&2
echo "out_dir=${OUT_DIR}" >&2
echo "filtered_barcodes=${FILTERED_BARCODES}" >&2
echo "emptydrops_results=${EMPTYDROPS_RESULTS}" >&2
