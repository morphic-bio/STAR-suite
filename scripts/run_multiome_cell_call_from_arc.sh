#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_multiome_cell_call_from_arc.sh [options]

Required:
  --out-dir DIR                    Output directory for all artifacts
  --filtered-barcodes PATH         GEX filtered_barcodes.txt path

ARC inputs (choose one):
  --arc-run-dir DIR                Cell Ranger ARC run dir containing outs/per_barcode_metrics.csv
  --arc-per-barcode-metrics PATH   Path to ARC outs/per_barcode_metrics.csv
  --arc-table PATH                 Prebuilt arc_barcode_table.tsv

Optional GEX inputs:
  --emptydrops-results PATH        Optional EmptyDrops/emptydrops_results.tsv
  --barcode-suffix SUFFIX          Append suffix when upstream barcodes lack one (example: -1)
  --gex-rescue-mode MODE           candidate_non_call|raw_p_non_call|fdr_non_call|none

Optional ATAC thresholds:
  --min-peak-region-cutsites N     Default: 1
  --min-peak-region-fragments N    Default: 0
  --min-atac-fragments N           Default: 0
  --min-peak-fraction X            Default: 0.0

Optional combiner settings:
  --combine-mode MODE              gex_atac_rescue|gex_only|union|intersection

Other options:
  --skip-compare                   Do not run comparison_vs_arc outputs
  --libscrna-tools-dir DIR         Override libscrna bin dir
  -h, --help                       Show this help

Outputs written under --out-dir:
  arc_barcode_table.tsv
  gex_evidence.tsv
  atac_evidence.tsv
  multiome_calls.tsv
  comparison_vs_arc.tsv
  comparison_summary.tsv
  run_manifest.txt
EOF
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

OUT_DIR=""
FILTERED_BARCODES=""
ARC_RUN_DIR=""
ARC_PER_BARCODE_METRICS=""
ARC_TABLE=""
EMPTYDROPS_RESULTS=""
BARCODE_SUFFIX=""
GEX_RESCUE_MODE="candidate_non_call"
MIN_PEAK_REGION_CUTSITES="1"
MIN_PEAK_REGION_FRAGMENTS="0"
MIN_ATAC_FRAGMENTS="0"
MIN_PEAK_FRACTION="0.0"
COMBINE_MODE="gex_atac_rescue"
RUN_COMPARE=1
LIBSCRNA_TOOLS_DIR="${REPO_ROOT}/core/features/libscrna/bin"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2;;
    --filtered-barcodes) FILTERED_BARCODES="$2"; shift 2;;
    --arc-run-dir) ARC_RUN_DIR="$2"; shift 2;;
    --arc-per-barcode-metrics) ARC_PER_BARCODE_METRICS="$2"; shift 2;;
    --arc-table) ARC_TABLE="$2"; shift 2;;
    --emptydrops-results) EMPTYDROPS_RESULTS="$2"; shift 2;;
    --barcode-suffix) BARCODE_SUFFIX="$2"; shift 2;;
    --gex-rescue-mode) GEX_RESCUE_MODE="$2"; shift 2;;
    --min-peak-region-cutsites) MIN_PEAK_REGION_CUTSITES="$2"; shift 2;;
    --min-peak-region-fragments) MIN_PEAK_REGION_FRAGMENTS="$2"; shift 2;;
    --min-atac-fragments) MIN_ATAC_FRAGMENTS="$2"; shift 2;;
    --min-peak-fraction) MIN_PEAK_FRACTION="$2"; shift 2;;
    --combine-mode) COMBINE_MODE="$2"; shift 2;;
    --skip-compare) RUN_COMPARE=0; shift 1;;
    --libscrna-tools-dir) LIBSCRNA_TOOLS_DIR="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage; exit 2;;
  esac
done

if [[ -z "${OUT_DIR}" || -z "${FILTERED_BARCODES}" ]]; then
  usage
  exit 2
fi

if [[ -n "${ARC_RUN_DIR}" && -n "${ARC_PER_BARCODE_METRICS}" ]]; then
  echo "Specify only one of --arc-run-dir or --arc-per-barcode-metrics" >&2
  exit 2
fi

if [[ -n "${ARC_RUN_DIR}" ]]; then
  ARC_PER_BARCODE_METRICS="${ARC_RUN_DIR}/outs/per_barcode_metrics.csv"
fi

if [[ -z "${ARC_TABLE}" && -z "${ARC_PER_BARCODE_METRICS}" ]]; then
  echo "One of --arc-table, --arc-per-barcode-metrics, or --arc-run-dir is required" >&2
  exit 2
fi

if [[ ! -f "${FILTERED_BARCODES}" ]]; then
  echo "Missing filtered_barcodes.txt: ${FILTERED_BARCODES}" >&2
  exit 2
fi

if [[ -n "${ARC_PER_BARCODE_METRICS}" && ! -f "${ARC_PER_BARCODE_METRICS}" ]]; then
  echo "Missing per_barcode_metrics.csv: ${ARC_PER_BARCODE_METRICS}" >&2
  exit 2
fi

if [[ -n "${ARC_TABLE}" && ! -f "${ARC_TABLE}" ]]; then
  echo "Missing arc_barcode_table.tsv: ${ARC_TABLE}" >&2
  exit 2
fi

if [[ -n "${EMPTYDROPS_RESULTS}" && ! -f "${EMPTYDROPS_RESULTS}" ]]; then
  echo "Missing emptydrops_results.tsv: ${EMPTYDROPS_RESULTS}" >&2
  exit 2
fi

ARC_TABLE_TOOL="${LIBSCRNA_TOOLS_DIR}/scrna_arc_barcode_table"
GEX_EVIDENCE_TOOL="${LIBSCRNA_TOOLS_DIR}/scrna_build_gex_evidence"
ATAC_EVIDENCE_TOOL="${LIBSCRNA_TOOLS_DIR}/scrna_build_atac_evidence"
COMBINE_TOOL="${LIBSCRNA_TOOLS_DIR}/scrna_multiome_combine"
COMPARE_TOOL="${LIBSCRNA_TOOLS_DIR}/scrna_compare_multiome_calls"

need_build=0
for tool in \
  "${ARC_TABLE_TOOL}" \
  "${GEX_EVIDENCE_TOOL}" \
  "${ATAC_EVIDENCE_TOOL}" \
  "${COMBINE_TOOL}" \
  "${COMPARE_TOOL}"; do
  if [[ ! -x "${tool}" ]]; then
    need_build=1
  fi
done

if [[ "${need_build}" -eq 1 ]]; then
  echo "Building libscrna tools..." >&2
  make -C "${REPO_ROOT}/core/features/libscrna" tools
fi

mkdir -p "${OUT_DIR}"

ARC_TABLE_OUT="${OUT_DIR}/arc_barcode_table.tsv"
if [[ -n "${ARC_TABLE}" ]]; then
  cp "${ARC_TABLE}" "${ARC_TABLE_OUT}"
else
  "${ARC_TABLE_TOOL}" \
    --per-barcode-metrics "${ARC_PER_BARCODE_METRICS}" \
    --out "${ARC_TABLE_OUT}"
fi

GEX_EVIDENCE_OUT="${OUT_DIR}/gex_evidence.tsv"
gex_args=(
  --arc-table "${ARC_TABLE_OUT}"
  --filtered-barcodes "${FILTERED_BARCODES}"
  --out "${GEX_EVIDENCE_OUT}"
  --rescue-mode "${GEX_RESCUE_MODE}"
)
if [[ -n "${EMPTYDROPS_RESULTS}" ]]; then
  gex_args+=(--emptydrops-results "${EMPTYDROPS_RESULTS}")
fi
if [[ -n "${BARCODE_SUFFIX}" ]]; then
  gex_args+=(--barcode-suffix "${BARCODE_SUFFIX}")
fi
"${GEX_EVIDENCE_TOOL}" "${gex_args[@]}"

ATAC_EVIDENCE_OUT="${OUT_DIR}/atac_evidence.tsv"
"${ATAC_EVIDENCE_TOOL}" \
  --arc-table "${ARC_TABLE_OUT}" \
  --out "${ATAC_EVIDENCE_OUT}" \
  --min-peak-region-cutsites "${MIN_PEAK_REGION_CUTSITES}" \
  --min-peak-region-fragments "${MIN_PEAK_REGION_FRAGMENTS}" \
  --min-atac-fragments "${MIN_ATAC_FRAGMENTS}" \
  --min-peak-fraction "${MIN_PEAK_FRACTION}"

MULTIOME_CALLS_OUT="${OUT_DIR}/multiome_calls.tsv"
"${COMBINE_TOOL}" \
  --arc-table "${ARC_TABLE_OUT}" \
  --gex-evidence "${GEX_EVIDENCE_OUT}" \
  --atac-evidence "${ATAC_EVIDENCE_OUT}" \
  --gex-rescue-col gex_rescue_eligible \
  --atac-low-targeting-col atac_low_targeting \
  --mode "${COMBINE_MODE}" \
  --out "${MULTIOME_CALLS_OUT}"

COMPARE_OUT="${OUT_DIR}/comparison_vs_arc.tsv"
SUMMARY_OUT="${OUT_DIR}/comparison_summary.tsv"
if [[ "${RUN_COMPARE}" -eq 1 ]]; then
  "${COMPARE_TOOL}" \
    --arc-table "${ARC_TABLE_OUT}" \
    --calls "${MULTIOME_CALLS_OUT}" \
    --out "${COMPARE_OUT}" \
    --summary "${SUMMARY_OUT}"
fi

MANIFEST_OUT="${OUT_DIR}/run_manifest.txt"
{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "repo_root=${REPO_ROOT}"
  echo "libscrna_tools_dir=${LIBSCRNA_TOOLS_DIR}"
  echo "arc_run_dir=${ARC_RUN_DIR}"
  echo "arc_per_barcode_metrics=${ARC_PER_BARCODE_METRICS}"
  echo "arc_table_input=${ARC_TABLE}"
  echo "arc_table_output=${ARC_TABLE_OUT}"
  echo "filtered_barcodes=${FILTERED_BARCODES}"
  echo "emptydrops_results=${EMPTYDROPS_RESULTS}"
  echo "barcode_suffix=${BARCODE_SUFFIX}"
  echo "gex_rescue_mode=${GEX_RESCUE_MODE}"
  echo "min_peak_region_cutsites=${MIN_PEAK_REGION_CUTSITES}"
  echo "min_peak_region_fragments=${MIN_PEAK_REGION_FRAGMENTS}"
  echo "min_atac_fragments=${MIN_ATAC_FRAGMENTS}"
  echo "min_peak_fraction=${MIN_PEAK_FRACTION}"
  echo "combine_mode=${COMBINE_MODE}"
  echo "gex_evidence=${GEX_EVIDENCE_OUT}"
  echo "atac_evidence=${ATAC_EVIDENCE_OUT}"
  echo "multiome_calls=${MULTIOME_CALLS_OUT}"
  echo "comparison_vs_arc=${COMPARE_OUT}"
  echo "comparison_summary=${SUMMARY_OUT}"
  echo "run_compare=${RUN_COMPARE}"
} > "${MANIFEST_OUT}"

echo "Multiome cell-calling run complete." >&2
echo "arc_table=${ARC_TABLE_OUT}" >&2
echo "gex_evidence=${GEX_EVIDENCE_OUT}" >&2
echo "atac_evidence=${ATAC_EVIDENCE_OUT}" >&2
echo "multiome_calls=${MULTIOME_CALLS_OUT}" >&2
if [[ "${RUN_COMPARE}" -eq 1 ]]; then
  echo "comparison_vs_arc=${COMPARE_OUT}" >&2
  echo "comparison_summary=${SUMMARY_OUT}" >&2
fi
echo "run_manifest=${MANIFEST_OUT}" >&2
