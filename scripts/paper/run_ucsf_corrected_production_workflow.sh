#!/usr/bin/env bash
# Paper-owned UCSF production workflow:
# corrected symlinked inputs -> STAR Y-removal + GeneFull/Velocyto -> CellBender/h5ad -> optional Globus.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BATCH_RUNNER="${REPO_ROOT}/scripts/run_ucsf_perturb_yremove_batch.sh"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR.release}"
DATASET_ROOT="${UCSF_DATASET_ROOT:-/mnt/pikachu/ucsf-perturb-seq-corrected}"
FEATURE_REF="${UCSF_FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
GENOME_DIR="${UCSF_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
SOLO_CB_WHITELIST="${UCSF_SOLO_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt}"
CR_WHITELIST="${UCSF_CR_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
THREADS="${UCSF_THREADS:-24}"
CELLBENDER_CPU_CORES="${UCSF_CELLBENDER_CPU_CORES:-${THREADS}}"
OUT_ROOT="${UCSF_OUT_ROOT:-/mnt/pikachu/ucsf-perturb-yremove_velocyto_cellbender_$(date +%Y%m%d_%H%M%S)}"

SAMPLES_CSV=""
ALL_SAMPLES=0
DOWNSAMPLE_READS=0
DOWNSAMPLE_SEED="${DOWNSAMPLE_SEED:-1}"
GLOBUS_SRC_ENDPOINT="${GLOBUS_SRC_ENDPOINT:-}"
GLOBUS_DST_ENDPOINT="${GLOBUS_DST_ENDPOINT:-}"
GLOBUS_DST_ROOT="${GLOBUS_DST_ROOT:-}"
GLOBUS_POLL_SECONDS="${GLOBUS_POLL_SECONDS:-30}"
CELLBENDER_GPU=0
TRIM_QC=0
DRY_RUN=0
STAR_ONLY=0

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Defaults:
  - release binary pinned to ${STAR_BIN}
  - corrected UCSF input root ${DATASET_ROOT}
  - GeneFull + Velocyto (no Gene), Y-removal enabled by the batch runner
  - downstream h5ad generation with CellBender
  - 24 threads / 24 CellBender CPU cores

Options:
  --samples CSV            Comma-separated sample IDs
  --all-samples            Run every corrected UCSF sample
  --downsample-reads N     Downsample each library before STAR
  --downsample-seed N      Downsampling seed (default: ${DOWNSAMPLE_SEED})
  --threads N              STAR threads (default: ${THREADS})
  --star-bin PATH          STAR binary (default: ${STAR_BIN})
  --cellbender-cpu-cores N CellBender CPU cores (default: ${CELLBENDER_CPU_CORES})
  --cellbender-gpu         Enable GPU CellBender
  --out-root DIR           Output root (default: ${OUT_ROOT})
  --dataset-root DIR       Corrected UCSF dataset root
  --feature-ref FILE       Feature reference CSV
  --genome-dir DIR         STAR genomeDir
  --solo-cb-whitelist P    Solo CB whitelist
  --cr-whitelist P         CR whitelist
  --globus-src-endpoint ID Source Globus collection ID
  --globus-dst-endpoint ID Destination Globus collection ID
  --globus-dst-root PATH   Destination root path
  --globus-poll-seconds N  Globus cleanup poll interval
  --trim-qc                Emit read-level trim-QC reports
  --star-only              Run STAR/Globus only; do not run downstream/CellBender
  --dry-run                Prepare manifests only
  -h, --help               Show help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES_CSV="$2"; shift 2 ;;
    --all-samples) ALL_SAMPLES=1; shift ;;
    --downsample-reads) DOWNSAMPLE_READS="$2"; shift 2 ;;
    --downsample-seed) DOWNSAMPLE_SEED="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --star-bin) STAR_BIN="$2"; shift 2 ;;
    --cellbender-cpu-cores) CELLBENDER_CPU_CORES="$2"; shift 2 ;;
    --cellbender-gpu) CELLBENDER_GPU=1; shift ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --dataset-root) DATASET_ROOT="$2"; shift 2 ;;
    --feature-ref) FEATURE_REF="$2"; shift 2 ;;
    --genome-dir) GENOME_DIR="$2"; shift 2 ;;
    --solo-cb-whitelist) SOLO_CB_WHITELIST="$2"; shift 2 ;;
    --cr-whitelist) CR_WHITELIST="$2"; shift 2 ;;
    --globus-src-endpoint) GLOBUS_SRC_ENDPOINT="$2"; shift 2 ;;
    --globus-dst-endpoint) GLOBUS_DST_ENDPOINT="$2"; shift 2 ;;
    --globus-dst-root) GLOBUS_DST_ROOT="$2"; shift 2 ;;
    --globus-poll-seconds) GLOBUS_POLL_SECONDS="$2"; shift 2 ;;
    --trim-qc) TRIM_QC=1; shift ;;
    --star-only) STAR_ONLY=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -x "${BATCH_RUNNER}" ]] || { echo "ERROR: missing batch runner ${BATCH_RUNNER}" >&2; exit 1; }
[[ -x "${STAR_BIN}" ]] || { echo "ERROR: missing release binary ${STAR_BIN}" >&2; exit 1; }

if [[ "${ALL_SAMPLES}" == "0" && -z "${SAMPLES_CSV}" ]]; then
  ALL_SAMPLES=1
fi

mkdir -p "${OUT_ROOT}"

ARGS=(
  "${BATCH_RUNNER}"
  --out-root "${OUT_ROOT}"
  --dataset-root "${DATASET_ROOT}"
  --feature-ref "${FEATURE_REF}"
  --genome-dir "${GENOME_DIR}"
  --solo-cb-whitelist "${SOLO_CB_WHITELIST}"
  --cr-whitelist "${CR_WHITELIST}"
  --star-bin "${STAR_BIN}"
  --threads "${THREADS}"
  --globus-poll-seconds "${GLOBUS_POLL_SECONDS}"
)

if [[ "${STAR_ONLY}" != "1" ]]; then
  ARGS+=(
    --run-downstream
    --adaptive-filter
    --run-cellbender
    --cellbender-cpu-cores "${CELLBENDER_CPU_CORES}"
  )
fi

if [[ "${ALL_SAMPLES}" == "1" ]]; then
  ARGS+=(--all-samples)
else
  ARGS+=(--samples "${SAMPLES_CSV}")
fi
if (( DOWNSAMPLE_READS > 0 )); then
  ARGS+=(--downsample-reads "${DOWNSAMPLE_READS}")
fi
if [[ -n "${GLOBUS_SRC_ENDPOINT}" ]]; then
  ARGS+=(--globus-src-endpoint "${GLOBUS_SRC_ENDPOINT}")
fi
if [[ -n "${GLOBUS_DST_ENDPOINT}" ]]; then
  ARGS+=(--globus-dst-endpoint "${GLOBUS_DST_ENDPOINT}")
fi
if [[ -n "${GLOBUS_DST_ROOT}" ]]; then
  ARGS+=(--globus-dst-root "${GLOBUS_DST_ROOT}")
fi
if [[ "${CELLBENDER_GPU}" == "1" && "${STAR_ONLY}" != "1" ]]; then
  ARGS+=(--cellbender-gpu)
fi
if [[ "${TRIM_QC}" == "1" ]]; then
  ARGS+=(--trim-qc)
fi
if [[ "${DRY_RUN}" == "1" ]]; then
  ARGS+=(--dry-run)
fi

{
  echo '#!/usr/bin/env bash'
  echo 'set -euo pipefail'
  printf 'export DOWNSAMPLE_SEED=%q\n' "${DOWNSAMPLE_SEED}"
  printf '%q ' "${ARGS[@]}"
  printf '\n'
} > "${OUT_ROOT}/RUN_WRAPPER_COMMAND.sh"
chmod +x "${OUT_ROOT}/RUN_WRAPPER_COMMAND.sh"

export DOWNSAMPLE_SEED
"${ARGS[@]}"
