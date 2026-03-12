#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
# shellcheck source=public_demo_sources.sh
source "${SCRIPT_DIR}/public_demo_sources.sh"

CR_INPUT_HELPER="${CR_INPUT_HELPER:-${REPO_ROOT}/scripts/ucsf_parity/render_star_inputs_from_cr_config.py}"
OUTDIR="${OUTDIR:-${RUN_ROOT}/perturb_public_demo}"
THREADS="${THREADS:-4}"
READ_LIMIT="${READ_LIMIT:-4000}"
SOURCE_READ_LIMIT="${SOURCE_READ_LIMIT:-200000}"
CR_MIN_UMI="${CR_MIN_UMI:-1}"
CR_CHEMISTRY="${CR_CHEMISTRY:-TRU}"
SOLO_CB_START="${SOLO_CB_START:-${CODESPACES_PUBLIC_10X_MALE_SOLO_CB_START:-1}}"
SOLO_CB_LEN="${SOLO_CB_LEN:-${CODESPACES_PUBLIC_10X_MALE_SOLO_CB_LEN:-16}}"
SOLO_UMI_START="${SOLO_UMI_START:-${CODESPACES_PUBLIC_10X_MALE_SOLO_UMI_START:-17}}"
SOLO_UMI_LEN="${SOLO_UMI_LEN:-${CODESPACES_PUBLIC_10X_MALE_SOLO_UMI_LEN:-10}}"
DRY_RUN=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Build a small public-data perturb demo from a public male 10x GEX run. The
script derives a chr22+chrY GEX fixture from public raw FASTQs, builds a mini
STAR index for chr22+chrY, and then generates a small guide-capture companion
library so the perturb wrapper can run inside Codespaces.

Options:
  --outdir DIR            Output directory. Default: ${RUN_ROOT}/perturb_public_demo
  --threads N             STAR threads. Default: ${THREADS}
  --read-limit N          Synthetic guide reads to generate. Default: ${READ_LIMIT}
  --source-read-limit N   Public GEX read pairs to inspect before chr22+chrY
                          filtering. Default: ${SOURCE_READ_LIMIT}
  --cr-min-umi N          CRISPR UMI threshold for the small demo. Default: ${CR_MIN_UMI}
  --cr-chemistry MODE     Barcode chemistry override for the demo whitelist.
                          Default: ${CR_CHEMISTRY}
  --solo-cb-len N         Cell barcode length for the public source. Default: ${SOLO_CB_LEN}
  --solo-cb-start N       Cell barcode start for the public source. Default: ${SOLO_CB_START}
  --solo-umi-start N      UMI start for the public source. Default: ${SOLO_UMI_START}
  --solo-umi-len N        UMI length for the public source. Default: ${SOLO_UMI_LEN}
  --dry-run               Prepare inputs and write RUN_COMMAND.sh without executing
  -h, --help              Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --read-limit) READ_LIMIT="$2"; shift 2 ;;
    --source-read-limit) SOURCE_READ_LIMIT="$2"; shift 2 ;;
    --cr-min-umi) CR_MIN_UMI="$2"; shift 2 ;;
    --cr-chemistry) CR_CHEMISTRY="$2"; shift 2 ;;
    --solo-cb-start) SOLO_CB_START="$2"; shift 2 ;;
    --solo-cb-len) SOLO_CB_LEN="$2"; shift 2 ;;
    --solo-umi-start) SOLO_UMI_START="$2"; shift 2 ;;
    --solo-umi-len) SOLO_UMI_LEN="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

ensure_star_built
need_cmd python3
"${SCRIPT_DIR}/build_public_chr22y_index.sh" --threads "${THREADS}"
DERIVE_ARGS=(--threads "${THREADS}" --read-limit "${SOURCE_READ_LIMIT}")
if [[ "${DRY_RUN}" == "1" ]]; then
  DERIVE_ARGS+=(--dry-run)
fi
"${SCRIPT_DIR}/derive_public_chr22y_gex_fixture.sh" "${DERIVE_ARGS[@]}"

FIXTURE_GEX_DIR="${DATA_ROOT}/public_chr22y_gex_fixture/gex/public_chr22y_demo"
REF_ROOT="${DATA_ROOT}/public_human_chr22y_ref"
INDEX_DIR="${INDEX_ROOT}/public_human_chr22y_star"
OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
ASSET_DIR="${OUTDIR}/assets"
RUN_DIR="${OUTDIR}/run"
python3 "${SCRIPT_DIR}/generate_public_chr22y_perturb_companion.py" \
  --gex-fastq-dir "${FIXTURE_GEX_DIR}" \
  --ref-root "${REF_ROOT}" \
  --outdir "${ASSET_DIR}" \
  --read-limit "${READ_LIMIT}"

WHITELIST_PATH="${ASSET_DIR}/whitelist_TRU.txt"
[[ -f "${WHITELIST_PATH}" ]] || die "Missing demo whitelist: ${WHITELIST_PATH}"

mkdir -p "${RUN_DIR}"
python3 "${CR_INPUT_HELPER}" \
  --config "${ASSET_DIR}/config.csv" \
  --pf-multi-out "${RUN_DIR}/pf_multi_config.from_cr.csv" \
  --emit-env > "${RUN_DIR}/cr_config_inputs.env"
# shellcheck disable=SC1090
source "${RUN_DIR}/cr_config_inputs.env"

CMD=(
  "${STAR_BIN}"
  --runThreadN "${THREADS}"
  --genomeDir "${INDEX_DIR}"
  --readFilesIn "${GEX_R2},${GUIDE_R2}" "${GEX_R1},${GUIDE_R1}"
  --readFilesCommand zcat
  --outFileNamePrefix "${RUN_DIR}/"
  --outSAMtype BAM SortedByCoordinate
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN gx gn sF sS sQ sM
  --clipAdapterType CellRanger4
  --alignEndsType Local
  --chimSegmentMin 1000000
  --soloType CB_UMI_Simple
  --soloCBstart "${SOLO_CB_START}"
  --soloCBlen "${SOLO_CB_LEN}"
  --soloUMIstart "${SOLO_UMI_START}"
  --soloUMIlen "${SOLO_UMI_LEN}"
  --soloBarcodeReadLength 0
  --soloCBwhitelist "${WHITELIST_PATH}"
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
  --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR
  --soloMultiMappers Unique
  --soloCellFilter EmptyDrops_CR
  --soloCbUbRequireTogether no
  --soloStrand Unstranded
  --soloFeatures GeneFull
  --soloCrGexFeature genefull
  --pfMultiConfig "${PF_MULTI_CONFIG}"
  --crChemistry "${CR_CHEMISTRY}"
  --crFeatureRef "${CR_FEATURE_REF}"
  --crWhitelist "${WHITELIST_PATH}"
  --crMinUmi "${CR_MIN_UMI}"
)

{
  printf 'date_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'star_bin=%s\n' "${STAR_BIN}"
  printf 'ref_root=%s\n' "${REF_ROOT}"
  printf 'index_dir=%s\n' "${INDEX_DIR}"
  printf 'asset_dir=%s\n' "${ASSET_DIR}"
  printf 'run_dir=%s\n' "${RUN_DIR}"
  printf 'cr_config=%s\n' "${ASSET_DIR}/config.csv"
  printf 'pf_multi_config=%s\n' "${PF_MULTI_CONFIG}"
  printf 'gex_r1=%s\n' "${GEX_R1}"
  printf 'gex_r2=%s\n' "${GEX_R2}"
  printf 'guide_r1=%s\n' "${GUIDE_R1}"
  printf 'guide_r2=%s\n' "${GUIDE_R2}"
  printf 'cr_feature_ref=%s\n' "${CR_FEATURE_REF}"
  printf 'cb_whitelist=%s\n' "${WHITELIST_PATH}"
  printf 'cr_chemistry=%s\n' "${CR_CHEMISTRY}"
  printf 'solo_cb_start=%s\n' "${SOLO_CB_START}"
  printf 'solo_cb_len=%s\n' "${SOLO_CB_LEN}"
  printf 'solo_umi_start=%s\n' "${SOLO_UMI_START}"
  printf 'solo_umi_len=%s\n' "${SOLO_UMI_LEN}"
  printf 'cr_min_umi=%s\n' "${CR_MIN_UMI}"
} > "${RUN_DIR}/RUN_MANIFEST.txt"

write_command_script "${OUTDIR}/RUN_COMMAND.sh" "${CMD[@]}"
write_command_script "${RUN_DIR}/RUN_COMMAND.sh" "${CMD[@]}"

if [[ "${DRY_RUN}" == "1" ]]; then
  log "Prepared ${OUTDIR}/RUN_COMMAND.sh"
  exit 0
fi

"${CMD[@]}"
