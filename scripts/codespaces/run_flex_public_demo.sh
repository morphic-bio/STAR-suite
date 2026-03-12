#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
# shellcheck source=public_demo_sources.sh
source "${SCRIPT_DIR}/public_demo_sources.sh"

OUTDIR="${OUTDIR:-${RUN_ROOT}/flex_public_demo}"
THREADS="${THREADS:-4}"
READ_LIMIT="${READ_LIMIT:-4000}"
SOURCE_READ_LIMIT="${SOURCE_READ_LIMIT:-200000}"
SOLO_CB_START="${SOLO_CB_START:-${CODESPACES_PUBLIC_10X_MALE_SOLO_CB_START:-1}}"
SOLO_CB_LEN="${SOLO_CB_LEN:-${CODESPACES_PUBLIC_10X_MALE_SOLO_CB_LEN:-16}}"
SOLO_UMI_START="${SOLO_UMI_START:-${CODESPACES_PUBLIC_10X_MALE_SOLO_UMI_START:-17}}"
SOLO_UMI_LEN="${SOLO_UMI_LEN:-${CODESPACES_PUBLIC_10X_MALE_SOLO_UMI_LEN:-10}}"
DRY_RUN=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Build a small public-data Flex demo from a public male 10x GEX run. The script
builds a public chr22+chrY mini-reference, derives a GEX fixture from public
raw FASTQs, then overlays a small chr22+chrY probe/sample-tag demo surface.

Options:
  --outdir DIR            Output directory. Default: ${RUN_ROOT}/flex_public_demo
  --threads N             STAR threads. Default: ${THREADS}
  --read-limit N          Synthetic Flex reads to generate. Default: ${READ_LIMIT}
  --source-read-limit N   Public GEX read pairs to inspect before chr22+chrY
                          filtering. Default: ${SOURCE_READ_LIMIT}
  --solo-cb-start N       Cell barcode start for the public source. Default: ${SOLO_CB_START}
  --solo-cb-len N         Cell barcode length for the public source. Default: ${SOLO_CB_LEN}
  --solo-umi-start N      UMI start for the public source. Default: ${SOLO_UMI_START}
  --solo-umi-len N        UMI length for the public source. Default: ${SOLO_UMI_LEN}
  --dry-run               Prepare assets and write RUN_COMMAND.sh without executing
  -h, --help              Show help
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --read-limit) READ_LIMIT="$2"; shift 2 ;;
    --source-read-limit) SOURCE_READ_LIMIT="$2"; shift 2 ;;
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
"${SCRIPT_DIR}/fetch_public_chr22y_reference.sh"
DERIVE_ARGS=(--threads "${THREADS}" --read-limit "${SOURCE_READ_LIMIT}")
if [[ "${DRY_RUN}" == "1" ]]; then
  DERIVE_ARGS+=(--dry-run)
fi
"${SCRIPT_DIR}/derive_public_chr22y_gex_fixture.sh" "${DERIVE_ARGS[@]}"

FIXTURE_GEX_DIR="${DATA_ROOT}/public_chr22y_gex_fixture/gex/public_chr22y_demo"
REF_ROOT="${DATA_ROOT}/public_human_chr22y_ref"
OUTDIR="$(mkdir -p "${OUTDIR}" && cd "${OUTDIR}" && pwd)"
ASSET_DIR="${OUTDIR}/assets"
python3 "${SCRIPT_DIR}/generate_public_chr22y_flex_demo.py" \
  --gex-fastq-dir "${FIXTURE_GEX_DIR}" \
  --ref-root "${REF_ROOT}" \
  --outdir "${ASSET_DIR}" \
  --read-limit "${READ_LIMIT}"

WHITELIST_PATH="${ASSET_DIR}/whitelist_TRU.txt"
[[ -f "${WHITELIST_PATH}" ]] || die "Missing demo whitelist: ${WHITELIST_PATH}"

python3 - <<'PY' "${ASSET_DIR}/probe_set.csv"
from pathlib import Path
import sys
path = Path(sys.argv[1])
text = path.read_text()
if text.endswith("\n"):
    path.write_text(text[:-1])
PY

FLEX_REF_DIR="${OUTDIR}/flex_reference"
"${REPO_ROOT}/flex/scripts/build_filtered_reference.sh" \
  --probe-set "${ASSET_DIR}/probe_set.csv" \
  --base-fasta "${REF_ROOT}/fasta/genome.fa" \
  --base-gtf "${REF_ROOT}/genes/genes.gtf" \
  --work-dir "${FLEX_REF_DIR}" \
  --skip-filter
"${REPO_ROOT}/flex/scripts/make_filtered_star_index.sh" \
  --filtered-reference "${FLEX_REF_DIR}/filtered_reference" \
  --output-dir "${OUTDIR}/star_index" \
  --threads "${THREADS}" \
  --sa-index-n-bases 11 \
  --star-bin "${STAR_BIN}"

CMD=(
  "${REPO_ROOT}/scripts/run_flex_cr_config.sh"
  --cr-config "${ASSET_DIR}/config.csv"
  --genome-dir "${OUTDIR}/star_index"
  --cb-whitelist "${WHITELIST_PATH}"
  --solo-cb-start "${SOLO_CB_START}"
  --solo-cb-len "${SOLO_CB_LEN}"
  --solo-umi-start "${SOLO_UMI_START}"
  --solo-umi-len "${SOLO_UMI_LEN}"
  --out-base "${OUTDIR}"
  --run-id run
  --threads "${THREADS}"
)
write_command_script "${OUTDIR}/RUN_COMMAND.sh" env STAR_BIN="${STAR_BIN}" "${CMD[@]}"

if [[ "${DRY_RUN}" == "1" ]]; then
  log "Prepared ${OUTDIR}/RUN_COMMAND.sh"
  exit 0
fi

env STAR_BIN="${STAR_BIN}" "${CMD[@]}"
