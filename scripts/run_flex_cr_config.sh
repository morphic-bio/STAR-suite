#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
FLEX_INPUT_HELPER="${FLEX_INPUT_HELPER:-${SCRIPT_DIR}/flex_compat/render_flex_inputs_from_cr_config.py}"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
THREADS="${THREADS:-8}"
GENOME_DIR="${GENOME_DIR:-/storage/flex_filtered_reference/star_index}"
SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST:-/storage/scRNAseq_output/whitelists/737K-fixed-rna-profiling.txt}"
SOLO_CB_START="${SOLO_CB_START:-1}"
SOLO_CB_LEN="${SOLO_CB_LEN:-16}"
SOLO_UMI_START="${SOLO_UMI_START:-17}"
SOLO_UMI_LEN="${SOLO_UMI_LEN:-12}"
SAMPLE_PROBE_CATALOG="${SAMPLE_PROBE_CATALOG:-/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt}"
SAMPLE_PROBE_OFFSET="${SAMPLE_PROBE_OFFSET:-68}"
CR_CONFIG="${CR_CONFIG:-}"
OUT_BASE="${OUT_BASE:-/tmp/flex_cr_config_runs}"
RUN_ID="${RUN_ID:-flex_cr_config_$(date +%Y%m%d_%H%M%S)}"
DRY_RUN="${DRY_RUN:-0}"
SAMPLE_PROBE_CATALOG_EXPLICIT=0
SAMPLE_PROBE_OFFSET_EXPLICIT=0

usage() {
  cat <<EOF
Usage: $(basename "$0") --cr-config FILE [options]

Options:
  --cr-config FILE           Cell Ranger Flex config.csv
  --threads N                STAR threads (default: ${THREADS})
  --genome-dir DIR           STAR genomeDir (default: ${GENOME_DIR})
  --cb-whitelist FILE        Fixed RNA CB whitelist (default: ${SOLO_CB_WHITELIST})
  --solo-cb-start N          Barcode start on R1 (default: ${SOLO_CB_START})
  --solo-cb-len N            Barcode length on R1 (default: ${SOLO_CB_LEN})
  --solo-umi-start N         UMI start on R1 (default: ${SOLO_UMI_START})
  --solo-umi-len N           UMI length on R1 (default: ${SOLO_UMI_LEN})
  --sample-probe-catalog P   Probe barcode catalog (default: ${SAMPLE_PROBE_CATALOG})
  --sample-probe-offset N    Probe offset relative to cDNA read (default: ${SAMPLE_PROBE_OFFSET})
  --out-base DIR             Output base directory (default: ${OUT_BASE})
  --run-id ID                Run directory name (default: ${RUN_ID})
  --dry-run                  Write manifest/command/helpers only
  -h, --help                 Show help
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cr-config) CR_CONFIG="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --genome-dir) GENOME_DIR="$2"; shift 2 ;;
    --cb-whitelist) SOLO_CB_WHITELIST="$2"; shift 2 ;;
    --solo-cb-start) SOLO_CB_START="$2"; shift 2 ;;
    --solo-cb-len) SOLO_CB_LEN="$2"; shift 2 ;;
    --solo-umi-start) SOLO_UMI_START="$2"; shift 2 ;;
    --solo-umi-len) SOLO_UMI_LEN="$2"; shift 2 ;;
    --sample-probe-catalog) SAMPLE_PROBE_CATALOG="$2"; SAMPLE_PROBE_CATALOG_EXPLICIT=1; shift 2 ;;
    --sample-probe-offset) SAMPLE_PROBE_OFFSET="$2"; SAMPLE_PROBE_OFFSET_EXPLICIT=1; shift 2 ;;
    --out-base) OUT_BASE="$2"; shift 2 ;;
    --run-id) RUN_ID="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1" ;;
  esac
done

[[ -n "${CR_CONFIG}" ]] || die "Specify --cr-config FILE"
[[ -x "${STAR_BIN}" ]] || die "Missing STAR binary: ${STAR_BIN}"
[[ -f "${CR_CONFIG}" ]] || die "Missing CR config: ${CR_CONFIG}"
[[ -x "${FLEX_INPUT_HELPER}" ]] || die "Missing helper: ${FLEX_INPUT_HELPER}"
[[ -d "${GENOME_DIR}" ]] || die "Missing genomeDir: ${GENOME_DIR}"
[[ -f "${SOLO_CB_WHITELIST}" ]] || die "Missing whitelist: ${SOLO_CB_WHITELIST}"

OUT_DIR="${OUT_BASE}/${RUN_ID}"
mkdir -p "${OUT_DIR}"

FLEX_ENV_FILE="${OUT_DIR}/cr_config_inputs.env"
SAMPLE_WHITELIST_OUT="${OUT_DIR}/sample_whitelist.from_cr.tsv"
SAMPLE_PROBES_OUT="${OUT_DIR}/sample_probes.from_cr.tsv"
PROBE_LIST_OUT="${OUT_DIR}/probe_list.from_cr.txt"

helper_args=(
  python3 "${FLEX_INPUT_HELPER}"
  --config "${CR_CONFIG}"
  --sample-whitelist-out "${SAMPLE_WHITELIST_OUT}"
  --sample-probes-out "${SAMPLE_PROBES_OUT}"
  --probe-list-out "${PROBE_LIST_OUT}"
  --probe-catalog "${SAMPLE_PROBE_CATALOG}"
)
helper_args+=(--emit-env)
"${helper_args[@]}" > "${FLEX_ENV_FILE}"
# shellcheck disable=SC1090
source "${FLEX_ENV_FILE}"

if [[ -n "${FLEX_SAMPLE_PROBE_CATALOG:-}" ]]; then
  SAMPLE_PROBE_CATALOG="${FLEX_SAMPLE_PROBE_CATALOG}"
fi
if [[ "${SAMPLE_PROBE_OFFSET_EXPLICIT}" != "1" && -n "${FLEX_SAMPLE_PROBE_OFFSET:-}" ]]; then
  SAMPLE_PROBE_OFFSET="${FLEX_SAMPLE_PROBE_OFFSET}"
fi

[[ -f "${SAMPLE_PROBE_CATALOG}" ]] || die "Missing sample probe catalog: ${SAMPLE_PROBE_CATALOG}"

CMD=(
  "${STAR_BIN}"
  --runThreadN "${THREADS}"
  --genomeDir "${GENOME_DIR}"
  --readFilesIn "${GEX_R2}" "${GEX_R1}"
  --readFilesCommand zcat
  --outFileNamePrefix "${OUT_DIR}/"
  --outSAMtype BAM Unsorted
  --soloType CB_UMI_Simple
  --soloCBstart "${SOLO_CB_START}"
  --soloCBlen "${SOLO_CB_LEN}"
  --soloUMIstart "${SOLO_UMI_START}"
  --soloUMIlen "${SOLO_UMI_LEN}"
  --soloBarcodeReadLength 0
  --soloCBwhitelist "${SOLO_CB_WHITELIST}"
  --flex yes
  --soloFlexExpectedCellsPerTag 3000
  --soloSampleWhitelist "${FLEX_SAMPLE_WHITELIST}"
  --soloProbeList "${FLEX_PROBE_LIST}"
  --soloSampleProbes "${FLEX_SAMPLE_PROBES}"
  --soloSampleProbeOffset "${SAMPLE_PROBE_OFFSET}"
  --soloFlexAllowedTags "${FLEX_ALLOWED_TAGS}"
  --soloFlexOutputPrefix "${OUT_DIR}/per_sample"
  --limitIObufferSize 50000000 50000000
  --outSJtype None
  --outBAMcompression 6
  --soloMultiMappers Rescue
  --alignIntronMax 500000
  --outFilterMismatchNmax 6
  --outFilterMismatchNoverReadLmax 1.0
  --outFilterMatchNmin 25
  --outSAMunmapped None
  --outFilterMatchNminOverLread 0
  --outFilterMultimapNmax 10000
  --outFilterMultimapScoreRange 4
  --outSAMmultNmax 10000
  --winAnchorMultimapNmax 200
  --outSAMprimaryFlag AllBestScore
  --outFilterScoreMin 0
  --outFilterScoreMinOverLread 0
  --outSAMattributes NH HI AS nM NM GX GN
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
  --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR
  --soloCellFilter None
  --clipAdapterType CellRanger4
  --soloFeatures Gene
  --alignEndsType Local
  --soloStrand Unstranded
  --chimSegmentMin 1000000
  --soloKeysCompat cr
  --soloSampleSearchNearby no
)

{
  printf 'run_id=%s\n' "${RUN_ID}"
  printf 'date_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'star_bin=%s\n' "${STAR_BIN}"
  printf 'threads=%s\n' "${THREADS}"
  printf 'out_dir=%s\n' "${OUT_DIR}"
  printf 'cr_config=%s\n' "${CR_CONFIG}"
  printf 'cr_gene_expression_reference=%s\n' "${CR_GENE_EXPRESSION_REFERENCE:-}"
  printf 'cr_gene_expression_probe_set=%s\n' "${CR_GENE_EXPRESSION_PROBE_SET:-}"
  printf 'gex_fastq_dirs=%s\n' "${GEX_FASTQ_DIRS:-}"
  printf 'solo_cb_start=%s\n' "${SOLO_CB_START}"
  printf 'solo_cb_len=%s\n' "${SOLO_CB_LEN}"
  printf 'solo_umi_start=%s\n' "${SOLO_UMI_START}"
  printf 'solo_umi_len=%s\n' "${SOLO_UMI_LEN}"
  printf 'sample_whitelist=%s\n' "${FLEX_SAMPLE_WHITELIST}"
  printf 'sample_probes=%s\n' "${FLEX_SAMPLE_PROBES}"
  printf 'probe_list=%s\n' "${FLEX_PROBE_LIST}"
  printf 'sample_probe_catalog=%s\n' "${SAMPLE_PROBE_CATALOG}"
  printf 'sample_probe_offset=%s\n' "${SAMPLE_PROBE_OFFSET}"
  printf 'sample_ids=%s\n' "${FLEX_SAMPLE_IDS:-}"
} > "${OUT_DIR}/RUN_MANIFEST.txt"

{
  echo '#!/usr/bin/env bash'
  echo 'set -euo pipefail'
  printf '%q ' "${CMD[@]}"
  printf '\n'
} > "${OUT_DIR}/RUN_COMMAND.sh"
chmod +x "${OUT_DIR}/RUN_COMMAND.sh"

echo "=== FLEX CR-config run parameters ==="
cat "${OUT_DIR}/RUN_MANIFEST.txt"
echo "=== command ==="
cat "${OUT_DIR}/RUN_COMMAND.sh"

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=1; not executing STAR."
  exit 0
fi

"${CMD[@]}"
