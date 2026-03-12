#!/usr/bin/env bash
set -euo pipefail

# Canonical UCSF full-sample STAR run (GEX + guides) with CR-compat rescue.
# This script exists to pin parameters and eliminate ambiguity across reruns.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CR_INPUT_HELPER="${CR_INPUT_HELPER:-${SCRIPT_DIR}/ucsf_parity/render_star_inputs_from_cr_config.py}"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
THREADS="${THREADS:-32}"
GENOME_DIR="${GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
CR_CONFIG="${CR_CONFIG:-}"

GEX_R2="${GEX_R2:-/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L001_R2_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L002_R2_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L003_R2_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L004_R2_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L005_R2_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L006_R2_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L007_R2_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L008_R2_001.fastq.gz}"
GEX_R1="${GEX_R1:-/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L001_R1_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L002_R1_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L003_R1_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L004_R1_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L005_R1_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L006_R1_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L007_R1_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/GEX/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L008_R1_001.fastq.gz}"
GUIDE_R2="${GUIDE_R2:-/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/guides/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L001_R2_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/guides/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L002_R2_001.fastq.gz}"
GUIDE_R1="${GUIDE_R1:-/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/guides/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L001_R1_001.fastq.gz,/storage/ucsf-full/bench_20260218_dynamic_first/staged_iPSC2_1_AALG2/guides/iPSC2_1_AALG2/iPSC2_1_AALG2_S30_L002_R1_001.fastq.gz}"

# Pinned to NXT by request for UCSF compatibility-mode runs.
SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
CR_WHITELIST="${CR_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
CR_FEATURE_REF="${CR_FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
PF_MULTI_CONFIG="${PF_MULTI_CONFIG:-/storage/ucsf-full/bench_20260218_dynamic_first/pf_multi_config.csv}"

OUT_BASE="${OUT_BASE:-/storage/ucsf-full/bench_20260218_dynamic_first/runs}"
RUN_ID="${RUN_ID:-star_full_iPSC2_1_AALG2_forward_rescue_guides_$(date +%Y%m%d_%H%M%S)}"
DRY_RUN="${DRY_RUN:-0}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --cr-config FILE   Cell Ranger config.csv to use as the source of truth for UCSF inputs
  --threads N        Override thread count
  --genome-dir DIR   Override STAR genomeDir
  --out-base DIR     Override output base directory
  --run-id ID        Override run identifier
  --dry-run          Write manifest/command only
  -h, --help         Show help
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cr-config)
      CR_CONFIG="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --genome-dir)
      GENOME_DIR="$2"
      shift 2
      ;;
    --out-base)
      OUT_BASE="$2"
      shift 2
      ;;
    --run-id)
      RUN_ID="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "Unknown option: $1"
      ;;
  esac
done

OUT_DIR="${OUT_BASE}/${RUN_ID}"
mkdir -p "${OUT_DIR}"

if [[ -n "${CR_CONFIG}" ]]; then
  [[ -f "${CR_CONFIG}" ]] || die "Missing --cr-config file: ${CR_CONFIG}"
  [[ -x "${CR_INPUT_HELPER}" ]] || die "Missing helper: ${CR_INPUT_HELPER}"
  CR_ENV_FILE="${OUT_DIR}/cr_config_inputs.env"
  GENERATED_PF_MULTI="${OUT_DIR}/pf_multi_config.from_cr.csv"
  python3 "${CR_INPUT_HELPER}" \
    --config "${CR_CONFIG}" \
    --pf-multi-out "${GENERATED_PF_MULTI}" \
    --emit-env > "${CR_ENV_FILE}"
  # shellcheck disable=SC1090
  source "${CR_ENV_FILE}"
  PF_MULTI_CONFIG="${GENERATED_PF_MULTI}"
fi

READS_R2="${GEX_R2},${GUIDE_R2}"
READS_R1="${GEX_R1},${GUIDE_R1}"

CMD=(
  "${STAR_BIN}"
  --runThreadN "${THREADS}"
  --genomeDir "${GENOME_DIR}"
  --readFilesIn "${READS_R2}" "${READS_R1}"
  --readFilesCommand zcat
  --outFileNamePrefix "${OUT_DIR}/"
  --outSAMtype None
  --clipAdapterType CellRanger4
  --alignEndsType Local
  --chimSegmentMin 1000000
  --soloType CB_UMI_Simple
  --soloCBstart 1
  --soloCBlen 16
  --soloUMIstart 17
  --soloUMIlen 12
  --soloBarcodeReadLength 0
  --soloCBwhitelist "${SOLO_CB_WHITELIST}"
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
  --soloUMIfiltering MultiGeneUMI_CR
  --soloUMIdedup 1MM_CR
  --soloMultiMappers Unique
  --soloCellFilter EmptyDrops_CR
  --soloCbUbRequireTogether no
  --soloStrand Forward
  --soloFeatures GeneFull
  --soloCrGexFeature genefull
  --soloCrMultimapRescue yes
  --pfMultiConfig "${PF_MULTI_CONFIG}"
  --crFeatureRef "${CR_FEATURE_REF}"
  --crWhitelist "${CR_WHITELIST}"
  --crMinUmi 3
  --crAssignMaxHamming 1
  --crAssignFeatureOffset 0
  --crAssignLimitSearch -1
  --crAssignMinCounts 0
  --crAssignMaxBarcodeMismatches 5
  --crAssignFeatureN 0
  --crAssignBarcodeN 1
  --crAssignConsumerThreads 4
  --crAssignSearchThreads 4
)

printf "run_id=%s\n" "${RUN_ID}" > "${OUT_DIR}/RUN_MANIFEST.txt"
printf "date_utc=%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "${OUT_DIR}/RUN_MANIFEST.txt"
printf "star_bin=%s\n" "${STAR_BIN}" >> "${OUT_DIR}/RUN_MANIFEST.txt"
printf "threads=%s\n" "${THREADS}" >> "${OUT_DIR}/RUN_MANIFEST.txt"
printf "out_dir=%s\n" "${OUT_DIR}" >> "${OUT_DIR}/RUN_MANIFEST.txt"
if [[ -n "${CR_CONFIG}" ]]; then
  printf "cr_config=%s\n" "${CR_CONFIG}" >> "${OUT_DIR}/RUN_MANIFEST.txt"
  printf "cr_gene_expression_reference=%s\n" "${CR_GENE_EXPRESSION_REFERENCE:-}" >> "${OUT_DIR}/RUN_MANIFEST.txt"
  printf "cr_gene_expression_chemistry=%s\n" "${CR_GENE_EXPRESSION_CHEMISTRY:-}" >> "${OUT_DIR}/RUN_MANIFEST.txt"
fi
{
  echo '#!/usr/bin/env bash'
  echo 'set -euo pipefail'
  printf "%q " "${CMD[@]}"
  printf "\n"
} > "${OUT_DIR}/RUN_COMMAND.sh"
chmod +x "${OUT_DIR}/RUN_COMMAND.sh"

echo "=== UCSF full run parameters (guide-optimized) ==="
cat "${OUT_DIR}/RUN_MANIFEST.txt"
echo "=== command ==="
cat "${OUT_DIR}/RUN_COMMAND.sh"

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "DRY_RUN=1; not executing STAR."
  exit 0
fi

"${CMD[@]}"
