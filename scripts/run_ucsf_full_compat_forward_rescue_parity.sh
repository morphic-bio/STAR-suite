#!/usr/bin/env bash
set -euo pipefail

# Canonical parity wrapper for UCSF full run (STAR compat forward+rescue vs CR).
# Pins comparison knobs so reruns are reproducible and auditable.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

STAR_RUN="${STAR_RUN:-}"
CR_RUN="${CR_RUN:-/storage/ucsf-full/bench_20260218_dynamic_first/cellranger_runs/cr_full_iPSC2_1_AALG2_crstar32_20260218_205804}"
OUT_DIR="${OUT_DIR:-/tmp/ucsf_full_compat_parity_$(date +%Y%m%d_%H%M%S)}"

# Pinned parity settings.
NORMALIZE_BARCODES="${NORMALIZE_BARCODES:-1}"
BARCODE_TRANSLATION="${BARCODE_TRANSLATION:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt.gz}"
TRANSLATION_DIRECTION="${TRANSLATION_DIRECTION:-left-to-right}"
TRANSLATE_SIDE="${TRANSLATE_SIDE:-both}"
GENE_CORR_MIN_COUNTS="${GENE_CORR_MIN_COUNTS:-20}"
GENE_CORR_MIN_CELLS_PCT="${GENE_CORR_MIN_CELLS_PCT:-0.01}"

if [[ -z "${STAR_RUN}" ]]; then
  cat >&2 <<'EOF'
ERROR: STAR_RUN is required.
Example:
  STAR_RUN=/storage/ucsf-full/bench_20260218_dynamic_first/runs/<your_star_run> \
  scripts/run_ucsf_full_compat_forward_rescue_parity.sh
EOF
  exit 2
fi

mkdir -p "${OUT_DIR}"
MANIFEST="${OUT_DIR}/PARITY_MANIFEST.txt"

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "star_run=${STAR_RUN}"
  echo "cr_run=${CR_RUN}"
  echo "normalize_barcodes=${NORMALIZE_BARCODES}"
  echo "barcode_translation=${BARCODE_TRANSLATION}"
  echo "translation_direction=${TRANSLATION_DIRECTION}"
  echo "translate_side=${TRANSLATE_SIDE}"
  echo "gene_corr_min_counts=${GENE_CORR_MIN_COUNTS}"
  echo "gene_corr_min_cells_pct=${GENE_CORR_MIN_CELLS_PCT}"
} > "${MANIFEST}"

if [[ "${NORMALIZE_BARCODES}" == "1" ]]; then
  CMD=(
    env
    STAR_RUN="${STAR_RUN}"
    CR_RUN="${CR_RUN}"
    OUT_DIR="${OUT_DIR}"
    BARCODE_TRANSLATION="${BARCODE_TRANSLATION}"
    TRANSLATION_DIRECTION="${TRANSLATION_DIRECTION}"
    TRANSLATE_SIDE="${TRANSLATE_SIDE}"
    GENE_CORR_MIN_COUNTS="${GENE_CORR_MIN_COUNTS}"
    GENE_CORR_MIN_CELLS_PCT="${GENE_CORR_MIN_CELLS_PCT}"
    "${REPO_ROOT}/scripts/run_gex_feature_parity_checks.sh"
  )
else
  CMD=(
    env
    STAR_RUN="${STAR_RUN}"
    CR_RUN="${CR_RUN}"
    OUT_DIR="${OUT_DIR}"
    GENE_CORR_MIN_COUNTS="${GENE_CORR_MIN_COUNTS}"
    GENE_CORR_MIN_CELLS_PCT="${GENE_CORR_MIN_CELLS_PCT}"
    "${REPO_ROOT}/scripts/run_gex_feature_parity_checks.sh"
    --no-barcode-normalization
  )
fi

{
  echo '#!/usr/bin/env bash'
  echo 'set -euo pipefail'
  printf "%q " "${CMD[@]}"
  printf "\n"
} > "${OUT_DIR}/PARITY_COMMAND.sh"
chmod +x "${OUT_DIR}/PARITY_COMMAND.sh"

echo "=== parity manifest ==="
cat "${MANIFEST}"
echo "=== command ==="
cat "${OUT_DIR}/PARITY_COMMAND.sh"

"${CMD[@]}"
