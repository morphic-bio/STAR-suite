#!/usr/bin/env bash
set -euo pipefail

# Canonical UCSF full-sample STAR run (GEX + guides) with CR-compat rescue.
# This script exists to pin parameters and eliminate ambiguity across reruns.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
THREADS="${THREADS:-32}"
GENOME_DIR="${GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"

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
OUT_DIR="${OUT_BASE}/${RUN_ID}"
DRY_RUN="${DRY_RUN:-0}"

mkdir -p "${OUT_DIR}"

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
  --crAssignFeatureN 1
  --crAssignBarcodeN 2
  --crAssignConsumerThreads 4
  --crAssignSearchThreads 4
)

printf "run_id=%s\n" "${RUN_ID}" > "${OUT_DIR}/RUN_MANIFEST.txt"
printf "date_utc=%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "${OUT_DIR}/RUN_MANIFEST.txt"
printf "star_bin=%s\n" "${STAR_BIN}" >> "${OUT_DIR}/RUN_MANIFEST.txt"
printf "threads=%s\n" "${THREADS}" >> "${OUT_DIR}/RUN_MANIFEST.txt"
printf "out_dir=%s\n" "${OUT_DIR}" >> "${OUT_DIR}/RUN_MANIFEST.txt"
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
