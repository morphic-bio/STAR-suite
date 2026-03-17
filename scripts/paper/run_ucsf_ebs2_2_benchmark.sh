#!/usr/bin/env bash
# Paper benchmark: UCSF EBs2_2 Perturb-seq (CRISPRa v2)
# Compares STAR CR-compat mode against CellRanger 9 reference output.
#
# Dataset:  UCSF EBs2_2 (AALG1 GEX + AALG2 CRISPRa guides)
# Chemistry: NXT namespace (2-column whitelist)
# Libraries: GEX + CRISPR Guide Capture (2 libraries)
#
# Modes:
#   --no-yremove  (default)  Standard GEX+feature run, no Y-chromosome removal
#   --yremove                GEX+feature run with Y-chromosome removal (BAM + Y/noY split)
#
# Both modes produce Solo GeneFull + outs/filtered_feature_bc_matrix + crispr_analysis.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# ── Configurable paths ──────────────────────────────────────────────
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
DATASET_ROOT="${UCSF_DATASET_ROOT:-/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2}"
GEX_DIR="${UCSF_GEX_DIR:-${DATASET_ROOT}/GEX}"
GUIDE_DIR="${UCSF_GUIDE_DIR:-${DATASET_ROOT}/guides}"
FEATURE_REF="${UCSF_FEATURE_REF:-/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv}"
SOLO_CB_WHITELIST="${UCSF_SOLO_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt}"
CR_WHITELIST="${UCSF_CR_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt}"
GENOME_DIR="${UCSF_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
THREADS="${UCSF_THREADS:-32}"
MIN_UMI="${UCSF_MIN_UMI:-3}"
YREMOVE=0
OUTDIR=""

usage() {
  cat <<EOF
Usage: $(basename "$0") [--yremove | --no-yremove] [--outdir DIR] [options]

Options:
  --yremove            Enable Y-chromosome removal (emits Y/noY BAM + FASTQs)
  --no-yremove         Standard run without Y-removal (default)
  --outdir DIR         Output directory (auto-generated if omitted)
  --threads N          Thread count (default: ${THREADS})
  --star-bin PATH      STAR binary path
  -h, --help           Show help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --yremove)       YREMOVE=1; shift ;;
    --no-yremove)    YREMOVE=0; shift ;;
    --outdir)        OUTDIR="$2"; shift 2 ;;
    --threads)       THREADS="$2"; shift 2 ;;
    --star-bin)      STAR_BIN="$2"; shift 2 ;;
    -h|--help)       usage; exit 0 ;;
    *)               echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

MODE_LABEL="standard"
[[ "${YREMOVE}" == "1" ]] && MODE_LABEL="yremove"

if [[ -z "${OUTDIR}" ]]; then
  OUTDIR="/storage/paper_bench_EBs2_2_${MODE_LABEL}_$(date +%Y%m%d_%H%M%S)"
fi

# ── Validate inputs ─────────────────────────────────────────────────
[[ -x "${STAR_BIN}" ]] || { echo "ERROR: STAR binary not found: ${STAR_BIN}" >&2; exit 1; }
for f in "${SOLO_CB_WHITELIST}" "${CR_WHITELIST}" "${FEATURE_REF}"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing file: $f" >&2; exit 1; }
done
for d in "${GENOME_DIR}" "${GEX_DIR}" "${GUIDE_DIR}"; do
  [[ -d "$d" ]] || { echo "ERROR: Missing directory: $d" >&2; exit 1; }
done

mkdir -p "${OUTDIR}"

echo "=========================================="
echo "Paper Benchmark: UCSF EBs2_2 (${MODE_LABEL})"
echo "=========================================="
echo "STAR binary:   ${STAR_BIN}"
echo "Output dir:    ${OUTDIR}"
echo "Threads:       ${THREADS}"
echo "Chemistry:     NXT (2-column whitelist, auto-detect)"
echo "Min UMI:       ${MIN_UMI}"
echo "Y-removal:     ${YREMOVE}"
echo ""

# ── Build FASTQ file lists ──────────────────────────────────────────
GEX_R1=$(ls "${GEX_DIR}/"*R1*.fastq.gz | sort | paste -sd, -)
GEX_R2=$(ls "${GEX_DIR}/"*R2*.fastq.gz | sort | paste -sd, -)
GUIDE_R1=$(ls "${GUIDE_DIR}/"*R1*.fastq.gz | sort | paste -sd, -)
GUIDE_R2=$(ls "${GUIDE_DIR}/"*R2*.fastq.gz | sort | paste -sd, -)

echo "GEX lanes:   $(echo "${GEX_R1}" | tr ',' '\n' | wc -l)"
echo "Guide lanes: $(echo "${GUIDE_R1}" | tr ',' '\n' | wc -l)"
echo ""

# ── Create multi-config ─────────────────────────────────────────────
MULTI_CONFIG="${OUTDIR}/pf_multi_config.csv"
cat > "${MULTI_CONFIG}" <<EOF
[libraries]
fastqs,sample,library_type,feature_types
${GEX_DIR},EBs2_2_AALG1,Gene Expression,Gene Expression
${GUIDE_DIR},EBs2_2_AALG2,CRISPR Guide Capture,CRISPR Guide Capture

[feature]
ref,${FEATURE_REF}
EOF

# ── Build STAR command ───────────────────────────────────────────────
CMD=(
  "${STAR_BIN}"
  --runThreadN "${THREADS}"
  --genomeDir "${GENOME_DIR}"
  --readFilesIn "${GEX_R2},${GUIDE_R2}" "${GEX_R1},${GUIDE_R1}"
  --readFilesCommand zcat
  --outFileNamePrefix "${OUTDIR}/"
  --clipAdapterType CellRanger4
  --clip3pPolyG yes
  --alignEndsType Local
  --chimSegmentMin 1000000
  --soloType CB_UMI_Simple
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
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
  --pfMultiConfig "${MULTI_CONFIG}"
  --crChemistry auto
  --crOutputChemistry TRU
  --crWhitelist "${CR_WHITELIST}"
  --crFeatureRef "${FEATURE_REF}"
  --crMinUmi "${MIN_UMI}"
  --crAssignMaxHamming 1
  --crAssignFeatureOffset -1
  --crAssignLimitSearch -1
  --crAssignMinCounts 0
  --crAssignMaxBarcodeMismatches 5
  --crAssignFeatureN 0
  --crAssignBarcodeN 1
  --crAssignConsumerThreads -1
  --crAssignSearchThreads 1
  --crAssignSkipQcOutputs 1
  --dynamicThreadInterface 1
  --dynamicThreadConstMapPermits "${THREADS}"
  --dynamicThreadTelemetry 1
)

if [[ "${YREMOVE}" == "1" ]]; then
  CMD+=(
    --outSAMtype BAM Unsorted
    --emitNoYBAM yes
    --emitYNoYFastq yes
    --emitYNoYFastqCompression gz
  )
else
  CMD+=(--outSAMtype None)
fi

# ── Save command for reproducibility ─────────────────────────────────
{
  echo '#!/usr/bin/env bash'
  echo 'set -euo pipefail'
  printf '%q ' "${CMD[@]}"
  printf '\n'
} > "${OUTDIR}/RUN_COMMAND.sh"
chmod +x "${OUTDIR}/RUN_COMMAND.sh"

echo "Running STAR (${MODE_LABEL} mode)..."
START_SEC=$SECONDS

"${CMD[@]}"

ELAPSED=$(( SECONDS - START_SEC ))
ELAPSED_MIN=$(awk "BEGIN {printf \"%.1f\", ${ELAPSED}/60}")

echo ""
echo "=========================================="
echo "STAR completed in ${ELAPSED_MIN} min (${ELAPSED} s)"
echo "=========================================="

# ── Validate expected outputs ────────────────────────────────────────
PASS=true
REQUIRED_OUTPUTS=(
  "${OUTDIR}/Solo.out/GeneFull/filtered/barcodes.tsv"
  "${OUTDIR}/Solo.out/GeneFull/filtered/matrix.mtx"
  "${OUTDIR}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
  "${OUTDIR}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
  "${OUTDIR}/outs/crispr_analysis/protospacer_calls_per_cell.csv"
)
if [[ "${YREMOVE}" == "1" ]]; then
  REQUIRED_OUTPUTS+=(
    "${OUTDIR}/Aligned.out_Y.bam"
    "${OUTDIR}/Aligned.out_noY.bam"
  )
fi

for f in "${REQUIRED_OUTPUTS[@]}"; do
  if [[ -f "$f" ]]; then
    echo "  OK: $(basename "$f")"
  else
    echo "  MISSING: $f"
    PASS=false
  fi
done

STAR_CELLS=$(wc -l < "${OUTDIR}/Solo.out/GeneFull/filtered/barcodes.tsv" 2>/dev/null || echo "?")
echo ""
echo "STAR filtered cells: ${STAR_CELLS}"
echo "Wall time: ${ELAPSED_MIN} min"

if [[ "${PASS}" == "true" ]]; then
  echo "TEST PASSED"
else
  echo "TEST FAILED: missing outputs"
  exit 1
fi

# ── Summary file ─────────────────────────────────────────────────────
cat > "${OUTDIR}/BENCHMARK_SUMMARY.txt" <<EOF
dataset=UCSF_EBs2_2
mode=${MODE_LABEL}
chemistry=NXT
star_cells=${STAR_CELLS}
wall_seconds=${ELAPSED}
wall_minutes=${ELAPSED_MIN}
threads=${THREADS}
bam=$(if [[ "${YREMOVE}" == "1" ]]; then echo "unsorted+yremove"; else echo "none"; fi)
dynamic_permits=${THREADS}
bootstrap_tiered_hash=yes
multimap_rescue=yes
outdir=${OUTDIR}
EOF

echo ""
echo "Summary: ${OUTDIR}/BENCHMARK_SUMMARY.txt"
echo "Done."
