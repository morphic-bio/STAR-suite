#!/usr/bin/env bash
# Paper benchmark: MSK 30polyKO 3-library Perturb-seq
# Compares STAR CR-compat mode against CellRanger 9 reference output.
#
# Dataset:  MSK scRNAseq 30polyKO (S5/S6/S6_1)
# Chemistry: Mixed — GEX (TRU), gRNA/PolyIII (NXT), LARRY (TRU)
# Libraries: GEX + CRISPR Guide Capture (PolyIII) + Custom (LARRY) — 3 libraries
#
# CellRanger requires two separate runs (GEX+gRNA and GEX+LARRY);
# STAR handles all three libraries in a single pass.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# ── Configurable paths ──────────────────────────────────────────────
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
MSK_ROOT="${MSK_ROOT:-/storage/MSK-perturb-comparison}"
FASTQ_ROOT="${MSK_FASTQ_ROOT:-${MSK_ROOT}/msk30ko_full_3lib_20260304_095911/fastqs}"
GEX_DIR="${MSK_GEX_DIR:-${FASTQ_ROOT}/mRNA}"
GRNA_DIR="${MSK_GRNA_DIR:-${FASTQ_ROOT}/PolyIII}"
LARRY_DIR="${MSK_LARRY_DIR:-${FASTQ_ROOT}/LARRY}"
GRNA_FEATURE_REF="${MSK_GRNA_REF:-/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv}"
LARRY_FEATURE_REF="${MSK_LARRY_REF:-/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv}"
SOLO_CB_WHITELIST="${MSK_SOLO_WHITELIST:-/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt}"
GRNA_WHITELIST="${MSK_GRNA_WHITELIST:-/storage/scRNAseq_output/whitelists/3M-february-2018_NXT.txt}"
LARRY_WHITELIST="${MSK_LARRY_WHITELIST:-/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt}"
GENOME_DIR="${MSK_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
OUTDIR="${MSK_OUTDIR:-/storage/MSK-perturb-comparison/paper_bench_$(date +%Y%m%d_%H%M%S)}"
THREADS="${MSK_THREADS:-32}"
MIN_UMI="${MSK_MIN_UMI:-2}"
WITH_BAM="${MSK_WITH_BAM:-1}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --outdir DIR         Output directory (auto-generated if omitted)
  --threads N          Thread count (default: ${THREADS})
  --no-bam             Skip BAM output (faster, but no CB/UB tags)
  --with-bam           Emit unsorted BAM (default)
  --star-bin PATH      STAR binary path
  -h, --help           Show help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)    OUTDIR="$2"; shift 2 ;;
    --threads)   THREADS="$2"; shift 2 ;;
    --no-bam)    WITH_BAM=0; shift ;;
    --with-bam)  WITH_BAM=1; shift ;;
    --star-bin)  STAR_BIN="$2"; shift 2 ;;
    -h|--help)   usage; exit 0 ;;
    *)           echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

# ── Validate inputs ─────────────────────────────────────────────────
[[ -x "${STAR_BIN}" ]] || { echo "ERROR: STAR binary not found: ${STAR_BIN}" >&2; exit 1; }
for f in "${SOLO_CB_WHITELIST}" "${GRNA_WHITELIST}" "${LARRY_WHITELIST}" "${GRNA_FEATURE_REF}" "${LARRY_FEATURE_REF}"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing file: $f" >&2; exit 1; }
done
for d in "${GENOME_DIR}" "${GEX_DIR}" "${GRNA_DIR}" "${LARRY_DIR}"; do
  [[ -d "$d" ]] || { echo "ERROR: Missing directory: $d" >&2; exit 1; }
done

mkdir -p "${OUTDIR}"

echo "=========================================="
echo "Paper Benchmark: MSK 30polyKO (3-library)"
echo "=========================================="
echo "STAR binary:   ${STAR_BIN}"
echo "Output dir:    ${OUTDIR}"
echo "Threads:       ${THREADS}"
echo "Chemistry:     Mixed (GEX=TRU, gRNA=NXT, LARRY=TRU)"
echo "Min UMI:       ${MIN_UMI}"
echo "BAM output:    $(if [[ "${WITH_BAM}" == "1" ]]; then echo "unsorted"; else echo "none"; fi)"
echo ""

# ── Build FASTQ file lists ──────────────────────────────────────────
GEX_R1=$(ls "${GEX_DIR}/"*R1*.fastq.gz | sort | paste -sd, -)
GEX_R2=$(ls "${GEX_DIR}/"*R2*.fastq.gz | sort | paste -sd, -)

echo "GEX lanes:   $(echo "${GEX_R1}" | tr ',' '\n' | wc -l)"
echo "gRNA lanes:  $(ls "${GRNA_DIR}/"*R1*.fastq.gz | wc -l)"
echo "LARRY lanes: $(ls "${LARRY_DIR}/"*R1*.fastq.gz | wc -l)"
echo ""

# ── Create multi-config (3-library with per-library chemistry) ───────
MULTI_CONFIG="${OUTDIR}/multi_config.csv"
cat > "${MULTI_CONFIG}" <<EOF
[libraries]
fastqs,sample,library_type,feature_types,star_chemistry,star_whitelist,star_feature_ref,star_library_id,star_max_hamming
${GEX_DIR},DE_30KO,Gene Expression,Gene Expression,TRU,,,gex_de,
${GRNA_DIR},DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,${GRNA_WHITELIST},${GRNA_FEATURE_REF},grna_de,1
${LARRY_DIR},DE_30KO,Custom,Custom,TRU,${LARRY_WHITELIST},${LARRY_FEATURE_REF},larry_de,1
EOF

echo "Multi-config: ${MULTI_CONFIG}"
cat "${MULTI_CONFIG}"
echo ""

# ── Build STAR command ───────────────────────────────────────────────
CMD=(
  "${STAR_BIN}"
  --runThreadN "${THREADS}"
  --genomeDir "${GENOME_DIR}"
  --readFilesIn "${GEX_R2}" "${GEX_R1}"
  --readFilesCommand zcat
  --outFileNamePrefix "${OUTDIR}/"
  --clip3pPolyG yes
  --soloType CB_UMI_Simple
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12
  --soloBarcodeReadLength 0
  --soloCBwhitelist "${SOLO_CB_WHITELIST}"
  --soloFeatures Gene GeneFull
  --soloUMIdedup 1MM_CR
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
  --soloCellFilter EmptyDrops_CR
  --soloUMIfiltering MultiGeneUMI_CR
  --soloMultiMappers Rescue
  --soloCbUbRequireTogether no
  --soloCrGexFeature GeneFull
  --soloCrMultimapRescue yes
  --pfMultiConfig "${MULTI_CONFIG}"
  --crChemistry auto
  --crMinUmi "${MIN_UMI}"
  --crAssignMaxHamming 5
  --crAssignConsumerThreads -1
  --crAssignSearchThreads 1
  --crAssignSkipQcOutputs 1
  --defaultCrCompat yes
  --dynamicThreadInterface 1
  --dynamicThreadConstMapPermits "${THREADS}"
)

if [[ "${WITH_BAM}" == "1" ]]; then
  CMD+=(--outSAMtype BAM Unsorted)
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

echo "Running STAR (3-library mode)..."
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
dataset=MSK_30polyKO
chemistry=mixed_TRU_NXT
libraries=3 (GEX + gRNA + LARRY)
star_cells=${STAR_CELLS}
wall_seconds=${ELAPSED}
wall_minutes=${ELAPSED_MIN}
threads=${THREADS}
bam=$(if [[ "${WITH_BAM}" == "1" ]]; then echo "unsorted"; else echo "none"; fi)
dynamic_permits=${THREADS}
bootstrap_tiered_hash=yes
multimap_rescue=yes
outdir=${OUTDIR}
EOF

echo ""
echo "Summary: ${OUTDIR}/BENCHMARK_SUMMARY.txt"
echo "Done."
