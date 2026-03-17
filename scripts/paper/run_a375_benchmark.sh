#!/usr/bin/env bash
# Paper benchmark: A375 1k CRISPR 5' GemX (10x public dataset)
# Compares STAR CR-compat mode against CellRanger 9 reference output.
#
# Dataset:  10x Genomics 1k_CRISPR_5p_gemx (SC5P-GEX + SC5P-CRISPR)
# Chemistry: 5' v2 GemX (TRU namespace, single-column whitelist)
# Libraries: GEX + CRISPR Guide Capture (2 libraries)
# Reference outputs: downloaded from 10x website (CellRanger count)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# ── Configurable paths ──────────────────────────────────────────────
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
ROOT="${A375_ROOT:-/storage/A375}"
FASTQ_ROOT="${A375_FASTQ_ROOT:-${ROOT}/fastqs/1k_CRISPR_5p_gemx_fastqs}"
GEX_DIR="${A375_GEX_DIR:-${FASTQ_ROOT}/gex}"
CRISPR_DIR="${A375_CRISPR_DIR:-${FASTQ_ROOT}/crispr}"
FEATURE_REF="${A375_FEATURE_REF:-${ROOT}/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv}"
WHITELIST="${A375_WHITELIST:-${ROOT}/3M-5pgex-jan-2023.txt}"
GENOME_DIR="${A375_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
OUTDIR="${A375_OUTDIR:-/storage/A375/paper_bench_$(date +%Y%m%d_%H%M%S)}"
THREADS="${A375_THREADS:-32}"
MIN_UMI="${A375_MIN_UMI:-10}"

# ── Validate inputs ─────────────────────────────────────────────────
for f in "${STAR_BIN}"; do
  [[ -x "$f" ]] || { echo "ERROR: STAR binary not found: $f" >&2; exit 1; }
done
for f in "${WHITELIST}" "${FEATURE_REF}"; do
  [[ -f "$f" ]] || { echo "ERROR: Missing file: $f" >&2; exit 1; }
done
for d in "${GENOME_DIR}" "${GEX_DIR}" "${CRISPR_DIR}"; do
  [[ -d "$d" ]] || { echo "ERROR: Missing directory: $d" >&2; exit 1; }
done

mkdir -p "${OUTDIR}"

echo "=========================================="
echo "Paper Benchmark: A375 1k CRISPR 5' GemX"
echo "=========================================="
echo "STAR binary:   ${STAR_BIN}"
echo "Output dir:    ${OUTDIR}"
echo "Threads:       ${THREADS}"
echo "Chemistry:     TRU (5' v2 GemX, single-column whitelist)"
echo "Min UMI:       ${MIN_UMI}"
echo ""

# ── Build FASTQ file lists ──────────────────────────────────────────
GEX_R1=$(ls "${GEX_DIR}/"*R1*.fastq.gz | paste -sd, -)
GEX_R2=$(ls "${GEX_DIR}/"*R2*.fastq.gz | paste -sd, -)
CRISPR_R1=$(ls "${CRISPR_DIR}/"*R1*.fastq.gz | paste -sd, -)
CRISPR_R2=$(ls "${CRISPR_DIR}/"*R2*.fastq.gz | paste -sd, -)

# ── Create multi-config ─────────────────────────────────────────────
MULTI_CONFIG="${OUTDIR}/multi_config.csv"
cat > "${MULTI_CONFIG}" <<EOF
[libraries]
fastqs,sample,library_type,feature_types
${GEX_DIR},A375,Gene Expression,Gene Expression
${CRISPR_DIR},A375,CRISPR Guide Capture,CRISPR Guide Capture

[feature]
ref,${FEATURE_REF}
EOF

echo "Multi-config: ${MULTI_CONFIG}"
cat "${MULTI_CONFIG}"
echo ""

# ── Run STAR ─────────────────────────────────────────────────────────
START_SEC=$SECONDS

"${STAR_BIN}" \
  --runThreadN "${THREADS}" \
  --genomeDir "${GENOME_DIR}" \
  --readFilesIn "${GEX_R2},${CRISPR_R2}" "${GEX_R1},${CRISPR_R1}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${OUTDIR}/" \
  --outSAMtype None \
  --clipAdapterType CellRanger4 \
  --alignEndsType Local \
  --chimSegmentMin 1000000 \
  --soloType CB_UMI_Simple \
  --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloCBwhitelist "${WHITELIST}" \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloMultiMappers Unique \
  --soloCellFilter EmptyDrops_CR \
  --soloCbUbRequireTogether no \
  --soloStrand Unstranded \
  --soloFeatures GeneFull \
  --soloCrGexFeature genefull \
  --soloCrMultimapRescue yes \
  --pfMultiConfig "${MULTI_CONFIG}" \
  --crChemistry TRU \
  --crWhitelist "${WHITELIST}" \
  --crFeatureRef "${FEATURE_REF}" \
  --crMinUmi "${MIN_UMI}" \
  --crAssignMaxHamming 1 \
  --crAssignFeatureOffset 0 \
  --crAssignLimitSearch -1 \
  --crAssignMinCounts 0 \
  --crAssignMaxBarcodeMismatches 5 \
  --crAssignFeatureN 0 \
  --crAssignBarcodeN 1 \
  --crAssignConsumerThreads -1 \
  --crAssignSearchThreads 1 \
  --crAssignSkipQcOutputs 1 \
  --dynamicThreadInterface 1 \
  --dynamicThreadConstMapPermits "${THREADS}" \
  --dynamicThreadTelemetry 1

ELAPSED=$(( SECONDS - START_SEC ))
ELAPSED_MIN=$(awk "BEGIN {printf \"%.1f\", ${ELAPSED}/60}")

echo ""
echo "=========================================="
echo "STAR completed in ${ELAPSED_MIN} min (${ELAPSED} s)"
echo "=========================================="

# ── Validate expected outputs ────────────────────────────────────────
PASS=true
for f in \
  "${OUTDIR}/Solo.out/GeneFull/filtered/barcodes.tsv" \
  "${OUTDIR}/Solo.out/GeneFull/filtered/matrix.mtx" \
  "${OUTDIR}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" \
  "${OUTDIR}/outs/filtered_feature_bc_matrix/matrix.mtx.gz" \
  "${OUTDIR}/outs/crispr_analysis/protospacer_calls_per_cell.csv"; do
  if [[ -f "$f" ]]; then
    echo "  OK: $f"
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
  echo ""
  echo "TEST PASSED"
else
  echo ""
  echo "TEST FAILED: missing outputs"
  exit 1
fi

# ── Summary file for downstream parity scripts ───────────────────────
cat > "${OUTDIR}/BENCHMARK_SUMMARY.txt" <<EOF
dataset=A375_1k_CRISPR_5p_gemx
chemistry=TRU
star_cells=${STAR_CELLS}
wall_seconds=${ELAPSED}
wall_minutes=${ELAPSED_MIN}
threads=${THREADS}
bam=none
dynamic_permits=${THREADS}
bootstrap_tiered_hash=yes
multimap_rescue=yes
outdir=${OUTDIR}
EOF

echo ""
echo "Summary: ${OUTDIR}/BENCHMARK_SUMMARY.txt"
echo "Done."
