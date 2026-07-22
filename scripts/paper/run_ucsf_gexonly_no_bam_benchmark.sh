#!/usr/bin/env bash
# UCSF EBs2_2 full GEX lanes, Solo only (no pfMulti / no guide FASTQs), --outSAMtype None.
#
# Parity: keep in sync with docs/HANDOFF_UCSF_GEXONLY_BENCHMARK_WRAPPER_20260326.md
# Invariants (do not regress):
#   - Exactly one of --historical-vanilla | --modern-optimized (MODE_FLAG_COUNT must be 1).
#   - Never mkdir "${OUTDIR}/tmp" before STAR; only rm -rf so STAR creates --outTmpDir.
#   - Historical: --readFilesCommand zcat for .fastq.gz on 7a7fb08-class binaries.
#   - Modern: no zcat (native gzip).
#
# Two command families (do not mix flags between them):
#   --historical-vanilla  Known-good Solo surface for STAR built from ~7a7fb08
#                         (matches tests/run_ucsf_solo_gex_100k_benchmark.sh shape;
#                         see docs/HANDOFF_SOLO_OPTIMIZATION_20260324.md).
#   --modern-optimized    CR-compat / bridge path for current Solo-optimized builds:
#                         STAR_SOLO_NONFLEX_HASH_BRIDGE=1 + --soloInlineHashMode yes
#                         plus perturb-style GeneFull flags (Unique multimappers, rescue, Forward, poly-G).
#
# Read compression:
#   --modern-optimized: native .gz via STAR (no --readFilesCommand zcat).
#   --historical-vanilla: STAR ~7a7fb08 expects plain FASTQ from readFilesIn; for .fastq.gz
#     inputs you must use --readFilesCommand zcat (otherwise gzip bytes are parsed as FASTQ).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
GEX_DIR="${UCSF_GEXONLY_GEX_DIR:-/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2/GEX}"
SOLO_CB_WHITELIST="${UCSF_GEXONLY_SOLO_WHITELIST:-/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/3M-february-2018_TRU.txt}"
GENOME_DIR="${UCSF_GEXONLY_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
THREADS="${UCSF_GEXONLY_THREADS:-32}"
OUTDIR=""
MODE=""
MODE_FLAG_COUNT=0

usage() {
  cat <<EOF
Usage: $(basename "$0") --outdir DIR (--historical-vanilla | --modern-optimized) [options]

  Exactly one of:
    --historical-vanilla   CLI compatible with historical STAR (e.g. commit 7a7fb08)
    --modern-optimized     Bridge + inline hash + current CR-compat Solo flags

  --outdir DIR       Output prefix directory (required)
  --star-bin PATH    STAR binary (default: \${REPO_ROOT}/core/legacy/source/STAR)
  --threads N        Thread count (default ${THREADS})

Env overrides: UCSF_GEXONLY_GEX_DIR, UCSF_GEXONLY_SOLO_WHITELIST, UCSF_GEXONLY_GENOME_DIR

Modern mode: no zcat. Historical mode: zcat required for .fastq.gz on 7a7fb08-class builds.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)        OUTDIR="$2"; shift 2 ;;
    --star-bin)      STAR_BIN="$2"; shift 2 ;;
    --threads)       THREADS="$2"; shift 2 ;;
    --historical-vanilla)
      MODE="historical-vanilla"
      MODE_FLAG_COUNT=$((MODE_FLAG_COUNT + 1))
      shift
      ;;
    --modern-optimized)
      MODE="modern-optimized"
      MODE_FLAG_COUNT=$((MODE_FLAG_COUNT + 1))
      shift
      ;;
    -h|--help)       usage; exit 0 ;;
    *)               echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -n "${OUTDIR}" ]] || { echo "ERROR: --outdir required" >&2; usage; exit 1; }
if [[ "${MODE_FLAG_COUNT}" -ne 1 ]]; then
  echo "ERROR: specify exactly one mode flag (--historical-vanilla or --modern-optimized); got count=${MODE_FLAG_COUNT}" >&2
  usage
  exit 1
fi
[[ -n "${MODE}" ]] || { echo "ERROR: internal: MODE empty" >&2; exit 1; }
[[ -x "${STAR_BIN}" ]] || { echo "ERROR: STAR not executable: ${STAR_BIN}" >&2; exit 1; }
[[ -d "${GEX_DIR}" ]] || { echo "ERROR: GEX dir missing: ${GEX_DIR}" >&2; exit 1; }
[[ -f "${SOLO_CB_WHITELIST}" ]] || { echo "ERROR: whitelist missing: ${SOLO_CB_WHITELIST}" >&2; exit 1; }
[[ -d "${GENOME_DIR}" ]] || { echo "ERROR: genome missing: ${GENOME_DIR}" >&2; exit 1; }

mkdir -p "${OUTDIR}"
# Do not pre-create --outTmpDir: some STAR builds (e.g. 7a7fb08) error if it already exists.
rm -rf "${OUTDIR}/tmp"

R1_FILES=$(find -L "${GEX_DIR}" -maxdepth 1 -type f -name '*R1*.fastq.gz' | sort | paste -sd, -)
R2_FILES=$(find -L "${GEX_DIR}" -maxdepth 1 -type f -name '*R2*.fastq.gz' | sort | paste -sd, -)
[[ -n "${R1_FILES}" && -n "${R2_FILES}" ]] || { echo "ERROR: no R1/R2 .fastq.gz under ${GEX_DIR}" >&2; exit 1; }

if [[ "${MODE}" == "historical-vanilla" ]]; then
  unset STAR_SOLO_NONFLEX_HASH_BRIDGE 2>/dev/null || true
  CMD=(
    "${STAR_BIN}"
    --runThreadN "${THREADS}"
    --genomeDir "${GENOME_DIR}"
    --readFilesIn "${R2_FILES}" "${R1_FILES}"
    --readFilesCommand zcat
    --outFileNamePrefix "${OUTDIR}/"
    --outTmpDir "${OUTDIR}/tmp"
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
    --soloMultiMappers Rescue
    --soloCbUbRequireTogether no
    --soloCellFilter EmptyDrops_CR
    --soloStrand Unstranded
    --soloFeatures GeneFull
  )
else
  export STAR_SOLO_NONFLEX_HASH_BRIDGE=1
  CMD=(
    "${STAR_BIN}"
    --runThreadN "${THREADS}"
    --genomeDir "${GENOME_DIR}"
    --readFilesIn "${R2_FILES}" "${R1_FILES}"
    --outFileNamePrefix "${OUTDIR}/"
    --outTmpDir "${OUTDIR}/tmp"
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
    --dynamicThreadInterface 1
    --dynamicThreadConstMapPermits "${THREADS}"
    --dynamicThreadTelemetry 1
    --outSAMtype None
    --soloInlineHashMode yes
  )
fi

{
  echo '#!/usr/bin/env bash'
  echo 'set -euo pipefail'
  [[ "${MODE}" == "modern-optimized" ]] && echo 'export STAR_SOLO_NONFLEX_HASH_BRIDGE=1'
  printf '%q ' "${CMD[@]}"
  echo
} > "${OUTDIR}/RUN_COMMAND.sh"
chmod +x "${OUTDIR}/RUN_COMMAND.sh"

echo "UCSF GEX-only no-BAM | mode=${MODE} | STAR=${STAR_BIN} | out=${OUTDIR}"
START_SEC=$SECONDS
"${CMD[@]}"
ELAPSED=$((SECONDS - START_SEC))
ELAPSED_MIN=$(awk "BEGIN {printf \"%.2f\", ${ELAPSED}/60}")

CELLS=$(wc -l < "${OUTDIR}/Solo.out/GeneFull/filtered/barcodes.tsv")
cat > "${OUTDIR}/BENCHMARK_SUMMARY.txt" <<EOF
dataset=UCSF_EBs2_2_GEX_only
mode=${MODE}
star_cells=${CELLS}
wall_seconds=${ELAPSED}
wall_minutes=${ELAPSED_MIN}
threads=${THREADS}
bam=none
read_decompress=$([[ "${MODE}" == "historical-vanilla" ]] && echo zcat || echo native_gzip)
outdir=${OUTDIR}
EOF

echo "Done: ${ELAPSED}s (${ELAPSED_MIN} min), filtered cells=${CELLS}"
