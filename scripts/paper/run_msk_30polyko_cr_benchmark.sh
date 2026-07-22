#!/usr/bin/env bash
# Paper benchmark: MSK 30polyKO 3-library Perturb-seq — CellRanger 9 reference run.
# Companion to scripts/paper/run_msk_30polyko_benchmark.sh (the STAR side).
#
# Dataset:   MSK scRNAseq 30polyKO
# Chemistry: GEX (TRU), gRNA/PolyIII (NXT), LARRY (TRU)
# Libraries: GEX + CRISPR Guide Capture (PolyIII) + Custom (LARRY)
#
# CellRanger 9 cannot combine 3 libraries in one `multi` invocation, so we run
# `cellranger multi` twice against the same staged fastqs:
#   1) GEX + PolyIII (CRISPR Guide Capture)
#   2) GEX + LARRY   (Custom)
#
# Both runs are run sequentially (per project policy: never parallel benchmarks).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# ── Configurable paths ──────────────────────────────────────────────
CR_BIN="${CR_BIN:-/home/lhhung/cellranger-9.0.1/bin/cellranger}"
MSK_ROOT="${MSK_ROOT:-/storage/MSK-perturb-comparison}"
FASTQ_ROOT="${MSK_FASTQ_ROOT:-${MSK_ROOT}/msk30ko_full_3lib_20260304_095911/fastqs}"
GEX_DIR="${MSK_GEX_DIR:-${FASTQ_ROOT}/mRNA}"
GRNA_DIR="${MSK_GRNA_DIR:-${FASTQ_ROOT}/PolyIII}"
LARRY_DIR="${MSK_LARRY_DIR:-${FASTQ_ROOT}/LARRY}"

# CR transcriptome (Cell Ranger-formatted GRCh38-2024-A).
# Use the stock 10x-built ref. Custom-built indexes (e.g.
# autoindex_110_44/refdata-gex-GRCh38-autoindex11044-cellranger) may include
# newer STAR `genomeParameters.txt` keys (e.g. `genomeType`) that CR9's
# bundled STAR rejects.
CR_TX_REF="${MSK_CR_TX_REF:-/storage/A375-CR-9.01/refdata-gex-GRCh38-2024-A}"

# Feature references (must use feature_type compatible with CR's `multi` config)
# gRNA feature ref must declare feature_type=CRISPR Guide Capture (CR rejects
# `Custom` here since multi.csv routes the library through CRISPR guide calling).
GRNA_FEATURE_REF="${MSK_GRNA_REF:-/mnt/pikachu/MSK-whitelists/ref_feature_geneBC_crispr.csv}"
LARRY_FEATURE_REF="${MSK_LARRY_REF:-/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv}"

# Sample IDs as referenced in fastq filenames (CR's `fastq_id` field)
GEX_SAMPLE_ID="${MSK_GEX_SAMPLE_ID:-mRNA}"
GRNA_SAMPLE_ID="${MSK_GRNA_SAMPLE_ID:-PolyIII}"
LARRY_SAMPLE_ID="${MSK_LARRY_SAMPLE_ID:-LARRY}"

OUTDIR="${MSK_OUTDIR:-${MSK_ROOT}/cr_paper_bench_$(date +%Y%m%d_%H%M%S)}"
THREADS="${MSK_THREADS:-32}"
LOCALMEM_GB="${MSK_LOCALMEM_GB:-110}"
WITH_BAM="${MSK_WITH_BAM:-0}"

# Which CR runs to perform; either or both.
RUN_GRNA="${MSK_CR_RUN_GRNA:-1}"
RUN_LARRY="${MSK_CR_RUN_LARRY:-1}"

# Disk-space preflight: CR `multi` keeps full pipestance state during the run.
# Final output dirs are typically ~5-7 GB each; transient working state is much larger.
REQUIRED_GB_NO_BAM="${MSK_CR_REQUIRED_GB_NO_BAM:-80}"
REQUIRED_GB_WITH_BAM="${MSK_CR_REQUIRED_GB_WITH_BAM:-300}"

require_free_gb() {
  local dir="$1"
  local need_gb="$2"
  local avail_kb
  avail_kb=$(df -Pk "${dir}" | awk 'NR==2 {print $4}')
  local avail_gb=$(( avail_kb / 1024 / 1024 ))
  if (( avail_gb < need_gb )); then
    echo "ERROR: Not enough free space under ${dir}: need >= ${need_gb}GB, have ${avail_gb}GB." >&2
    echo "SOLUTION: choose a different --outdir on a larger filesystem, or run --no-bam." >&2
    exit 1
  fi
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Runs CellRanger 9 \`multi\` twice against the MSK 30polyKO fastqs:
  - cr_gex_grna   (GEX + PolyIII as CRISPR Guide Capture)
  - cr_gex_larry  (GEX + LARRY  as Custom)

Options:
  --outdir DIR         Output directory (auto-generated if omitted)
  --threads N          CR --localcores (default: ${THREADS})
  --localmem GB        CR --localmem GB (default: ${LOCALMEM_GB})
  --no-bam             create-bam,false (default)
  --with-bam           create-bam,true
  --grna-only          Only run the gRNA + GEX config
  --larry-only         Only run the LARRY + GEX config
  --cr-bin PATH        CellRanger binary (default: ${CR_BIN})
  -h, --help           Show help

Environment overrides (selected):
  MSK_FASTQ_ROOT, MSK_GEX_DIR, MSK_GRNA_DIR, MSK_LARRY_DIR
  MSK_GEX_SAMPLE_ID, MSK_GRNA_SAMPLE_ID, MSK_LARRY_SAMPLE_ID
  MSK_CR_TX_REF, MSK_GRNA_REF, MSK_LARRY_REF, MSK_OUTDIR
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir)      OUTDIR="$2"; shift 2 ;;
    --threads)     THREADS="$2"; shift 2 ;;
    --localmem)    LOCALMEM_GB="$2"; shift 2 ;;
    --no-bam)      WITH_BAM=0; shift ;;
    --with-bam)    WITH_BAM=1; shift ;;
    --grna-only)   RUN_GRNA=1; RUN_LARRY=0; shift ;;
    --larry-only)  RUN_GRNA=0; RUN_LARRY=1; shift ;;
    --cr-bin)      CR_BIN="$2"; shift 2 ;;
    -h|--help)     usage; exit 0 ;;
    *)             echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

# ── Validate inputs ─────────────────────────────────────────────────
[[ -x "${CR_BIN}" ]] || { echo "ERROR: cellranger binary not found or not executable: ${CR_BIN}" >&2; exit 1; }
[[ -d "${CR_TX_REF}" ]] || { echo "ERROR: Missing CR transcriptome reference: ${CR_TX_REF}" >&2; exit 1; }
[[ -f "${CR_TX_REF}/reference.json" ]] || { echo "ERROR: Not a CR-format ref (missing reference.json): ${CR_TX_REF}" >&2; exit 1; }
[[ -d "${GEX_DIR}" ]] || { echo "ERROR: Missing GEX dir: ${GEX_DIR}" >&2; exit 1; }
if [[ "${RUN_GRNA}" == "1" ]]; then
  [[ -d "${GRNA_DIR}" ]]            || { echo "ERROR: Missing gRNA dir: ${GRNA_DIR}" >&2; exit 1; }
  [[ -f "${GRNA_FEATURE_REF}" ]]    || { echo "ERROR: Missing gRNA feature ref: ${GRNA_FEATURE_REF}" >&2; exit 1; }
fi
if [[ "${RUN_LARRY}" == "1" ]]; then
  [[ -d "${LARRY_DIR}" ]]           || { echo "ERROR: Missing LARRY dir: ${LARRY_DIR}" >&2; exit 1; }
  [[ -f "${LARRY_FEATURE_REF}" ]]   || { echo "ERROR: Missing LARRY feature ref: ${LARRY_FEATURE_REF}" >&2; exit 1; }
fi

mkdir -p "${OUTDIR}"

# Disk space preflight on the output filesystem.
if [[ "${WITH_BAM}" == "1" ]]; then
  require_free_gb "${OUTDIR}" "${REQUIRED_GB_WITH_BAM}"
else
  require_free_gb "${OUTDIR}" "${REQUIRED_GB_NO_BAM}"
fi

CREATE_BAM_VAL=$([[ "${WITH_BAM}" == "1" ]] && echo "true" || echo "false")

echo "=========================================="
echo "Paper Benchmark: MSK 30polyKO — CellRanger 9"
echo "=========================================="
echo "cellranger:   ${CR_BIN}"
echo "Output dir:   ${OUTDIR}"
echo "Cores/mem:    ${THREADS} cores / ${LOCALMEM_GB} GB"
echo "create-bam:   ${CREATE_BAM_VAL}"
echo "Tx ref:       ${CR_TX_REF}"
echo "GEX dir:      ${GEX_DIR}"
[[ "${RUN_GRNA}" == "1" ]]  && echo "gRNA dir:     ${GRNA_DIR}    (ref: ${GRNA_FEATURE_REF})"
[[ "${RUN_LARRY}" == "1" ]] && echo "LARRY dir:    ${LARRY_DIR}   (ref: ${LARRY_FEATURE_REF})"
echo ""

# ── Helpers ─────────────────────────────────────────────────────────
write_multi_csv() {
  local out_csv="$1"
  local feature_ref="$2"
  local feature_dir="$3"
  local feature_sample_id="$4"
  local feature_type="$5"   # "CRISPR Guide Capture" or "Custom"

  cat > "${out_csv}" <<EOF
[gene-expression]
reference,${CR_TX_REF}
create-bam,${CREATE_BAM_VAL}
no-secondary-analysis,true

[feature]
reference,${feature_ref}

[libraries]
fastq_id,fastqs,feature_types
${GEX_SAMPLE_ID},${GEX_DIR},Gene Expression
${feature_sample_id},${feature_dir},${feature_type}
EOF
}

run_cr_multi() {
  local id="$1"
  local csv="$2"
  local run_dir="$3"

  mkdir -p "${run_dir}"
  cd "${run_dir}"

  # Persist exact invocation for paper repro.
  {
    echo '#!/usr/bin/env bash'
    echo 'set -euo pipefail'
    printf '%q ' \
      "${CR_BIN}" multi \
      --id "${id}" \
      --csv "${csv}" \
      --localcores "${THREADS}" \
      --localmem "${LOCALMEM_GB}"
    printf '\n'
  } > "${run_dir}/RUN_COMMAND.sh"
  chmod +x "${run_dir}/RUN_COMMAND.sh"

  echo "[${id}] cellranger multi starting at $(date -Is)"
  local start_sec=$SECONDS
  /usr/bin/time -v -o "${run_dir}/cellranger.time.log" \
    "${CR_BIN}" multi \
      --id "${id}" \
      --csv "${csv}" \
      --localcores "${THREADS}" \
      --localmem "${LOCALMEM_GB}" \
      2>&1 | tee "${run_dir}/cellranger.stdout.log"

  local elapsed=$(( SECONDS - start_sec ))
  local elapsed_min
  elapsed_min=$(awk "BEGIN {printf \"%.1f\", ${elapsed}/60}")
  echo "[${id}] completed in ${elapsed_min} min (${elapsed} s)"
  echo "${elapsed}" > "${run_dir}/wall_seconds.txt"
}

validate_cr_outputs() {
  local id="$1"
  local run_dir="$2"
  local outs="${run_dir}/${id}/outs"
  local pass=true
  local required=(
    "${outs}/per_sample_outs"
    "${outs}/multi"
  )
  for p in "${required[@]}"; do
    if [[ -e "$p" ]]; then
      echo "  OK: ${p#${run_dir}/}"
    else
      echo "  MISSING: ${p}"
      pass=false
    fi
  done
  if ! ${pass}; then
    return 1
  fi
}

# ── Run #1: GEX + PolyIII (CRISPR Guide Capture) ─────────────────────
GRNA_ELAPSED=
LARRY_ELAPSED=

if [[ "${RUN_GRNA}" == "1" ]]; then
  GRNA_RUN_DIR="${OUTDIR}/cr_gex_grna"
  GRNA_CSV="${GRNA_RUN_DIR}/multi.csv"
  mkdir -p "${GRNA_RUN_DIR}"
  write_multi_csv "${GRNA_CSV}" "${GRNA_FEATURE_REF}" "${GRNA_DIR}" "${GRNA_SAMPLE_ID}" "CRISPR Guide Capture"
  echo "── GEX + gRNA config (${GRNA_CSV}) ──"
  cat "${GRNA_CSV}"
  echo ""
  run_cr_multi "cr9_gex_grna" "${GRNA_CSV}" "${GRNA_RUN_DIR}"
  GRNA_ELAPSED=$(<"${GRNA_RUN_DIR}/wall_seconds.txt")
  validate_cr_outputs "cr9_gex_grna" "${GRNA_RUN_DIR}" || true
fi

# ── Run #2: GEX + LARRY (Custom) ────────────────────────────────────
if [[ "${RUN_LARRY}" == "1" ]]; then
  LARRY_RUN_DIR="${OUTDIR}/cr_gex_larry"
  LARRY_CSV="${LARRY_RUN_DIR}/multi.csv"
  mkdir -p "${LARRY_RUN_DIR}"
  write_multi_csv "${LARRY_CSV}" "${LARRY_FEATURE_REF}" "${LARRY_DIR}" "${LARRY_SAMPLE_ID}" "Custom"
  echo "── GEX + LARRY config (${LARRY_CSV}) ──"
  cat "${LARRY_CSV}"
  echo ""
  run_cr_multi "cr9_gex_larry" "${LARRY_CSV}" "${LARRY_RUN_DIR}"
  LARRY_ELAPSED=$(<"${LARRY_RUN_DIR}/wall_seconds.txt")
  validate_cr_outputs "cr9_gex_larry" "${LARRY_RUN_DIR}" || true
fi

# ── Combined summary ────────────────────────────────────────────────
TOTAL_ELAPSED=$(( ${GRNA_ELAPSED:-0} + ${LARRY_ELAPSED:-0} ))
TOTAL_MIN=$(awk "BEGIN {printf \"%.1f\", ${TOTAL_ELAPSED}/60}")

cat > "${OUTDIR}/BENCHMARK_SUMMARY.txt" <<EOF
dataset=MSK_30polyKO
tool=cellranger-9.0.1
chemistry=mixed_TRU_NXT
runs=$(( (RUN_GRNA == 1) + (RUN_LARRY == 1) ))
gex_grna_wall_seconds=${GRNA_ELAPSED:-NA}
gex_larry_wall_seconds=${LARRY_ELAPSED:-NA}
total_wall_seconds=${TOTAL_ELAPSED}
total_wall_minutes=${TOTAL_MIN}
threads=${THREADS}
localmem_gb=${LOCALMEM_GB}
create_bam=${CREATE_BAM_VAL}
tx_ref=${CR_TX_REF}
grna_feature_ref=${GRNA_FEATURE_REF}
larry_feature_ref=${LARRY_FEATURE_REF}
gex_dir=${GEX_DIR}
grna_dir=${GRNA_DIR}
larry_dir=${LARRY_DIR}
outdir=${OUTDIR}
EOF

echo ""
echo "=========================================="
echo "All requested CR runs complete"
echo "Total wall time: ${TOTAL_MIN} min (${TOTAL_ELAPSED} s)"
echo "Summary: ${OUTDIR}/BENCHMARK_SUMMARY.txt"
echo "=========================================="
