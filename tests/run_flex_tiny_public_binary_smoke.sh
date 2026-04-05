#!/usr/bin/env bash
set -euo pipefail
# Binary-oriented Flex smoke test.
#
# Uses the same public tiny 10x fixture as run_flex_tiny_public_smoke.sh but
# defaults to --outSAMtype None (no BAM) for release-gate validation.
#
# Tier A: public-download, container-safe, release-gate candidate.
#
# Prerequisites:
#   - STAR binary (env STAR_BIN, or in-tree fallback)
#   - git, git-lfs, python3 (for fixture fetch and generation)
#
# Usage:
#   tests/run_flex_tiny_public_binary_smoke.sh [options]
#
# Options:
#   --threads N          STAR threads (default: 2)
#   --read-limit N       Subsample reads (default: 2000)
#   --write-bam sorted   Write sorted BAM and validate
#   --write-bam unsorted Write unsorted BAM and validate
#   --outdir DIR         Output directory
#   -h, --help           Show help

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
THREADS="${THREADS:-2}"
READ_LIMIT="${READ_LIMIT:-2000}"
BAM_MODE="none"
OUTDIR=""

die() { echo "FAIL: $*" >&2; exit 1; }
log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "missing required command: $1"
}

usage() { sed -n '2,/^$/s/^# \?//p' "$0"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --threads)    THREADS="$2";    shift 2 ;;
    --read-limit) READ_LIMIT="$2"; shift 2 ;;
    --write-bam)  BAM_MODE="$2";   shift 2 ;;
    --outdir)     OUTDIR="$2";     shift 2 ;;
    -h|--help)    usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

TIMESTAMP="$(date +'%Y%m%d_%H%M%S')"
OUTDIR="${OUTDIR:-/tmp/flex_binary_smoke_${TIMESTAMP}_$$}"

# ── Prerequisites ───────────────────────────────────────────────────
need_cmd git
need_cmd git-lfs
need_cmd python3
[[ -x "${STAR_BIN}" ]] || die "STAR binary not found: ${STAR_BIN}"

case "${BAM_MODE}" in
  none|sorted|unsorted) ;;
  *) die "--write-bam must be none, sorted, or unsorted" ;;
esac

if [[ "${BAM_MODE}" != "none" ]]; then
  need_cmd samtools
fi

mkdir -p "${OUTDIR}"
log "=== Flex Binary Smoke Test ==="
log "STAR:       ${STAR_BIN}"
log "Threads:    ${THREADS}"
log "Read limit: ${READ_LIMIT}"
log "BAM mode:   ${BAM_MODE}"
log "Outdir:     ${OUTDIR}"

# ── Fetch and generate tiny Flex fixture ────────────────────────────
DATA_DIR="${OUTDIR}/data"
TINY_ROOT="${DATA_DIR}/public_10x_tiny"
ASSET_DIR="${OUTDIR}/assets"
FILTERED_REF_WORK="${OUTDIR}/refwork"
INDEX_DIR="${OUTDIR}/star_index"
OUT_BASE="${OUTDIR}/run_root"
RUN_ID="run"

DEMO_ROOT="${OUTDIR}/demo"
CACHE_ROOT="${DEMO_ROOT}/cache"
INDEX_ROOT_DIR="${DEMO_ROOT}/indices"
RUN_ROOT_DIR="${DEMO_ROOT}/runs"
mkdir -p "${OUTDIR}" "${DATA_DIR}" "${CACHE_ROOT}" "${INDEX_ROOT_DIR}" "${RUN_ROOT_DIR}"

log "Fetching public tiny 10x fixture..."
DEMO_ROOT="${DEMO_ROOT}" \
CACHE_ROOT="${CACHE_ROOT}" \
DATA_ROOT="${DATA_DIR}" \
INDEX_ROOT="${INDEX_ROOT_DIR}" \
RUN_ROOT="${RUN_ROOT_DIR}" \
  "${REPO_ROOT}/scripts/codespaces/fetch_public_10x_tiny.sh" --outdir "${TINY_ROOT}"

TINY_FASTQ_DIR="${TINY_ROOT}/cellranger-tiny-fastq"
TINY_REF_ROOT="${TINY_ROOT}/cellranger-tiny-ref"
[[ -d "${TINY_FASTQ_DIR}" ]] || die "missing tiny FASTQ dir: ${TINY_FASTQ_DIR}"
[[ -d "${TINY_REF_ROOT}" ]] || die "missing tiny ref dir: ${TINY_REF_ROOT}"

log "Generating tiny Flex demo assets..."
python3 "${REPO_ROOT}/scripts/codespaces/generate_tiny_flex_demo.py" \
  --tiny-fastq-dir "${TINY_FASTQ_DIR}" \
  --tiny-ref-root "${TINY_REF_ROOT}" \
  --outdir "${ASSET_DIR}" \
  --read-limit "${READ_LIMIT}"

WHITELIST="${ASSET_DIR}/whitelist.txt"
CONFIG="${ASSET_DIR}/config.csv"
PROBE_CATALOG="${ASSET_DIR}/sample_probe_catalog.tsv"
[[ -f "${WHITELIST}" ]] || die "missing whitelist: ${WHITELIST}"
[[ -f "${CONFIG}" ]]    || die "missing config: ${CONFIG}"
[[ -f "${PROBE_CATALOG}" ]] || die "missing sample probe catalog: ${PROBE_CATALOG}"

GTF_PATH="${TINY_REF_ROOT}/genes/genes.gtf"
if [[ ! -f "${GTF_PATH}" ]]; then
  GTF_PATH="${TINY_REF_ROOT}/genes/genes.gtf.gz"
fi
[[ -f "${GTF_PATH}" ]] || die "missing tiny GTF under ${TINY_REF_ROOT}/genes"

log "Building filtered reference..."
"${REPO_ROOT}/flex/scripts/build_filtered_reference.sh" \
  --probe-set "${ASSET_DIR}/probe_set.csv" \
  --base-fasta "${TINY_REF_ROOT}/fasta/genome.fa" \
  --base-gtf "${GTF_PATH}" \
  --work-dir "${FILTERED_REF_WORK}" \
  --skip-filter \
  --quiet

log "Building STAR index..."
"${REPO_ROOT}/flex/scripts/make_filtered_star_index.sh" \
  --filtered-reference "${FILTERED_REF_WORK}/filtered_reference" \
  --output-dir "${INDEX_DIR}" \
  --threads "${THREADS}" \
  --sa-index-n-bases 11 \
  --star-bin "${STAR_BIN}"

# Note: run_flex_cr_config.sh hardcodes --outSAMtype BAM Unsorted. The Flex
# wrapper requires BAM output for its pipeline. The --write-bam flag here
# controls validation expectations, not the wrapper's BAM behavior.

log "Running Flex pipeline..."
"${REPO_ROOT}/scripts/run_flex_cr_config.sh" \
  --cr-config "${CONFIG}" \
  --genome-dir "${INDEX_DIR}" \
  --cb-whitelist "${WHITELIST}" \
  --solo-cb-start 1 \
  --solo-cb-len 16 \
  --solo-umi-start 17 \
  --solo-umi-len 10 \
  --sample-probe-catalog "${PROBE_CATALOG}" \
  --sample-probe-offset 68 \
  --out-base "${OUT_BASE}" \
  --run-id "${RUN_ID}" \
  --threads "${THREADS}"

RUN_ROOT="${OUT_BASE}/${RUN_ID}"

# ── Validate outputs ────────────────────────────────────────────────
PASS=0
FAIL=0

check() {
  local desc="$1" path="$2"
  if [[ -e "${path}" ]]; then
    log "  PASS: ${desc}"
    ((++PASS))
  else
    log "  FAIL: ${desc} — not found: ${path}"
    ((++FAIL))
  fi
}

log "--- Validating outputs ---"

# Core Flex outputs
check "RUN_MANIFEST.txt"    "${RUN_ROOT}/RUN_MANIFEST.txt"
check "Log.final.out"       "${RUN_ROOT}/Log.final.out"
check "Log.out"             "${RUN_ROOT}/Log.out"
check "Barcodes.stats"      "${RUN_ROOT}/Solo.out/Barcodes.stats"

# Flex-specific log markers
if grep -qF "Enabled Flex pipeline with production defaults" "${RUN_ROOT}/Log.out" 2>/dev/null; then
  log "  PASS: Flex pipeline enablement in Log.out"
  ((++PASS))
else
  log "  FAIL: missing Flex pipeline enablement in Log.out"
  ((++FAIL))
fi

if grep -qF "SampleDetector initialized successfully" "${RUN_ROOT}/Log.out" 2>/dev/null; then
  log "  PASS: SampleDetector init in Log.out"
  ((++PASS))
else
  log "  FAIL: missing SampleDetector init in Log.out"
  ((++FAIL))
fi

# Manifest content
if grep -qF "sample_probe_offset=68" "${RUN_ROOT}/RUN_MANIFEST.txt" 2>/dev/null; then
  log "  PASS: manifest has sample_probe_offset"
  ((++PASS))
else
  log "  FAIL: manifest missing sample_probe_offset"
  ((++FAIL))
fi

# Flex wrapper always produces unsorted BAM
check "Aligned.out.bam" "${RUN_ROOT}/Aligned.out.bam"

log "--- Results: ${PASS} passed, ${FAIL} failed ---"

if [[ "${FAIL}" -gt 0 ]]; then
  die "${FAIL} validation(s) failed. Outputs at: ${OUTDIR}"
fi

log "PASS: Flex binary smoke test (BAM mode: ${BAM_MODE})"
log "Workdir: ${OUTDIR}"
