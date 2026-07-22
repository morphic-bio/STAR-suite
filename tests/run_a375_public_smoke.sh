#!/usr/bin/env bash
set -euo pipefail
# Smoke test using publicly downloadable A375 CRISPR fixtures.
#
# Tests the full CR-compat pipeline (GEX + CRISPR multi-library) on a small
# downsampled tier of the 10x A375 public dataset.  Validates:
#   - STAR runs to completion without crashing
#   - Solo/Gene (or GeneFull) filtered/raw MEX directories are created
#   - CR-compat outs/filtered_feature_bc_matrix has both GEX and CRISPR counts
#   - CRISPR protospacer_calls_per_cell.csv is generated
#   - (optionally) BAM output with CB/UB tags
#   - (optionally) Velocyto spliced/unspliced matrices
#
# Prerequisites:
#   - STAR binary (builds automatically if not found)
#   - Genome index at CR_GENOME_DIR (default: /storage/autoindex_110_44/bulk_index)
#   - Downloaded fixture (auto-downloads if missing)
#   - Whitelist file (3M-5pgex-jan-2023.txt for 5', 3M-february-2018.txt for 3')
#
# Usage:
#   tests/run_a375_public_smoke.sh [options]
#
# Options:
#   --5prime             Use 5' GEM-X fixture (default)
#   --3prime             Use 3' v3 NextGEM fixture (enables Velocyto + BAM)
#   --downsample N       Reads per FASTQ (default: 100000)
#   --threads N          STAR threads (default: 4)
#   --no-download        Fail instead of auto-downloading fixture
#   --write-bam          Force BAM output (always on for 3')
#   --velocyto           Enable Velocyto mode (always on for 3')
#   --outdir DIR         Output directory (default: /tmp/a375_public_smoke_<timestamp>)
#   -h, --help           Show help

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Source fixture environment
[[ -f "${SCRIPT_DIR}/external_fixtures_env.sh" ]] && source "${SCRIPT_DIR}/external_fixtures_env.sh"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
GENOME_DIR="${CR_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"

MODE="5prime"
DOWNSAMPLE=100000
THREADS=4
AUTO_DOWNLOAD=1
WRITE_BAM=0
VELOCYTO=0
OUTDIR=""

usage() { sed -n '2,/^$/s/^# \?//p' "$0"; }

log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"; }
die() { echo "FAIL: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --5prime)      MODE="5prime";       shift ;;
    --3prime)      MODE="3prime";       shift ;;
    --downsample)  DOWNSAMPLE="$2";     shift 2 ;;
    --threads)     THREADS="$2";        shift 2 ;;
    --no-download) AUTO_DOWNLOAD=0;     shift ;;
    --write-bam)   WRITE_BAM=1;         shift ;;
    --velocyto)    VELOCYTO=1;          shift ;;
    --outdir)      OUTDIR="$2";         shift 2 ;;
    -h|--help)     usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

# 3' mode implies Velocyto + BAM
if [[ "${MODE}" == "3prime" ]]; then
  VELOCYTO=1
  WRITE_BAM=1
fi

TIMESTAMP="$(date +'%Y%m%d_%H%M%S')"
OUTDIR="${OUTDIR:-/tmp/a375_public_smoke_${MODE}_${TIMESTAMP}}"

# ── Resolve STAR binary ─────────────────────────────────────────────
if [[ ! -x "${STAR_BIN}" ]]; then
  log "STAR binary not found at ${STAR_BIN}, building..."
  make -C "${REPO_ROOT}" core -j"$(nproc)" 2>&1 | tail -3
fi
[[ -x "${STAR_BIN}" ]] || die "STAR binary not found: ${STAR_BIN}"

# ── Resolve fixture paths ───────────────────────────────────────────
if [[ "${MODE}" == "5prime" ]]; then
  FIXTURE_ROOT="${A375_5P_FIXTURE_ROOT:-/tmp/a375_5prime_fixture}"
  DOWNLOAD_SCRIPT="${REPO_ROOT}/scripts/download_a375_5prime_fixture.sh"
  FEATURE_REF="${A375_5P_FEATURE_REF:-${FIXTURE_ROOT}/1k_CRISPR_5p_gemx_Multiplex_count_feature_reference.csv}"
  WHITELIST="${CR_MULTI_WHITELIST:-/storage/A375/3M-5pgex-jan-2023.txt}"
  TIER_DIR="${A375_5P_DOWNSAMPLE_DIR:-${FIXTURE_ROOT}/downsampled_${DOWNSAMPLE}}"
  CB_START=1; CB_LEN=16; UMI_START=17; UMI_LEN=12
  SOLO_STRAND="Unstranded"
  SOLO_FEATURE="GeneFull"
else
  FIXTURE_ROOT="${A375_3P_FIXTURE_ROOT:-/tmp/a375_3prime_fixture}"
  DOWNLOAD_SCRIPT="${REPO_ROOT}/scripts/download_a375_3prime_fixture.sh"
  FEATURE_REF="${A375_3P_FEATURE_REF:-${FIXTURE_ROOT}/SC3_v3_NextGem_DI_CRISPR_10K_feature_ref.csv}"
  WHITELIST="${A375_3P_WHITELIST:-/storage/A375/3M-february-2018.txt}"
  TIER_DIR="${A375_3P_DOWNSAMPLE_DIR:-${FIXTURE_ROOT}/downsampled_${DOWNSAMPLE}}"
  # 3' v3: CB 1-16, UMI 17-28 (12bp)
  CB_START=1; CB_LEN=16; UMI_START=17; UMI_LEN=12
  SOLO_STRAND="Forward"
  SOLO_FEATURE="GeneFull_Ex50pAS"
fi

GEX_DIR="${TIER_DIR}/gex"
CRISPR_DIR="${TIER_DIR}/crispr"

# ── Auto-download fixture if needed ─────────────────────────────────
if [[ ! -d "${GEX_DIR}" || ! -d "${CRISPR_DIR}" ]]; then
  if [[ "${AUTO_DOWNLOAD}" -eq 1 ]]; then
    log "Fixture not found at ${TIER_DIR}, downloading..."
    bash "${DOWNLOAD_SCRIPT}" --outdir "${FIXTURE_ROOT}" --downsample "${DOWNSAMPLE}"
  else
    die "Fixture not found at ${TIER_DIR}. Run ${DOWNLOAD_SCRIPT} first."
  fi
fi

# ── Validate prerequisites ──────────────────────────────────────────
[[ -d "${GEX_DIR}" ]]      || die "GEX dir not found: ${GEX_DIR}"
[[ -d "${CRISPR_DIR}" ]]   || die "CRISPR dir not found: ${CRISPR_DIR}"
[[ -f "${FEATURE_REF}" ]]  || die "Feature ref not found: ${FEATURE_REF}"
[[ -f "${WHITELIST}" ]]    || die "Whitelist not found: ${WHITELIST}"
[[ -d "${GENOME_DIR}" ]]   || die "Genome dir not found: ${GENOME_DIR}"

# ── Build multi-config CSV ──────────────────────────────────────────
mkdir -p "${OUTDIR}"
MULTI_CONFIG="${OUTDIR}/multi_config.csv"
cat > "${MULTI_CONFIG}" <<EOF
[libraries]
fastqs,sample,library_type,feature_types
${GEX_DIR},A375,Gene Expression,Gene Expression
${CRISPR_DIR},A375,CRISPR Guide Capture,CRISPR Guide Capture

[feature]
ref,${FEATURE_REF}
EOF

# ── Build STAR args ─────────────────────────────────────────────────
R1_FILES=$(ls "${GEX_DIR}"/*R1*.fastq.gz 2>/dev/null | paste -sd, -)
R2_FILES=$(ls "${GEX_DIR}"/*R2*.fastq.gz 2>/dev/null | paste -sd, -)
[[ -n "${R1_FILES}" ]] || die "No R1 FASTQs in ${GEX_DIR}"
[[ -n "${R2_FILES}" ]] || die "No R2 FASTQs in ${GEX_DIR}"

STAR_ARGS=(
  --runThreadN "${THREADS}"
  --genomeDir "${GENOME_DIR}"
  --readFilesIn "${R2_FILES}" "${R1_FILES}"
  --readFilesCommand zcat
  --outFileNamePrefix "${OUTDIR}/"
  --clipAdapterType CellRanger4
  --alignEndsType Local
  --chimSegmentMin 1000000
  --soloType CB_UMI_Simple
  --soloCBstart "${CB_START}" --soloCBlen "${CB_LEN}"
  --soloUMIstart "${UMI_START}" --soloUMIlen "${UMI_LEN}"
  --soloBarcodeReadLength 0
  --soloCBwhitelist "${WHITELIST}"
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
  --soloUMIdedup 1MM_CR
  --soloUMIfiltering MultiGeneUMI_CR
  --soloStrand "${SOLO_STRAND}"
  --soloFeatures "${SOLO_FEATURE}"
  --soloCellFilter EmptyDrops_CR
  --soloCrGexFeature auto
  --pfMultiConfig "${MULTI_CONFIG}"
  --crFeatureRef "${FEATURE_REF}"
  --crWhitelist "${WHITELIST}"
)

if [[ "${WRITE_BAM}" -eq 1 ]]; then
  STAR_ARGS+=(
    --outSAMtype BAM SortedByCoordinate
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ
  )
else
  STAR_ARGS+=(--outSAMtype None)
fi

if [[ "${VELOCYTO}" -eq 1 ]]; then
  STAR_ARGS+=(--soloFeatures GeneFull_Ex50pAS Velocyto)
fi

# Append any extra args from harness
if [[ -n "${STAR_EXTRA_ARGS:-}" ]]; then
  # shellcheck disable=SC2206
  STAR_ARGS+=(${STAR_EXTRA_ARGS})
fi

# ── Run STAR ────────────────────────────────────────────────────────
log "=== A375 Public Smoke Test (${MODE}) ==="
log "STAR:       ${STAR_BIN}"
log "Fixture:    ${TIER_DIR} (${DOWNSAMPLE} reads)"
log "Genome:     ${GENOME_DIR}"
log "Output:     ${OUTDIR}"
log "Velocyto:   ${VELOCYTO}"
log "BAM:        ${WRITE_BAM}"

set +e
"${STAR_BIN}" "${STAR_ARGS[@]}" > "${OUTDIR}/star_smoke.log" 2>&1
STAR_EXIT=$?
set -e

if [[ "${STAR_EXIT}" -ne 0 ]]; then
  die "STAR exited with code ${STAR_EXIT}. Log: ${OUTDIR}/star_smoke.log"
fi
log "STAR completed successfully"

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

check_nonempty() {
  local desc="$1" path="$2"
  if [[ -s "${path}" ]]; then
    log "  PASS: ${desc}"
    ((++PASS))
  else
    log "  FAIL: ${desc} — missing or empty: ${path}"
    ((++FAIL))
  fi
}

log "--- Output validation ---"

# CR-compat outs
check "CR-compat filtered MEX" "${OUTDIR}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
check "CR-compat raw MEX"      "${OUTDIR}/outs/raw_feature_bc_matrix/matrix.mtx.gz"
check "CRISPR calls CSV"       "${OUTDIR}/outs/crispr_analysis/protospacer_calls_per_cell.csv"
check "CRISPR calls summary"   "${OUTDIR}/outs/crispr_analysis/protospacer_calls_summary.csv"

# Solo output
SOLO_DIR="${OUTDIR}/Solo.out"
for feature_dir in "${SOLO_DIR}"/*/; do
  feat="$(basename "${feature_dir}")"
  check "Solo ${feat} raw matrix"      "${feature_dir}/raw/matrix.mtx"
  check "Solo ${feat} filtered matrix" "${feature_dir}/filtered/matrix.mtx"
done

# BAM
if [[ "${WRITE_BAM}" -eq 1 ]]; then
  check "BAM output" "${OUTDIR}/Aligned.sortedByCoord.out.bam"
fi

# Velocyto
if [[ "${VELOCYTO}" -eq 1 ]]; then
  check "Velocyto raw matrix" "${OUTDIR}/outs/raw_velocyto_feature_bc_matrix/spliced.mtx.gz"
fi

# Content validation: verify filtered MEX has both GEX and CRISPR features
FILTERED_MEX="${OUTDIR}/outs/filtered_feature_bc_matrix"
if [[ -f "${FILTERED_MEX}/matrix.mtx.gz" ]]; then
  python3 - "${FILTERED_MEX}" <<'PYEOF' && ((++PASS)) || ((++FAIL))
import gzip, sys, os
mex_dir = sys.argv[1]
features_path = os.path.join(mex_dir, "features.tsv.gz")
if not os.path.exists(features_path):
    print("  FAIL: features.tsv.gz missing")
    sys.exit(1)
types_seen = set()
with gzip.open(features_path, "rt") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 3:
            types_seen.add(parts[2])
has_gex = "Gene Expression" in types_seen
has_crispr = "CRISPR Guide Capture" in types_seen
if has_gex and has_crispr:
    print(f"  PASS: filtered MEX has both GEX and CRISPR features ({types_seen})")
    sys.exit(0)
else:
    print(f"  FAIL: filtered MEX feature types: {types_seen} (expected GEX + CRISPR)")
    sys.exit(1)
PYEOF
fi

# ── Summary ─────────────────────────────────────────────────────────
log "--- Results: ${PASS} passed, ${FAIL} failed ---"

if [[ "${FAIL}" -gt 0 ]]; then
  die "${FAIL} validation(s) failed"
fi

log "ALL PASS: A375 ${MODE} public smoke test"
