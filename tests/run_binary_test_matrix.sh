#!/usr/bin/env bash
set -euo pipefail
# Binary test matrix runner.
#
# Runs the STAR-suite binary validation suite against an installed or in-tree
# STAR binary.  Designed for release gating: every test uses public-download
# or synthetic fixtures and can run in a clean Ubuntu 22.04/24.04 container.
#
# Usage:
#   tests/run_binary_test_matrix.sh [options]
#
# Options:
#   --tier a          Run Tier A tests only (default)
#   --tier b          Run Tier B tests only (artifact-backed, host-dependent)
#   --tier all        Run both Tier A and Tier B
#   --threads N       STAR threads for individual tests (default: 4)
#   --outdir DIR      Root output directory (default: /tmp/binary_test_matrix_<ts>)
#   --keep            Keep output on success
#   --stop-on-fail    Stop at first failure
#   -h, --help        Show help
#
# Environment:
#   STAR_BIN          Path to STAR binary (default: in-tree)
#
# Tier A (release-gate, public/download-based):
#   1. scRNA GEX no-BAM smoke          (run_gex_binary_smoke.sh)
#   2. Perturb A375 5' public smoke    (run_a375_public_smoke.sh --5prime)
#   3. Flex tiny public smoke           (run_flex_tiny_public_binary_smoke.sh)
#   4. PE bulk public smoke             (run_public_bulk_pe_smoke.sh)
#   5. SLAM SNP-mask build smoke        (slam/test_snp_mask_build_smoke.sh)
#   6. GEX BAM sorted smoke            (run_gex_binary_smoke.sh --write-bam sorted)
#   7. GEX BAM unsorted smoke          (run_gex_binary_smoke.sh --write-bam unsorted)
#   8. genomeGenerate validation        (run_genome_generate_validation.sh)
#   9. Y-removal unit tests             (unit_test_y_split.sh)
#  10. Y-removal comprehensive          (run_y_removal_comprehensive_test.sh)
#  11. TRU/NXT chemistry autodetect     (autodetect_nxt_tru/test_autodetect.sh)
#  12. TRU/NXT config column parsing    (autodetect_nxt_tru/test_star_chemistry_column.sh)
#  13. Trim QC FASTQ smoke              (run_trim_qc_fastq_smoke.sh)
#  14. Adapter clip synthetic           (run_adapter_clip_synthetic_test.sh)
#  15. Cutadapt 3.2 parity             (test_cutadapt32_parity.sh)
#
# Tier B (artifact-backed, extended parity):
#  16. Sorted BAM CB/UB non-Flex        (run_sorted_bam_cbub_nonflex_test.sh)
#  17. Unsorted BAM CB/UB               (run_unsorted_bam_cbub_test.sh)
#  18. Solo Y-removal smoke             (run_solo_yremove_smoke.sh)
#  19. Flex CR-config smoke             (run_flex_cr_config_smoke.sh)
#  20. SLAM end-to-end fixture          (run_slam_end_to_end_fixture.sh)
#  21. TranscriptVB quick test          (transcriptvb/quick_test.sh)
#  22. TranscriptVB edge cases          (transcriptvb/edge_case_tests.sh)
#  23. Salmon parity                    (transcriptvb/salmon_parity_test.sh)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
TIER="a"
THREADS="${THREADS:-4}"
OUTDIR=""
STOP_ON_FAIL=0

die() { echo "FATAL: $*" >&2; exit 1; }
log() { printf '[%s] %s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')" "$*"; }

usage() { sed -n '2,/^$/s/^# \?//p' "$0"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tier)          TIER="$2";    shift 2 ;;
    --threads)       THREADS="$2"; shift 2 ;;
    --outdir)        OUTDIR="$2";  shift 2 ;;
    --keep)          shift ;;
    --stop-on-fail)  STOP_ON_FAIL=1; shift ;;
    -h|--help)       usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

TIMESTAMP="$(date +'%Y%m%d_%H%M%S')"
OUTDIR="${OUTDIR:-/tmp/binary_test_matrix_${TIMESTAMP}_$$}"
mkdir -p "${OUTDIR}"

[[ -x "${STAR_BIN}" ]] || die "STAR binary not found: ${STAR_BIN}"

export STAR_BIN

log "=========================================="
log "  STAR-suite Binary Test Matrix"
log "=========================================="
log "STAR:    ${STAR_BIN}"
log "Version: $("${STAR_BIN}" --version 2>/dev/null || echo unknown)"
log "Tier:    ${TIER}"
log "Threads: ${THREADS}"
log "Outdir:  ${OUTDIR}"
log ""

TOTAL=0
PASSED=0
FAILED=0
SKIPPED=0
FAILED_NAMES=()

run_test() {
  local name="$1"
  shift
  local cmd=("$@")

  ((++TOTAL))
  log "--- [${TOTAL}] ${name} ---"

  local test_out="${OUTDIR}/${name//\//_}.log"

  set +e
  "${cmd[@]}" > "${test_out}" 2>&1
  local rc=$?
  set -e

  if [[ ${rc} -eq 0 ]]; then
    log "  PASS: ${name}"
    ((++PASSED))
  else
    log "  FAIL: ${name} (exit ${rc}, log: ${test_out})"
    ((++FAILED))
    FAILED_NAMES+=("${name}")
    if [[ ${STOP_ON_FAIL} -eq 1 ]]; then
      log "Stopping on first failure (--stop-on-fail)."
      print_summary
      exit 1
    fi
  fi
}

skip_test() {
  local name="$1"
  local reason="$2"
  ((++TOTAL))
  ((++SKIPPED))
  log "--- [${TOTAL}] ${name} ---"
  log "  SKIP: ${reason}"
}

print_summary() {
  log ""
  log "=========================================="
  log "  Summary: ${PASSED} passed, ${FAILED} failed, ${SKIPPED} skipped / ${TOTAL} total"
  log "=========================================="
  if [[ ${FAILED} -gt 0 ]]; then
    log "Failed tests:"
    for name in "${FAILED_NAMES[@]}"; do
      log "  - ${name}"
    done
  fi
  log "Logs: ${OUTDIR}"
}

# ── Tier A: release-gate tests ──────────────────────────────────────

run_tier_a() {
  log "=== Tier A: Release-gate binary workflow tests ==="

  # 1. scRNA GEX no-BAM smoke
  run_test "gex-binary-smoke-no-bam" \
    bash "${SCRIPT_DIR}/run_gex_binary_smoke.sh" \
      --threads "${THREADS}" \
      --read-limit 5000 \
      --outdir "${OUTDIR}/gex_no_bam"

  # 2. Perturb A375 5' public smoke
  if [[ -d "${CR_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}" ]]; then
    run_test "perturb-a375-5prime" \
      bash "${SCRIPT_DIR}/run_a375_public_smoke.sh" \
        --5prime \
        --downsample 100000 \
        --threads "${THREADS}" \
        --outdir "${OUTDIR}/a375_5prime"
  else
    skip_test "perturb-a375-5prime" "genome dir not available (needs CR_GENOME_DIR)"
  fi

  # 3. Flex tiny public smoke
  run_test "flex-tiny-public-binary" \
    bash "${SCRIPT_DIR}/run_flex_tiny_public_binary_smoke.sh" \
      --threads "${THREADS}" \
      --outdir "${OUTDIR}/flex_tiny"

  # 4. PE bulk public smoke
  #    Requires pre-downloaded SRR4422207 fixture, genome index, salmon, and
  #    trimvalidate. Skip if any prerequisite is missing.
  if [[ -f "${SCRIPT_DIR}/run_public_bulk_pe_smoke.sh" ]]; then
    local pe_genome="${PUBLIC_BULK_GENOME_DIR:-/storage/autoindex_110_44/bulk_index}"
    local pe_fixture="${PUBLIC_BULK_FIXTURE_ROOT:-/tmp/public_bulk_fixture_SRR4422207}"
    local pe_r1="${PUBLIC_BULK_R1:-${pe_fixture}/SRR4422207_1.fastq.gz}"
    if [[ ! -d "${pe_genome}" ]]; then
      skip_test "pe-bulk-public" "genome dir not available (${pe_genome})"
    elif [[ ! -f "${pe_r1}" ]]; then
      skip_test "pe-bulk-public" "fixture not downloaded (${pe_r1})"
    elif ! command -v salmon >/dev/null 2>&1; then
      skip_test "pe-bulk-public" "salmon not installed"
    else
      run_test "pe-bulk-public" \
        bash "${SCRIPT_DIR}/run_public_bulk_pe_smoke.sh"
    fi
  else
    skip_test "pe-bulk-public" "run_public_bulk_pe_smoke.sh not found"
  fi

  # 5. SLAM SNP-mask build smoke
  run_test "slam-snp-mask-build" \
    bash "${SCRIPT_DIR}/slam/test_snp_mask_build_smoke.sh"

  # 6. GEX BAM sorted smoke
  run_test "gex-binary-smoke-bam-sorted" \
    bash "${SCRIPT_DIR}/run_gex_binary_smoke.sh" \
      --threads "${THREADS}" \
      --read-limit 5000 \
      --write-bam sorted \
      --outdir "${OUTDIR}/gex_bam_sorted"

  # 7. GEX BAM unsorted smoke
  run_test "gex-binary-smoke-bam-unsorted" \
    bash "${SCRIPT_DIR}/run_gex_binary_smoke.sh" \
      --threads "${THREADS}" \
      --read-limit 5000 \
      --write-bam unsorted \
      --outdir "${OUTDIR}/gex_bam_unsorted"

  # 8. genomeGenerate validation
  run_test "genome-generate-validation" \
    bash "${SCRIPT_DIR}/run_genome_generate_validation.sh" \
      --threads "${THREADS}" \
      --outdir "${OUTDIR}/genome_generate"

  # 9. Y-removal unit tests
  run_test "y-removal-unit" \
    bash "${SCRIPT_DIR}/unit_test_y_split.sh"

  # 10. Y-removal comprehensive
  local remove_y_bin="${REPO_ROOT}/core/features/yremove_fastq/tools/remove_y_reads/remove_y_reads"
  if [[ -x "${remove_y_bin}" ]]; then
    run_test "y-removal-comprehensive" \
      bash "${SCRIPT_DIR}/run_y_removal_comprehensive_test.sh"
  else
    skip_test "y-removal-comprehensive" "remove_y_reads binary not built"
  fi

  # 11. TRU/NXT chemistry autodetect
  local assign_bin="${REPO_ROOT}/core/features/process_features/assignBarcodes"
  if [[ -x "${assign_bin}" ]]; then
    run_test "tru-nxt-autodetect" \
      bash "${SCRIPT_DIR}/autodetect_nxt_tru/test_autodetect.sh"
  else
    skip_test "tru-nxt-autodetect" "assignBarcodes binary not built"
  fi

  # 12. TRU/NXT config column parsing
  if command -v g++ >/dev/null 2>&1 && [[ -f "${REPO_ROOT}/core/legacy/source/PfMultiConfig.o" ]]; then
    run_test "tru-nxt-config-parsing" \
      bash "${SCRIPT_DIR}/autodetect_nxt_tru/test_star_chemistry_column.sh"
  else
    skip_test "tru-nxt-config-parsing" "g++ or PfMultiConfig.o not available"
  fi

  # 13. Trim QC FASTQ smoke
  run_test "trim-qc-fastq" \
    bash "${SCRIPT_DIR}/run_trim_qc_fastq_smoke.sh"

  # 14. Adapter clip synthetic
  run_test "adapter-clip-synthetic" \
    bash "${SCRIPT_DIR}/run_adapter_clip_synthetic_test.sh" \
      --threads "${THREADS}" \
      --outdir "${OUTDIR}/adapter_clip"

  # 15. Cutadapt 3.2 parity
  local trimvalidate_bin="${REPO_ROOT}/core/features/vbem/tools/trimvalidate/trimvalidate"
  if [[ -x "${trimvalidate_bin}" ]] && [[ -f "${SCRIPT_DIR}/fixtures/cutadapt32_adapter_fixtures.json" ]]; then
    run_test "cutadapt-parity" \
      bash "${SCRIPT_DIR}/test_cutadapt32_parity.sh"
  else
    skip_test "cutadapt-parity" "trimvalidate binary or fixtures not available"
  fi
}

# ── Tier B: extended parity / artifact-backed ───────────────────────

run_tier_b() {
  log "=== Tier B: Extended binary parity tests (artifact-backed) ==="

  local scripts=(
    "sorted-bam-cbub-nonflex:run_sorted_bam_cbub_nonflex_test.sh"
    "unsorted-bam-cbub:run_unsorted_bam_cbub_test.sh"
    "solo-yremove:run_solo_yremove_smoke.sh"
    "flex-cr-config:run_flex_cr_config_smoke.sh"
    "slam-e2e-fixture:run_slam_end_to_end_fixture.sh"
  )

  for entry in "${scripts[@]}"; do
    local name="${entry%%:*}"
    local script="${entry##*:}"
    local path="${SCRIPT_DIR}/${script}"
    if [[ -f "${path}" ]]; then
      run_test "${name}" bash "${path}"
    else
      skip_test "${name}" "${script} not found"
    fi
  done

  # TranscriptVB / Salmon (need pre-built genome index + test reads)
  local vb_genome="${VB_GENOME_DIR:-/tmp/star_vb_test/star_new_index}"
  local vb_reads="/mnt/pikachu/test-datasets-rnaseq/testdata/GSE110004/SRR6357070_1.fastq.gz"

  if [[ -d "${vb_genome}" ]] && [[ -f "${vb_reads}" ]]; then
    run_test "transcriptvb-quick" \
      bash "${SCRIPT_DIR}/transcriptvb/quick_test.sh" \
        "${STAR_BIN}" "${vb_genome}"
    run_test "transcriptvb-edge-cases" \
      bash "${SCRIPT_DIR}/transcriptvb/edge_case_tests.sh" \
        "${STAR_BIN}" "${vb_genome}"
  else
    skip_test "transcriptvb-quick" "VB genome index or GSE110004 reads not available"
    skip_test "transcriptvb-edge-cases" "VB genome index or GSE110004 reads not available"
  fi

  local vb_txome="${VB_TRANSCRIPTOME:-/mnt/pikachu/test-datasets-rnaseq/reference/transcriptome.fasta}"
  if [[ -d "${vb_genome}" ]] && [[ -f "${vb_reads}" ]] && [[ -f "${vb_txome}" ]] && command -v salmon >/dev/null 2>&1; then
    run_test "salmon-parity" \
      bash "${SCRIPT_DIR}/transcriptvb/salmon_parity_test.sh" \
        "${vb_genome}" "${vb_txome}"
  else
    skip_test "salmon-parity" "VB genome, transcriptome, reads, or salmon not available"
  fi
}

# ── Execute selected tiers ──────────────────────────────────────────

case "${TIER}" in
  a)   run_tier_a ;;
  b)   run_tier_b ;;
  all) run_tier_a; run_tier_b ;;
  *)   die "Unknown tier: ${TIER}. Use a, b, or all." ;;
esac

print_summary

if [[ ${FAILED} -gt 0 ]]; then
  exit 1
fi
