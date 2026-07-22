#!/usr/bin/env bash
# Run canonical STAR (Velocyto) once and copy parity-relevant artifacts into BASELINE_DIR.
# Use the STAR binary from the commit you want frozen as the independent reference, then point
#   export UCSF_VELOCYTO_BASELINE_OUTDIR="<baseline-dir>/star_run"
# at tests/run_ucsf_velocyto_exact_*.sh so refactored stream_t1 is diffed against this tree.
#
# Usage:
#   save_velocyto_baseline.sh --profile 100k|2m --threads N --baseline-dir /path/to/baseline
#
# Requires the native STAR Velocyto outs writer so UCSF_VELOCYTO_BASELINE_OUTDIR is --mode-all safe.
#
# Requires: scripts/run_star_velocyto_canonical.sh and a writable baseline directory.
# star_run must be fresh (no Solo.out) unless UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1 — see run_star_velocyto_canonical.sh.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CANONICAL="${SCRIPT_DIR}/run_star_velocyto_canonical.sh"
# Same default as run_star_velocyto_canonical.sh; recorded in BASELINE_MANIFEST for audit.
STAR_BIN="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
export STAR_BIN

PROFILE=""
THREADS=""
BASELINE_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile)
      PROFILE="${2:?}"
      shift 2
      ;;
    --threads)
      THREADS="${2:?}"
      shift 2
      ;;
    --baseline-dir)
      BASELINE_DIR="${2:?}"
      shift 2
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 2
      ;;
  esac
done

die() {
  echo "ERROR: $*" >&2
  exit 1
}

[[ -f "${CANONICAL}" ]] || die "Missing ${CANONICAL}"
[[ -n "${PROFILE}" ]] || die "--profile required"
[[ -n "${THREADS}" ]] || die "--threads required"
[[ -n "${BASELINE_DIR}" ]] || die "--baseline-dir required"

RUN_DIR="${BASELINE_DIR}/star_run"
RUN_DIR="${RUN_DIR%/}"
mkdir -p "${BASELINE_DIR}"

if [[ "${UCSF_VELOCYTO_REUSE_STAR_OUTDIR:-0}" != "1" ]]; then
  if [[ -e "${RUN_DIR}/Solo.out" || -e "${RUN_DIR}/Log.out" || -e "${RUN_DIR}/Aligned.out.bam" || -e "${RUN_DIR}/Aligned.out.sam" ]]; then
    die "Refusing non-fresh star_run ${RUN_DIR}. Remove it, use a new --baseline-dir, or UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1"
  fi
  if [[ -d "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix" || -d "${RUN_DIR}/outs/filtered_velocyto_feature_bc_matrix" ]]; then
    die "Refusing non-fresh star_run (outs/*velocyto*). Remove or use UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1"
  fi
fi

"${CANONICAL}" --profile "${PROFILE}" --threads "${THREADS}" --out-prefix "${RUN_DIR}" --prepare-mex

[[ -f "${RUN_DIR}/outs/velocyto_feature_bc_matrix_manifest.json" ]] || die "Missing native Velocyto manifest"
[[ -f "${RUN_DIR}/outs/raw_velocyto_feature_bc_matrix/matrix.mtx.gz" ]] || die "Missing native raw Velocyto MEX"
[[ -f "${RUN_DIR}/outs/filtered_velocyto_feature_bc_matrix/matrix.mtx.gz" ]] || die "Missing native filtered Velocyto MEX"

STAMP="${BASELINE_DIR}/BASELINE_MANIFEST.txt"
STAR_ABS="$(readlink -f "${STAR_BIN}" 2>/dev/null || echo "${STAR_BIN}")"
STAR_SHA256="$(sha256sum "${STAR_BIN}" 2>/dev/null | awk '{print $1}')"
STAR_VER_LINE="$("${STAR_BIN}" --version 2>&1 | head -n1 || true)"
GIT_HEAD="$(git -C "${REPO_ROOT}" rev-parse HEAD 2>/dev/null || echo "unknown")"
GIT_DESCRIBE="$(git -C "${REPO_ROOT}" describe --always --dirty 2>/dev/null || echo "unknown")"
STAR_STAT="$(stat -c '%y %s' "${STAR_BIN}" 2>/dev/null || stat -f '%Sm %z' "${STAR_BIN}" 2>/dev/null || echo "unknown")"
{
  echo "date=$(date -Iseconds)"
  echo "profile=${PROFILE}"
  echo "threads=${THREADS}"
  echo "STAR_VELOCYTO_DETERMINISTIC_REPLAY=${STAR_VELOCYTO_DETERMINISTIC_REPLAY:-}"
  echo "run_dir=${RUN_DIR}"
  echo "STAR_BIN=${STAR_BIN}"
  echo "STAR_BIN_CANONICAL=${STAR_ABS}"
  echo "STAR_BIN_SHA256=${STAR_SHA256}"
  echo "STAR_VERSION_LINE=${STAR_VER_LINE}"
  echo "STAR_BIN_STAT=${STAR_STAT}"
  echo "repo_root=${REPO_ROOT}"
  echo "git_commit=${GIT_HEAD}"
  echo "git_describe=${GIT_DESCRIBE}"
} > "${STAMP}"

mkdir -p "${BASELINE_DIR}/solo_velocyto_snapshot"
cp -a "${RUN_DIR}/Solo.out/Velocyto" "${BASELINE_DIR}/solo_velocyto_snapshot/" || true
mkdir -p "${BASELINE_DIR}/logs"
for f in Log.out Log.final.out Log.progress.out; do
  [[ -f "${RUN_DIR}/${f}" ]] && cp -a "${RUN_DIR}/${f}" "${BASELINE_DIR}/logs/" || true
done

echo "Baseline saved under ${BASELINE_DIR}"
echo "Manifest: ${STAMP}"
echo "Exact-parity harness: export UCSF_VELOCYTO_BASELINE_OUTDIR=${RUN_DIR}"
