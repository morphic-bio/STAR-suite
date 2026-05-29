#!/usr/bin/env bash

set -euo pipefail

MODE=""
TARBALL=""
BUNDLE=""
EXPECTED_LABEL=""
EXPECTED_VERSION="1.0.2"
REPO_ROOT=""
PROFILE="core"

usage() {
  cat <<USAGE
Usage:
  $0 --mode tarball --tarball <path> --repo-root <path> [options]
  $0 --mode bundle --bundle <path> --expected-label <label> --repo-root <path> [options]

Options:
  --expected-version VER   STAR-suite version to check (default: 1.0.2)
  --profile PROFILE        Smoke profile to run (default: core)

Profiles:
  core                     Original 3-test smoke (solo, SLAM snp-mask, Flex tiny)
  binary-workflows-tier-a  Full Tier A binary workflow matrix
  all                      Run both core and binary-workflows-tier-a
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --tarball) TARBALL="$2"; shift 2 ;;
    --bundle) BUNDLE="$2"; shift 2 ;;
    --expected-label) EXPECTED_LABEL="$2"; shift 2 ;;
    --expected-version) EXPECTED_VERSION="$2"; shift 2 ;;
    --repo-root) REPO_ROOT="$2"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -d "${REPO_ROOT}" ]] || { echo "ERROR: missing --repo-root" >&2; exit 1; }

workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT
prefix="${workdir}/prefix"
repo_copy="${workdir}/repo"

install_from_tarball() {
  [[ -f "${TARBALL}" ]] || { echo "ERROR: tarball not found: ${TARBALL}" >&2; exit 1; }
  mkdir -p "${workdir}/tarball"
  tar -xzf "${TARBALL}" -C "${workdir}/tarball"
  "${workdir}/tarball/install.sh" --prefix "${prefix}" --force >/dev/null
}

install_from_bundle() {
  [[ -f "${BUNDLE}" ]] || { echo "ERROR: bundle not found: ${BUNDLE}" >&2; exit 1; }
  [[ -n "${EXPECTED_LABEL}" ]] || { echo "ERROR: --expected-label required for bundle mode" >&2; exit 1; }
  mkdir -p "${workdir}/bundle"
  tar -xzf "${BUNDLE}" -C "${workdir}/bundle"
  selection="$(${workdir}/bundle/install.sh --print-selection)"
  IFS=$'\t' read -r selected_label _rest <<< "${selection}"
  [[ "${selected_label}" == "${EXPECTED_LABEL}" ]] || {
    echo "ERROR: expected bundle to select ${EXPECTED_LABEL}, got ${selected_label}" >&2
    exit 1
  }
  "${workdir}/bundle/install.sh" --prefix "${prefix}" --force >/dev/null
}

case "${MODE}" in
  tarball) install_from_tarball ;;
  bundle) install_from_bundle ;;
  *) echo "ERROR: --mode must be tarball or bundle" >&2; exit 1 ;;
esac

STAR_BIN="${prefix}/bin/STAR"
[[ -x "${STAR_BIN}" ]] || { echo "ERROR: installed STAR missing: ${STAR_BIN}" >&2; exit 1; }
[[ "$("${STAR_BIN}" --version)" == "${EXPECTED_VERSION}" ]] || {
  echo "ERROR: installed STAR-suite version mismatch" >&2
  exit 1
}

cp -a "${REPO_ROOT}/." "${repo_copy}/"

run_core() {
  echo "=== Profile: core ==="
  (
    cd "${repo_copy}"
    export STAR_BIN
    bash tests/run_solo_smoke.sh
    bash tests/slam/test_snp_mask_build_smoke.sh
    bash tests/run_flex_tiny_public_smoke.sh
  )
  echo "PASS: core profile"
}

run_binary_workflows_tier_a() {
  echo "=== Profile: binary-workflows-tier-a ==="
  (
    cd "${repo_copy}"
    export STAR_BIN
    bash tests/run_binary_test_matrix.sh --tier a
  )
  echo "PASS: binary-workflows-tier-a profile"
}

case "${PROFILE}" in
  core)
    run_core
    ;;
  binary-workflows-tier-a)
    run_binary_workflows_tier_a
    ;;
  all)
    run_core
    run_binary_workflows_tier_a
    ;;
  *)
    echo "ERROR: unknown profile: ${PROFILE}" >&2
    exit 1
    ;;
esac

echo "PASS: release smoke checks (${MODE}, profile: ${PROFILE})"
