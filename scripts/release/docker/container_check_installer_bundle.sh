#!/usr/bin/env bash

set -euo pipefail

BUNDLE=""
EXPECTED_LABEL=""
EXPECTED_VERSION="2.7.11b"
MANIFEST_OUT=""

usage() {
  cat <<USAGE
Usage: $0 --bundle <path> --expected-label <label> [--expected-version <version>] [--manifest-out <path>]
USAGE
}

detect_glibc() {
  if command -v getconf >/dev/null 2>&1; then
    local out
    out="$(getconf GNU_LIBC_VERSION 2>/dev/null || true)"
    if [[ "$out" =~ ([0-9]+\.[0-9]+) ]]; then
      printf '%s\n' "${BASH_REMATCH[1]}"
      return 0
    fi
  fi

  if command -v ldd >/dev/null 2>&1; then
    local first_line
    first_line="$(ldd --version 2>&1 | head -n1)"
    if [[ "$first_line" =~ ([0-9]+\.[0-9]+) ]]; then
      printf '%s\n' "${BASH_REMATCH[1]}"
      return 0
    fi
  fi

  printf 'unknown\n'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bundle)
      BUNDLE="$2"
      shift 2
      ;;
    --expected-label)
      EXPECTED_LABEL="$2"
      shift 2
      ;;
    --expected-version)
      EXPECTED_VERSION="$2"
      shift 2
      ;;
    --manifest-out)
      MANIFEST_OUT="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$BUNDLE" || -z "$EXPECTED_LABEL" ]]; then
  echo "ERROR: --bundle and --expected-label are required" >&2
  exit 1
fi

if [[ ! -f "$BUNDLE" ]]; then
  echo "ERROR: bundle not found: $BUNDLE" >&2
  exit 1
fi

workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

mkdir -p "$workdir/unpack"
tar -xzf "$BUNDLE" -C "$workdir/unpack"
selection="$($workdir/unpack/install.sh --print-selection)"
IFS=$'\t' read -r selected_label selected_baseline selected_relpath selected_description <<< "$selection"
if [[ "$selected_label" != "$EXPECTED_LABEL" ]]; then
  echo "ERROR: expected installer to select $EXPECTED_LABEL, got $selected_label" >&2
  exit 1
fi

prefix="$workdir/prefix"
"$workdir/unpack/install.sh" --prefix "$prefix" --force
version_output="$($prefix/bin/STAR --version)"
if [[ "$version_output" != "$EXPECTED_VERSION" ]]; then
  echo "ERROR: expected STAR version $EXPECTED_VERSION, got $version_output" >&2
  exit 1
fi

manifest="Bundle: $(basename "$BUNDLE")
Container glibc: $(detect_glibc)
Selected label: ${selected_label}
Selected baseline: ${selected_baseline}
Selected description: ${selected_description}
STAR version: ${version_output}
"

printf '%s' "$manifest"
if [[ -n "$MANIFEST_OUT" ]]; then
  mkdir -p "$(dirname "$MANIFEST_OUT")"
  printf '%s' "$manifest" > "$MANIFEST_OUT"
fi
