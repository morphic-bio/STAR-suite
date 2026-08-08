#!/usr/bin/env bash

set -euo pipefail

BUNDLE=""
EXPECTED_LABEL=""
EXPECTED_VERSION="1.7.0"
EXPECTED_COMMIT=""
MANIFEST_OUT=""

usage() {
  cat <<USAGE
Usage: $0 --bundle <path> --expected-label <label> [--expected-version <version>] [--expected-commit <sha>] [--manifest-out <path>]
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
    --expected-commit)
      EXPECTED_COMMIT="$2"
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
  echo "ERROR: expected STAR-suite version $EXPECTED_VERSION, got $version_output" >&2
  exit 1
fi
source_revision="$($prefix/bin/STAR --source-revision)"
if [[ ! "$source_revision" =~ ^[0-9a-f]{40}$ ]]; then
  echo "ERROR: invalid STAR-suite source revision: $source_revision" >&2
  exit 1
fi
variant_metadata="$workdir/unpack/variants/${selected_label}/release-metadata.env"
if [[ ! -f "$variant_metadata" ]]; then
  echo "ERROR: selected installer variant lacks release metadata" >&2
  exit 1
fi
unset COMMIT_SHA
# shellcheck disable=SC1090
source "$variant_metadata"
metadata_commit="${COMMIT_SHA:-}"
if [[ ! "$metadata_commit" =~ ^[0-9a-f]{40}$ \
    || "$source_revision" != "$metadata_commit" ]]; then
  echo "ERROR: selected installer binary source revision does not match metadata" >&2
  exit 1
fi
if [[ -n "$EXPECTED_COMMIT" && "$source_revision" != "${EXPECTED_COMMIT,,}" ]]; then
  echo "ERROR: binary source revision $source_revision does not match ${EXPECTED_COMMIT,,}" >&2
  exit 1
fi
if [[ -n "$EXPECTED_COMMIT" ]]; then
  /usr/local/bin/check_spatial_release_binary.sh \
    --binary "$prefix/bin/STAR" \
    --expected-version "$EXPECTED_VERSION" \
    --expected-commit "$EXPECTED_COMMIT"
fi
for tool in molecule_first_resolver molecule_first_bam_ledger molecule_first_materialize; do
  if [[ ! -x "$prefix/bin/$tool" ]] || [[ "$($prefix/bin/$tool --version)" != "$EXPECTED_VERSION" ]]; then
    echo "ERROR: installed $tool missing or version-mismatched" >&2
    exit 1
  fi
done
for tool in molecule_first_resolver molecule_first_bam_ledger molecule_first_materialize \
            transcriptvb_finalize trim_qc_fastq trim_qc_merge; do
  if [[ ! -x "$prefix/bin/$tool" ]]; then
    echo "ERROR: installed release companion $tool missing" >&2
    exit 1
  fi
  if ! ldd_output="$(ldd "$prefix/bin/$tool" 2>&1)"; then
    if [[ "$ldd_output" == *"not a dynamic executable"* || "$ldd_output" == *"statically linked"* ]]; then
      continue
    fi
    echo "ERROR: could not inspect runtime dependencies for $tool: $ldd_output" >&2
    exit 1
  fi
  if [[ "$ldd_output" == *'=> not found'* ]]; then
    echo "ERROR: unresolved runtime dependency for $tool:" >&2
    printf '%s\n' "$ldd_output" >&2
    exit 1
  fi
done
for data_file in \
  share/star-suite/SNAPSHOTS.json \
  share/star-suite/catalogs/official/catalog.yaml \
  share/star-suite/evidence/official/schema/record-v1.schema.json
do
  if [[ ! -f "$prefix/$data_file" ]]; then
    echo "ERROR: installed official release data missing: $data_file" >&2
    exit 1
  fi
  if ! cmp -s "$workdir/unpack/$data_file" "$prefix/$data_file"; then
    echo "ERROR: installed official release data differs: $data_file" >&2
    exit 1
  fi
done

manifest="Bundle: $(basename "$BUNDLE")
Container glibc: $(detect_glibc)
Selected label: ${selected_label}
Selected baseline: ${selected_baseline}
Selected description: ${selected_description}
STAR-suite version: ${version_output}
STAR-suite source revision: ${source_revision}
"

printf '%s' "$manifest"
if [[ -n "$MANIFEST_OUT" ]]; then
  mkdir -p "$(dirname "$MANIFEST_OUT")"
  printf '%s' "$manifest" > "$MANIFEST_OUT"
fi
