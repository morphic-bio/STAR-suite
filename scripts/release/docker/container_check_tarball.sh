#!/usr/bin/env bash

set -euo pipefail

TARBALL=""
EXPECTED_VERSION="1.0.0"
MANIFEST_OUT=""

usage() {
  cat <<USAGE
Usage: $0 --tarball <path> [--expected-version <version>] [--manifest-out <path>]
USAGE
}

trim() {
  local value="$1"
  value="${value#"${value%%[![:space:]]*}"}"
  value="${value%"${value##*[![:space:]]}"}"
  printf '%s' "$value"
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

runtime_package_for() {
  local path="$1"
  local resolved
  resolved="$(readlink -f "$path" 2>/dev/null || printf '%s' "$path")"
  local candidate
  local query

  if command -v dpkg-query >/dev/null 2>&1; then
    for candidate in \
      "$path" \
      "$resolved" \
      "${path/#\/lib\//\/usr\/lib/}" \
      "${resolved/#\/lib\//\/usr\/lib/}" \
      "${path/#\/usr\/lib\//\/lib/}" \
      "${resolved/#\/usr\/lib\//\/lib/}"; do
      [[ -z "$candidate" ]] && continue
      query="$(dpkg-query -S "$candidate" 2>/dev/null | head -n1 || true)"
      if [[ -n "$query" ]]; then
        if [[ "$query" == diversion\ by\ * ]]; then
          query="${query#diversion by }"
          query="${query%% from*}"
        else
          query="${query%%:*}"
        fi
        printf '%s\n' "$query"
        return 0
      fi
    done
  fi

  printf 'unknown\n'
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tarball)
      TARBALL="$2"
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

if [[ -z "$TARBALL" ]]; then
  echo "ERROR: --tarball is required" >&2
  exit 1
fi

if [[ ! -f "$TARBALL" ]]; then
  echo "ERROR: tarball not found: $TARBALL" >&2
  exit 1
fi

workdir="$(mktemp -d)"
trap 'rm -rf "$workdir"' EXIT

mkdir -p "$workdir/unpack"
tar -xzf "$TARBALL" -C "$workdir/unpack"

installer="$workdir/unpack/install.sh"
metadata="$workdir/unpack/release-metadata.env"
if [[ ! -x "$installer" ]]; then
  echo "ERROR: installer missing from tarball: $installer" >&2
  exit 1
fi

compat_label=""
glibc_baseline=""
asset_name="$(basename "$TARBALL")"
if [[ -f "$metadata" ]]; then
  # shellcheck disable=SC1090
  source "$metadata"
  compat_label="${COMPAT_LABEL:-}"
  glibc_baseline="${GLIBC_BASELINE:-}"
  asset_name="${ASSET_NAME:-$asset_name}"
fi

prefix="$workdir/prefix"
"$installer" --prefix "$prefix" --force >/dev/null
binary="$prefix/bin/STAR"
if [[ ! -x "$binary" ]]; then
  echo "ERROR: installed STAR binary missing: $binary" >&2
  exit 1
fi

version_output="$($binary --version)"
if [[ "$version_output" != "$EXPECTED_VERSION" ]]; then
  echo "ERROR: expected STAR-suite version $EXPECTED_VERSION, got $version_output" >&2
  exit 1
fi

container_glibc="$(detect_glibc)"
max_glibc_symbol="$(grep -aoE 'GLIBC_[0-9]+\.[0-9]+' "$binary" | sort -Vu | tail -n1 || true)"
ldd_output="$({ ldd "$binary"; } 2>&1)"

declare -A unique_packages=()
declare -a lib_rows=()
while IFS= read -r raw_line; do
  line="$(trim "$raw_line")"
  [[ -z "$line" ]] && continue
  [[ "$line" == linux-vdso.so.* ]] && continue
  if [[ "$line" == *'=> not found'* ]]; then
    echo "ERROR: unresolved runtime dependency: $line" >&2
    exit 1
  fi

  soname=""
  path=""
  if [[ "$line" == *'=>'* ]]; then
    soname="${line%% => *}"
    rest="${line#*=> }"
    path="${rest%% *}"
  else
    path="${line%% *}"
    soname="$(basename "$path")"
  fi

  package="$(runtime_package_for "$path")"
  if [[ -n "$package" ]]; then
    unique_packages["$package"]=1
  fi
  lib_rows+=("$soname|$path|${package:-unknown}")
done <<< "$ldd_output"

package_list=""
if [[ "${#unique_packages[@]}" -gt 0 ]]; then
  package_list="$(printf '%s\n' "${!unique_packages[@]}" | sort | paste -sd ',' -)"
fi

manifest="Asset: $asset_name
Compatibility label: ${compat_label:-unknown}
Documented glibc baseline: ${glibc_baseline:-unknown}
Container glibc: ${container_glibc}
Maximum referenced GLIBC symbol: ${max_glibc_symbol:-unknown}
STAR-suite version: ${version_output}
Minimum runtime packages (Ubuntu/Debian package names): ${package_list:-unknown}
Dynamic libraries:
"
for row in "${lib_rows[@]}"; do
  IFS='|' read -r soname path package <<< "$row"
  manifest+="  - ${soname} => ${path}"
  if [[ -n "$package" ]]; then
    manifest+=" [${package}]"
  fi
  manifest+=$'\n'
done

printf '%s' "$manifest"
if [[ -n "$MANIFEST_OUT" ]]; then
  mkdir -p "$(dirname "$MANIFEST_OUT")"
  printf '%s' "$manifest" > "$MANIFEST_OUT"
fi
