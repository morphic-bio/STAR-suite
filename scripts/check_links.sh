#!/usr/bin/env bash
set -euo pipefail

MANIFEST=""
ROOT=""

usage() {
  cat <<'USAGE'
Usage: check_links.sh [--manifest PATH] [--root PATH]

Validates symlinks listed in the manifest:
- link exists
- link is a symlink (not a regular file)
- target exists
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --root) ROOT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
MANIFEST="${MANIFEST:-$ROOT/links/manifest.tsv}"

if [[ ! -f "$MANIFEST" ]]; then
  echo "ERROR: manifest not found: $MANIFEST" >&2
  exit 1
fi

fail=0

while IFS=$'\t' read -r link_path target_path; do
  [[ -z "$link_path" ]] && continue
  [[ "$link_path" =~ ^# ]] && continue

  link_abs="$ROOT/$link_path"
  target_abs="$ROOT/$target_path"

  if [[ ! -L "$link_abs" ]]; then
    if [[ -e "$link_abs" ]]; then
      echo "ERROR: not a symlink: $link_abs" >&2
    else
      echo "ERROR: missing link: $link_abs" >&2
    fi
    fail=1
    continue
  fi
  if [[ ! -e "$target_abs" ]]; then
    echo "ERROR: missing target: $target_abs (from $link_abs)" >&2
    fail=1
  fi
done < "$MANIFEST"

exit $fail
