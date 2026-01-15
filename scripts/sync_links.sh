#!/usr/bin/env bash
set -euo pipefail

MANIFEST=""
ROOT=""

usage() {
  cat <<'USAGE'
Usage: sync_links.sh [--manifest PATH] [--root PATH]

Creates symlinks listed in the manifest. Lines are:
  link_path<TAB>target_path
Comments (#) and blank lines are ignored.
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

while IFS=$'\t' read -r link_path target_path; do
  [[ -z "$link_path" ]] && continue
  [[ "$link_path" =~ ^# ]] && continue

  link_abs="$ROOT/$link_path"
  target_abs="$ROOT/$target_path"

  link_dir="$(dirname "$link_abs")"
  mkdir -p "$link_dir"

  if [[ -L "$link_abs" || -e "$link_abs" ]]; then
    rm -f "$link_abs"
  fi

  rel_target=$(python3 - <<PY
import os
link_dir = os.path.abspath("$link_dir")
target = os.path.abspath("$target_abs")
print(os.path.relpath(target, link_dir))
PY
)

  ln -s "$rel_target" "$link_abs"

  if [[ ! -e "$target_abs" ]]; then
    echo "WARNING: target does not exist yet: $target_abs" >&2
  fi
done < "$MANIFEST"
