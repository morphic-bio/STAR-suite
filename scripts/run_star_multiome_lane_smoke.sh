#!/usr/bin/env bash
# Compatibility launcher for the canonical Morphic recipe copy.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export STAR_SUITE_ROOT="${STAR_SUITE_ROOT:-$(cd -- "${SCRIPT_DIR}/.." && pwd)}"
SCRIPT_NAME="$(basename -- "$0")"
RECIPE_ROOT="${MORPHIC_RECIPES_ROOT:-/mnt/pikachu/morphic-recipes}"
TARGET="${RECIPE_ROOT}/scripts/${SCRIPT_NAME}"

if [[ ! -f "${TARGET}" ]]; then
  {
    printf 'ERROR: canonical Morphic recipe not found: %s\n' "${TARGET}"
    printf 'Set MORPHIC_RECIPES_ROOT=/path/to/morphic-recipes or run the recipe there directly.\n'
  } >&2
  exit 1
fi

if [[ -x "${TARGET}" ]]; then
  exec "${TARGET}" "$@"
fi

exec bash "${TARGET}" "$@"
