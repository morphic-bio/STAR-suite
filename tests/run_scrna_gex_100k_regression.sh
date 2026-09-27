#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
args=()
if [[ -n "${OUT_ROOT:-}" ]]; then
    args+=(--outdir "$OUT_ROOT")
fi
exec python3 "$here/run_scrna_gex_100k_regression.py" "${args[@]}" "$@"
