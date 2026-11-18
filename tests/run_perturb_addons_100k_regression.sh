#!/usr/bin/env bash
# Pinned 100K-read regression of the split-pool and PIP RNA modes on public
# reads. PERTURB_ADDONS_PROFILES selects "splitseq", "pipseq" or both (default).
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
args=()
if [[ -n "${OUT_ROOT:-}" ]]; then
    args+=(--outdir "$OUT_ROOT")
fi
exec python3 "$here/run_perturb_addons_100k_regression.py" "${args[@]}" "$@"
