#!/usr/bin/env bash
set -euo pipefail

WORKSPACE_ROOT="/workspaces/STAR-suite"
if [[ ! -d "${WORKSPACE_ROOT}" ]]; then
  WORKSPACE_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
fi

cd "${WORKSPACE_ROOT}"
mkdir -p .codespaces-demo/{cache,data,indices,runs}

EXT_DIR="${PWD}/tools/vscode-star-suite-walkthrough"
VSIX_DIR="${PWD}/.codespaces-demo/cache/extensions"
VSIX_PATH="${VSIX_DIR}/star-suite-walkthrough.vsix"
mkdir -p "${VSIX_DIR}"

if [[ "${SKIP_WALKTHROUGH_INSTALL:-0}" != "1" ]] && command -v code >/dev/null 2>&1 && command -v npx >/dev/null 2>&1 && [[ -d "${EXT_DIR}" ]]; then
  (
    cd "${EXT_DIR}"
    npx --yes @vscode/vsce package --out "${VSIX_PATH}" >/dev/null
  )
  timeout 60s code --install-extension "${VSIX_PATH}" --force >/dev/null 2>&1 || true
fi

printf '\nCodespaces demo helpers are ready. Suggested starting points:\n'
printf '  1. Command Palette -> STAR-suite: Open Codespaces Walkthrough\n'
printf '  2. Read docs/codespaces/00_overview.md\n'
printf '  3. Optional foundation: bash scripts/codespaces/build_public_chr22y_index.sh --threads 4\n'
printf '  4. Bulk chapter: bash scripts/codespaces/run_bulk_public_demo.sh --dry-run --genome-dir /path/to/index\n'
printf '  5. SLAM chapter: bash scripts/codespaces/run_slam_public_demo.sh --dry-run --genome-dir /path/to/index\n\n'
