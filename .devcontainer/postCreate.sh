#!/usr/bin/env bash
set -euo pipefail

cd /workspaces/STAR-suite
mkdir -p .codespaces-demo/{cache,data,indices,runs}

printf '\nCodespaces demo helpers are ready. Suggested starting points:\n'
printf '  1. bash scripts/codespaces/run_bulk_public_demo.sh --dry-run --genome-dir /path/to/index\n'
printf '  2. bash scripts/codespaces/run_perturb_public_demo.sh --dry-run\n'
printf '  3. bash scripts/codespaces/run_flex_public_demo.sh --dry-run\n'
printf '  4. bash scripts/codespaces/run_slam_public_demo.sh --dry-run --genome-dir /path/to/index\n\n'
