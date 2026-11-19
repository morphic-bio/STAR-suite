#!/usr/bin/env bash
set -euo pipefail
exec python3 "$(dirname "$0")/run_pipseq_ipm_regression.py" "$@"
