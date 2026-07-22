#!/usr/bin/env bash
# BLOCKED: external velocyto.py comparison (plan Phase 4 / runbook “Real external velocyto.py comparison”).
#
# This requires full retained aligned BAMs on the same UCSF surface; header-only stubs are explicitly out of scope.
# Do not treat exit 0 as “passed validation” — it means the check is intentionally not run.
#
# When BAMs exist: implement normalization to shared feature/barcode axes, then diff vs STAR outputs.
# Reference: docs/RUNBOOK_VELOCYTO_BRIDGE_IMPLEMENTATION_20260326.md (Phase 4 section).
set -euo pipefail
cat >&2 <<EOF
BLOCKED: run_ucsf_velocyto_external_compare.sh — external velocyto.py workflow not run (BAM availability / wiring).
This is a distinct validation target from legacy STAR parity; document the blocker in run notes rather than skipping silently.
EOF
exit 0
