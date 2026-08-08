#!/usr/bin/env bash
set -euo pipefail
# === MINIMAL 10x multiome recipe (the compose-up floor) =====================
# The smallest functional multiome pipeline. It emits the analysis-ready CORE
# only and nothing else:
#   - GeneFull gene x cell MEX (raw + filtered)
#   - ATAC fragments (possorted BAM) + binary sidecar
#   - MACS3 narrow peaks + peak x cell MEX
# This is apples-to-apples with Cell Ranger ARC --no-bam followed by a
# Signac/MACS peak re-call — the common "matrices + re-called peaks" workflow.
#
# COMPOSE UP from here. If your target ALSO needs RNA-velocity, a GEX BAM, or the
# remote CellBender/MuData downstream, do NOT fork this script: run the full
# recipe (run_star_multiome_lane_smoke.sh) and ADD only the layer you need — see
# its COMPOSITION block and this catalog's compose-up policy.
# Parameter VALUES for any added layer come from the provenance oracle
# (public or project provenance records), not invented.
#
# This wrapper is intentionally thin — it just sets --profile matrices-peaks on
# the one engine — so the minimal floor and the full recipe can never drift.
# Preview first:  run_multiome_minimal.sh ... --dry-run
# All flags of run_star_multiome_lane_smoke.sh are accepted and passed through.
# ============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${SCRIPT_DIR}/run_star_multiome_lane_smoke.sh" --profile matrices-peaks "$@"
