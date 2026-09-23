#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
caller="${CALL_FEATURES_BIN:-${repo}/core/features/process_features/call_features}"
[[ -x "$caller" ]] || { echo "Build process_features/call_features first: $caller" >&2; exit 1; }

work=$(mktemp -d /tmp/larry_dominant_calling.XXXXXX)
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/mex"
printf '%s\n' LARRY_A LARRY_B > "$work/mex/features.txt"
printf '%s\n' CELL_1 CELL_2 CELL_3 CELL_4 CELL_5 > "$work/mex/barcodes.txt"
printf '%s\n' \
  '%%MatrixMarket matrix coordinate integer general' \
  '%' \
  '2 5 9' \
  '1 1 2' '2 1 1' \
  '1 2 3' '2 2 2' \
  '1 3 2' '2 3 2' \
  '1 4 1' \
  '1 5 4' '2 5 2' > "$work/mex/matrix.mtx"

"$caller" --guide-caller dominant --min_counts 2 --fraction 0 \
  --margin 1 --min-ratio 2 "$work/mex" "$work/calls" > "$work/caller.log"

awk -F, '
  NR == 1 { if ($1 != "barcode" || $2 != "feature_call") exit 1; next }
  $1 == "CELL_1" && $2 == "LARRY_A" { ok1 = 1 }
  $1 == "CELL_2" && $2 == "Multiplet" { ok2 = 1 }
  $1 == "CELL_3" && $2 == "Multiplet" { ok3 = 1 }
  $1 == "CELL_4" && $2 == "Unassigned" { ok4 = 1 }
  $1 == "CELL_5" && $2 == "LARRY_A" { ok5 = 1 }
  END { if (!(ok1 && ok2 && ok3 && ok4 && ok5)) exit 1 }
' "$work/calls/feature_calls.csv"

"$caller" --guide-caller dominant --min_counts 2 --fraction 0 \
  --margin 1 --min-ratio 3 "$work/mex" "$work/strict_calls" > "$work/strict.log"
awk -F, '
  $1 == "CELL_1" && $2 == "Multiplet" { low_runner_up = 1 }
  $1 == "CELL_5" && $2 == "Multiplet" { passing_runner_up = 1 }
  END { if (!(low_runner_up && passing_runner_up)) exit 1 }
' "$work/strict_calls/feature_calls.csv"

"$caller" --guide-caller dominant --min_counts 1 --fraction 0 \
  --margin 1 --min-ratio 1 "$work/mex" "$work/production_calls" > "$work/production.log"
awk -F, '
  $1 == "CELL_1" && $2 == "LARRY_A" { c1 = 1 }
  $1 == "CELL_2" && $2 == "LARRY_A" { c2 = 1 }
  $1 == "CELL_3" && $2 == "Multiplet" { c3 = 1 }
  $1 == "CELL_4" && $2 == "LARRY_A" { c4 = 1 }
  $1 == "CELL_5" && $2 == "LARRY_A" { c5 = 1 }
  END { if (!(c1 && c2 && c3 && c4 && c5)) exit 1 }
' "$work/production_calls/feature_calls.csv"

echo 'PASS: production top-count dominance and optional 2-UMI/2:1 thresholds'
