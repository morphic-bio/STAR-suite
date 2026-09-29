#!/usr/bin/env bash
# Builds and runs the per-mate FASTX reader harness
# (core/legacy/source/input/fastx_mate_reader_harness.cpp): the reader threads'
# chunk fill against a copy of STAR's single-threaded fill, plus the intended
# behaviour for mates with different read counts. Needs no genome or data.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$ROOT/core/legacy/source"
OUT="${OUT_ROOT:-$(mktemp -d "${TMPDIR:-/tmp}/fastx-mate-reader-harness.XXXXXX")}"
mkdir -p "$OUT"

make -C "$SRC" fastx-mate-reader-harness > "$OUT/build.log" 2>&1 || {
    echo "ERROR: harness build failed; see $OUT/build.log" >&2
    exit 1
}
"$SRC/fastx_mate_reader_harness" "$OUT" | tee "$OUT/harness.log"
