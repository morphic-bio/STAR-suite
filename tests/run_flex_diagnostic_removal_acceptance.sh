#!/usr/bin/env bash
# Serial, source-tree acceptance gates for withdrawal of the optional diagnostic.
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
export STAR_BIN="${STAR_BIN:-$root/core/legacy/source/STAR}"
out="${OUT_ROOT:-$(mktemp -d /tmp/star-flex-removal.XXXXXX)}"
mkdir -p "$out"
[[ ! -e "$out/summary.tsv" ]] || { echo "Use a fresh OUT_ROOT" >&2; exit 2; }
printf 'test\tstatus\tlog\n' > "$out/summary.tsv"
cd "$root"
run() {
    local name="$1"
    shift
    echo "Testing $name; log: $out/$name.log"
    if "$@" > "$out/$name.log" 2>&1; then
        printf '%s\tPASS\t%s\n' "$name" "$out/$name.log" >> "$out/summary.tsv"
    else
        printf '%s\tFAIL\t%s\n' "$name" "$out/$name.log" >> "$out/summary.tsv"
        tail -n 50 "$out/$name.log"
        return 1
    fi
}
run removed python3 tests/test_flex_gdna_removed.py
run packed_counts make -C core/legacy/source test_FlexProbeRegion
run cache_formats bash tests/test_flex_khash_storage.sh "$out/cache"
run htslib_discovery python3 tests/test_htslib_build_discovery.py
run scrna_exact python3 tests/test_scrna_gex_counts.py --outdir "$out/scrna"
run solo_smoke bash tests/run_solo_smoke.sh
run cr_golden bash tests/run_scrna_sidecar_off_golden.sh
run flex_tiny env WORKDIR="$out/flex_tiny" bash tests/run_flex_tiny_public_smoke.sh
echo "PASS: diagnostic-removal acceptance; $out/summary.tsv"
