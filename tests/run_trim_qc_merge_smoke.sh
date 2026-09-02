#!/usr/bin/env bash
# Merged per-shard trim QC must equal a single pass over the same reads.
#
# The shards deliberately carry different read lengths: merge has to grow its
# per-position arrays, which is the case a naive accumulator gets wrong.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/.." && pwd)"
src="$root/core/legacy/source"
TRIM_QC_FASTQ_BIN="${TRIM_QC_FASTQ_BIN:-$src/trim_qc_fastq}"
TRIM_QC_MERGE_BIN="${TRIM_QC_MERGE_BIN:-$src/trim_qc_merge}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

[ -x "$TRIM_QC_FASTQ_BIN" ] || make -C "$src" trim_qc_fastq WITH_CHROMAP=0 >/dev/null
[ -x "$TRIM_QC_MERGE_BIN" ] || make -C "$src" trim_qc_merge WITH_CHROMAP=0 >/dev/null

python3 - "$work" <<'PY'
import random, os, sys
random.seed(20260807)
d = sys.argv[1]
specs = [(1500, 75), (1200, 150), (900, 100), (1100, 60)]
allr = []
for si, (n, L) in enumerate(specs):
    reads = []
    for i in range(n):
        ln = L if i % 7 else max(20, L - 13)
        seq = ''.join(random.choice("ACGTN" if i % 97 == 0 else "ACGT") for _ in range(ln))
        qual = ''.join(chr(33 + random.randint(2, 40)) for _ in range(ln))
        rec = f"@s{si}_r{i}\n{seq}\n+\n{qual}\n"
        reads.append(rec); allr.append(rec)
    open(os.path.join(d, f"shard{si}.fastq"), "w").write(''.join(reads))
open(os.path.join(d, "all.fastq"), "w").write(''.join(allr))
PY

for i in 0 1 2 3; do
    "$TRIM_QC_FASTQ_BIN" --input "$work/shard$i.fastq" --report "$work/s$i" \
        --shard-out "$work/s$i.trimqc" >/dev/null
done

"$TRIM_QC_MERGE_BIN" --out-prefix "$work/merged" --stage fastq_replay "$work"/s*.trimqc >/dev/null
"$TRIM_QC_FASTQ_BIN" --input "$work/all.fastq" --report "$work/ref" >/dev/null

if ! cmp -s "$work/merged.trim_qc.json" "$work/ref.trim_qc.json"; then
    echo "FAIL: merged shard QC differs from a single pass" >&2
    diff "$work/merged.trim_qc.json" "$work/ref.trim_qc.json" | head -20 >&2
    exit 1
fi

# Merge is a sum, so shard order must not matter.
"$TRIM_QC_MERGE_BIN" --out-prefix "$work/rev" --stage fastq_replay \
    "$work/s3.trimqc" "$work/s1.trimqc" "$work/s0.trimqc" "$work/s2.trimqc" >/dev/null
if ! cmp -s "$work/rev.trim_qc.json" "$work/ref.trim_qc.json"; then
    echo "FAIL: merge is order-dependent" >&2
    exit 1
fi

[ -s "$work/merged.trim_qc.html" ] || { echo "FAIL: no HTML emitted" >&2; exit 1; }

echo "PASS: 4 shards (lengths 75/150/100/60, 4700 reads) merge byte-identically"
echo "PASS: merge is order-independent"
