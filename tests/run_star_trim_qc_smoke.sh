#!/usr/bin/env bash
# STAR's trim QC must actually record reads, and per-shard dumps must merge
# back to what a single run produces.
#
# This exists because --trimQcReport silently emitted a well-formed report with
# "total_reads": 0 -- the collector was constructed and merged across chunk
# threads, but nothing ever fed it a read. The files were produced, so the
# failure was invisible unless you opened them.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/.." && pwd)"
src="$root/core/legacy/source"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

for tool in STAR trim_qc_merge; do
    [ -x "$src/$tool" ] || make -C "$src" "$tool" WITH_CHROMAP=0 >/dev/null
done

READS=2000
python3 - "$work" "$READS" <<'PY'
import random, sys, os
random.seed(7)
d, n = sys.argv[1], int(sys.argv[2])
ref = ''.join(random.choice("ACGT") for _ in range(200000))
with open(os.path.join(d, "ref.fa"), "w") as f:
    f.write(">chrT\n")
    for i in range(0, len(ref), 60): f.write(ref[i:i+60] + "\n")
comp = {'A':'T','C':'G','G':'C','T':'A'}
r1, r2 = [], []
for i in range(n):
    p = random.randint(0, len(ref) - 400)
    a = ref[p:p+100]
    b = ''.join(comp[c] for c in reversed(ref[p+200:p+300]))
    q = 'I' * 100
    r1.append("@rd%d\n%s\n+\n%s\n" % (i, a, q))
    r2.append("@rd%d\n%s\n+\n%s\n" % (i, b, q))
open(os.path.join(d, "r1.fastq"), "w").write(''.join(r1))
open(os.path.join(d, "r2.fastq"), "w").write(''.join(r2))
PY

mkdir -p "$work/genome"
"$src/STAR" --runMode genomeGenerate --genomeDir "$work/genome" \
    --genomeFastaFiles "$work/ref.fa" --genomeSAindexNbases 8 --runThreadN 4 \
    --outFileNamePrefix "$work/gg_" >/dev/null 2>&1

mkdir -p "$work/full"
"$src/STAR" --runMode alignReads --genomeDir "$work/genome" \
    --readFilesIn "$work/r1.fastq" "$work/r2.fastq" \
    --outSAMtype None --runThreadN 4 \
    --trimQcReport "$work/full/qc" --outFileNamePrefix "$work/full/" >/dev/null 2>&1

python3 - "$work/full/qc.trim_qc.json" "$READS" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); n = int(sys.argv[2])
if d["total_reads"] != 2 * n:
    print("FAIL: total_reads %s, expected %s" % (d["total_reads"], 2 * n)); sys.exit(1)
for m in d["mates"]:
    if m["reads"] != n:
        print("FAIL: mate %s has %s reads, expected %s" % (m["mate"], m["reads"], n)); sys.exit(1)
    if m["length_hist"][100] != n:
        print("FAIL: mate %s length histogram missing 100bp reads" % m["mate"]); sys.exit(1)
    p0 = m["positions"][0]
    if abs(p0["mean_qual"] - 40.0) > 1e-9:
        print("FAIL: mean quality %s, expected 40.0" % p0["mean_qual"]); sys.exit(1)
    if p0["base_counts"]["N"] != 0 or sum(p0["base_counts"].values()) != n:
        print("FAIL: base counts at position 1 do not sum to %s" % n); sys.exit(1)
print("PASS: STAR recorded %s reads with correct lengths, quality and composition" % d["total_reads"])
PY

# Split, run per shard, merge, and require the same histograms.
head -$(( READS * 2 )) "$work/r1.fastq" > "$work/a_r1.fastq"
head -$(( READS * 2 )) "$work/r2.fastq" > "$work/a_r2.fastq"
tail -$(( READS * 2 )) "$work/r1.fastq" > "$work/b_r1.fastq"
tail -$(( READS * 2 )) "$work/r2.fastq" > "$work/b_r2.fastq"
for s in a b; do
    mkdir -p "$work/$s"
    "$src/STAR" --runMode alignReads --genomeDir "$work/genome" \
        --readFilesIn "$work/${s}_r1.fastq" "$work/${s}_r2.fastq" \
        --outSAMtype None --runThreadN 2 \
        --trimQcShardOut "$work/$s.trimqc" \
        --outFileNamePrefix "$work/$s/" >/dev/null 2>&1
done
"$src/trim_qc_merge" --out-prefix "$work/merged" --stage raw \
    "$work/a.trimqc" "$work/b.trimqc" >/dev/null

python3 - "$work/full/qc.trim_qc.json" "$work/merged.trim_qc.json" <<'PY'
import json, sys
full = json.load(open(sys.argv[1])); mrg = json.load(open(sys.argv[2]))
if full["total_reads"] != mrg["total_reads"]:
    print("FAIL: merged total_reads %s vs full %s" % (mrg["total_reads"], full["total_reads"])); sys.exit(1)
for a, b in zip(full["mates"], mrg["mates"]):
    for field in ("reads", "length_hist", "gc_hist", "positions"):
        if a[field] != b[field]:
            print("FAIL: mate %s differs in %s" % (a["mate"], field)); sys.exit(1)
print("PASS: two STAR shards merge to the single-run histograms exactly")
PY
