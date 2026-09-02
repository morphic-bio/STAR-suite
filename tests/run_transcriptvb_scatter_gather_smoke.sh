#!/usr/bin/env bash
# Scatter-gather quantification must equal in-process quantification.
#
# STAR writes one evidence sidecar per shard with --quantVBSidecarOnly, and
# transcriptvb_finalize gathers them into quant.sf without loading the genome
# index. That is what lets quantification run as an HPC scatter-gather job, so
# the round trip is checked against a single-process run over the same reads.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/.." && pwd)"
src="$root/core/legacy/source"
STAR_BIN="${STAR_BIN:-$src/STAR}"
TRANSCRIPTVB_FINALIZE_BIN="${TRANSCRIPTVB_FINALIZE_BIN:-$src/transcriptvb_finalize}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

run_quiet() {
    local stage="$1"
    shift
    local log="$work/${stage}.log"
    local rc
    set +e
    "$@" >"$log" 2>&1
    rc=$?
    set -e
    if [ "$rc" -ne 0 ]; then
        echo "FAIL: $stage exited with status $rc" >&2
        tail -n 80 "$log" >&2
        exit "$rc"
    fi
}

[ -x "$STAR_BIN" ] || make -C "$src" STAR WITH_CHROMAP=0 >/dev/null
[ -x "$TRANSCRIPTVB_FINALIZE_BIN" ] || make -C "$src" transcriptvb-finalize WITH_CHROMAP=0 >/dev/null

python3 - "$work" <<'PY'
import random, sys, os
random.seed(31); d = sys.argv[1]
tx = [''.join(random.choice("ACGT") for _ in range(1200)) for _ in range(40)]
ref = ''.join(tx)
open(os.path.join(d, "ref.fa"), "w").write(">chrT\n" + "\n".join(ref[i:i+60] for i in range(0, len(ref), 60)) + "\n")
rows = []
for t in range(40):
    s = t*1200 + 1
    rows.append('chrT\ttest\texon\t%d\t%d\t.\t+\t.\tgene_id "gene%d"; transcript_id "tx%d";' % (s, s+1199, t//2, t))
open(os.path.join(d, "genes.gtf"), "w").write("\n".join(rows) + "\n")
with open(os.path.join(d, "txome.fa"), "w") as f:
    for t in range(40):
        f.write(">tx%d\n%s\n" % (t, tx[t]))
comp = {'A':'T','C':'G','G':'C','T':'A'}
r1, r2 = [], []
for i in range(4000):
    t = random.randint(0, 39); p = t*1200 + random.randint(0, 900)
    r1.append("@r%d\n%s\n+\n%s\n" % (i, ref[p:p+100], 'I'*100))
    r2.append("@r%d\n%s\n+\n%s\n" % (i, ''.join(comp[c] for c in reversed(ref[p+200:p+300])), 'I'*100))
open(os.path.join(d, "r1.fq"), "w").write(''.join(r1))
open(os.path.join(d, "r2.fq"), "w").write(''.join(r2))
PY

mkdir -p "$work/g"
run_quiet genome_generate "$STAR_BIN" --runMode genomeGenerate --genomeDir "$work/g" --genomeFastaFiles "$work/ref.fa" \
    --sjdbGTFfile "$work/genes.gtf" --genomeSAindexNbases 8 --runThreadN 4 \
    --outFileNamePrefix "$work/gg_"

common=(--genomeDir "$work/g" --outSAMtype None --runThreadN 4 --quantMode TranscriptVB
        --quantVBLibType IU --transcriptomeFasta "$work/txome.fa")

mkdir -p "$work/inproc"
run_quiet in_process_quant "$STAR_BIN" --runMode alignReads "${common[@]}" \
    --readFilesIn "$work/r1.fq" "$work/r2.fq" --outFileNamePrefix "$work/inproc/"

head -8000 "$work/r1.fq" > "$work/a1.fq"; head -8000 "$work/r2.fq" > "$work/a2.fq"
tail -8000 "$work/r1.fq" > "$work/b1.fq"; tail -8000 "$work/r2.fq" > "$work/b2.fq"
for i in 0 1; do
    mkdir -p "$work/s$i"
    if [ "$i" = 0 ]; then R1="$work/a1.fq"; R2="$work/a2.fq"; else R1="$work/b1.fq"; R2="$work/b2.fq"; fi
    run_quiet "sidecar_shard_${i}" "$STAR_BIN" --runMode alignReads "${common[@]}" --readFilesIn "$R1" "$R2" \
        --outFileNamePrefix "$work/s$i/" \
        --quantVBSidecarOut "$work/s$i/evidence.stvb" --quantVBSidecarOnly 1 \
        --quantVBSidecarSampleId smoke --quantVBSidecarInputId smoke-input \
        --quantVBSidecarShardOrdinal "$i" --quantVBSidecarShardCount 2
    [ -s "$work/s$i/evidence.stvb" ] || { echo "FAIL: shard $i wrote no sidecar" >&2; exit 1; }
done

# --no-gc matches the in-process default; the dynamic GC correction is opt-in
# there, and comparing a corrected gather against an uncorrected run is not a
# like-for-like test.
run_quiet transcriptvb_finalize "$TRANSCRIPTVB_FINALIZE_BIN" \
    --genome-dir "$work/g" --transcriptome "$work/txome.fa" \
    --out-prefix "$work/gathered_" --no-gc --threads 4 \
    "$work/s0/evidence.stvb" "$work/s1/evidence.stvb"

python3 - "$work/inproc/quant.sf" "$work/gathered_quant.sf" <<'PY'
import sys
def rd(p):
    rows = open(p).read().splitlines()[1:]
    return {l.split('\t')[0]: tuple(float(x) for x in l.split('\t')[1:]) for l in rows}
a, b = rd(sys.argv[1]), rd(sys.argv[2])
if set(a) != set(b):
    print("FAIL: transcript sets differ"); sys.exit(1)
ks = sorted(a)
for i, lbl, tol in [(3, "NumReads", 1e-6), (2, "TPM", 1e-6), (1, "EffLength", 1e-3)]:
    worst = max(abs(a[k][i] - b[k][i]) / max(1.0, abs(a[k][i])) for k in ks)
    if worst > tol:
        print("FAIL: %s differs by %.3g relative (tol %g)" % (lbl, worst, tol)); sys.exit(1)
tot_a = sum(a[k][3] for k in ks); tot_b = sum(b[k][3] for k in ks)
if abs(tot_a - tot_b) > 1e-6:
    print("FAIL: total NumReads %.6f vs %.6f" % (tot_a, tot_b)); sys.exit(1)
print("PASS: 2 shards gathered == in-process over %d transcripts, %.0f reads" % (len(ks), tot_a))
PY

[ -s "$work/gathered_quant.genes.sf" ] || { echo "FAIL: no gene output" >&2; exit 1; }
echo "PASS: gene-level output produced"
