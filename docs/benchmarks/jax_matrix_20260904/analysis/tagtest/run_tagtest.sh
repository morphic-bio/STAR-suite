#!/bin/bash
# Small-sample test of the sample-tag tolerance: five known tag classes, each run
# alone, fix on and fix off. Reads Log.final.out tallies. Then the full JAX validation.
set -u; W=/mnt/pikachu/STAR-suite-fix-sampletag-20260905; T=/mnt/pikachu/nvme_jax_v184; D="$T/tagtest"; REF=/home/lhhung/jax_stage_20260903/ref
STAR="$W/core/legacy/source/STAR"; say(){ echo "[$(date -u +%FT%TZ)] $*"; }
until [ -f "$W/build.done" ]; do sleep 15; done; say "initial build finished: $(tail -1 $W/build_star.log)"
( cd "$W/core/legacy/source" && make -j32 STAR > "$W/build_star2.log" 2>&1 ); rc=$?; say "incremental rebuild exit=$rc"; grep -nE "error:" "$W/build_star2.log" | head -5
[ $rc -eq 0 ] || { echo "REBUILD FAILED" > "$D/tagtest.done"; exit 1; }
say "binary $($STAR --version | head -1) sha=$(sha256sum $STAR | cut -c1-16)"
python3 - "$D" <<'PY'
import gzip,sys,collections; D=sys.argv[1]
r1=gzip.open(f"{D}/test_R1.fastq.gz","rt"); r2=gzip.open(f"{D}/test_R2.fastq.gz","rt"); outs={}
for (h1,s1,p1,q1),(h2,s2,p2,q2) in zip(zip(r1,r1,r1,r1),zip(r2,r2,r2,r2)):
    c=h1[1:].split(":")[0]
    if c not in outs: outs[c]=(gzip.open(f"{D}/{c}_R1.fastq.gz","wt"),gzip.open(f"{D}/{c}_R2.fastq.gz","wt"))
    outs[c][0].write(h1+s1+p1+q1); outs[c][1].write(h2+s2+p2+q2)
for a,b in outs.values(): a.close(); b.close()
print("split classes:", list(outs))
PY
common=(--runThreadN 4 --genomeDir "$REF/star_index" --soloType CB_UMI_Simple --soloCBstart 1 --soloUMIstart 17 --soloCBlen 16 --soloUMIlen 12 --soloBarcodeReadLength 0
  --soloCBwhitelist "$REF/737K-fixed-rna-profiling.txt" --flex yes --soloFlexExpectedCellsPerTag 3000 --soloSampleWhitelist "$REF/sample_whitelist_full_16.tsv"
  --soloProbeList "$REF/star_index/flex_probe_artifacts/probe_list.txt" --soloSampleProbes "$REF/probe-barcodes-fixed-rna-profiling-rna.txt" --soloSampleProbeOffset 68
  --soloFlexAllowedTags "$REF/sample_whitelist_full_16.tsv" --outSJtype None --outSAMtype None --outSAMattributes None --soloFeatures Gene --soloCellFilter None
  --soloMultiMappers Rescue --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloStrand Unstranded
  --clipAdapterType CellRanger4 --alignEndsType Local --chimSegmentMin 0 --soloKeysCompat cr --soloHashScreenFile "$REF/h01_cache.bin" --soloBucketMode ram --soloBucketCount 256
  --flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0 --flexNoAlign 1 --readFilesCommand zcat)
printf "%-12s %-8s %8s %8s %8s %8s   %s\n" class mode KEEP DENY UNMATCH PASS policy-line > "$D/tally.txt"
for mode in on off; do
  extra=(); [ $mode = off ] && extra=(--soloSampleTagMismatch 0 --soloSampleSearchNearby no)
  for c in exact variant hd1_unique hd1_ambig shifted; do
    o="$D/run_${c}_$mode"; rm -rf "$o"; mkdir -p "$o"
    "$STAR" "${common[@]}" "${extra[@]}" --readFilesIn "$D/${c}_R2.fastq.gz" "$D/${c}_R1.fastq.gz" --soloFlexOutputPrefix "$o/per_sample" --outTmpDir "$o/_STARtmp" --outFileNamePrefix "$o/" > "$o/stdout" 2>&1
    k=$(grep -aoP 'Hash screen: KEEP \|\t\K\d+' "$o/Log.final.out"); d=$(grep -aoP 'Hash screen: DENY \|\t\K\d+' "$o/Log.final.out"); u=$(grep -aoP 'unmatched sample tag \|\t\K\d+' "$o/Log.final.out"); p=$(grep -aoP 'Hash screen: PASS \|\t\K\d+' "$o/Log.final.out")
    pol=$(grep -a "Flex sample tag detection" "$o/Log.out" | sed 's/.*detection: //')
    printf "%-12s %-8s %8s %8s %8s %8s   %s\n" "$c" "$mode" "${k:-?}" "${d:-?}" "${u:-?}" "${p:-?}" "$pol" >> "$D/tally.txt"
  done
done
cat "$D/tally.txt"
# pass criteria: on -> unmatched ~0 for exact/variant/hd1_unique/shifted and ~2000 for hd1_ambig; off -> unmatched ~2000 for hd1_unique/hd1_ambig/shifted, ~0 for exact/variant
python3 - "$D/tally.txt" <<'PY'
import sys; rows={(r.split()[0],r.split()[1]):int(r.split()[4]) for r in open(sys.argv[1]).read().splitlines()[1:] if r.split()[4].isdigit()}
exp={("exact","on"):0,("variant","on"):0,("hd1_unique","on"):0,("shifted","on"):0,("hd1_ambig","on"):2000,("exact","off"):0,("variant","off"):0,("hd1_unique","off"):2000,("hd1_ambig","off"):2000,("shifted","off"):2000}
bad=[(k,rows.get(k),v) for k,v in exp.items() if rows.get(k) is None or abs(rows[k]-v)>20]
print("SMALL-SAMPLE TEST:", "PASS" if not bad else f"FAIL {bad}")
open(sys.argv[1]+".verdict","w").write("PASS" if not bad else "FAIL")
PY
if grep -q PASS "$D/tally.txt.verdict"; then say "small test passed; running full JAX validation"; "$T/run_tagfix_test.sh" > "$T/tagfix_driver.log" 2>&1; cat "$T/tagfix_driver.log"; else say "small test FAILED; full validation not run"; fi
echo "finished=$(date -u +%FT%TZ)" > "$D/tagtest.done"
