#!/bin/bash
# Edge-position DENY fix: regenerate the H0/H1 cache with the patched generator (terminal variants inherit
# the probe's gene), verify the DENY-by-position profile, then rerun fix-on JAX CBQ with the new cache.
set -u; export PATH="$HOME/.cargo/bin:$HOME/.local/bin:$PATH"
T=/mnt/pikachu/nvme_jax_v184; W=/mnt/pikachu/STAR-suite-fix-sampletag-20260905; STAR=$W/core/legacy/source/STAR
REF=/home/lhhung/jax_stage_20260903/ref; CBQ=/home/lhhung/jax_stage_20260903/cbq; OUT=/storage/jax_tagfix; NEWCACHE=$REF/h01_cache_edgefix.bin
CONC=/mnt/pikachu/STAR-suite/docs/benchmarks/jax_matrix_20260904/concordance_vs_cr.py; SENT=$T/edgefix.done; rm -f "$SENT"
say(){ echo "[$(date -u +%FT%TZ)] $*"; }
manifest(){ (cd "$1" && find . -path '*/Gene/filtered/*' -type f ! -name gdna_metrics.json | sort | xargs sha256sum | sha256sum | cut -c1-16); }
say "binary sha=$(sha256sum $STAR | cut -c1-16)  generating $NEWCACHE"
G=$OUT/cachegen_edgefix; rm -rf $G; mkdir -p $G
/usr/bin/time -v -o $T/cachegen_edgefix.time.txt "$STAR" --runMode hashCacheGenerate --runThreadN 30 --genomeDir $REF/star_index \
  --soloType CB_UMI_Simple --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 --soloBarcodeReadLength 0 \
  --soloCBwhitelist $REF/737K-fixed-rna-profiling.txt --flex yes --soloProbeList $REF/star_index/probe_gene_list.txt --soloFlexExpectedCellsPerTag 3000 \
  --soloSampleWhitelist $REF/sample_whitelist_full_16.tsv --soloSampleProbes $REF/probe-barcodes-fixed-rna-profiling-rna.txt --soloSampleProbeOffset 68 \
  --soloFlexAllowedTags $REF/sample_whitelist_full_16.tsv --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR \
  --soloCellFilter None --soloFeatures Gene --alignEndsType Local --soloStrand Unstranded --soloKeysCompat cr --outSAMtype None \
  --hashCacheTiers "H0,H1" --hashCacheOutput $NEWCACHE --outFileNamePrefix $G/ --outTmpDir $G/_STARtmp > $T/cachegen_edgefix.stdout 2>&1
say "generation exit=$? wall=$(grep -aoP 'Elapsed.*: \K.*' $T/cachegen_edgefix.time.txt) size=$(stat -c %s $NEWCACHE 2>/dev/null)"; grep -a "HASH-CACHE-GEN" $G/Log.out | tail -3
python3 - $NEWCACHE <<'PY'
import sys, numpy as np, csv
from collections import Counter
hdr=np.dtype([("magic","S8"),("version","<u2"),("k","<u2"),("recsize","<u4"),("n","<u8")]); rec=np.dtype([("lo","<u8"),("hi","<u8"),("gene","<u4"),("cls","u1"),("neg","u1"),("res","<u2")])
R=np.fromfile(sys.argv[1],dtype=rec,offset=hdr.itemsize); print(f"new cache records {len(R):,}; class counts {dict(Counter(R['cls'].tolist()))}; DENY {int((R['neg']!=0).sum()):,}")
keys=(R["hi"].astype(object)<<64)|R["lo"].astype(object); deny=set(k for k,ng in zip(keys.tolist(),R["neg"].tolist()) if ng); present=set(keys.tolist())
tbl=str.maketrans("ACGT","0123"); pos=Counter(); absent=0
for row in csv.reader(l for l in open("/home/lhhung/jax_stage_20260903/ref/star_index/flex_probe_artifacts/filtered_probe_set.csv") if not l.startswith("#")):
    if row[0]=="gene_id" or row[3]!="TRUE": continue
    s=row[1]
    for i in (0,1,2,25,47,48,49):
        for b in "ACGT":
            if b==s[i]: continue
            k=int((s[:i]+b+s[i+1:]).translate(tbl),4)
            if k in deny: pos[i]+=1
            elif k not in present: absent+=1
print("DENY variants by position (0,1,2,25,47,48,49):", dict(pos), " absent:", absent, "  [old cache: 0/1/48/49 ~148-155k each, interior ~6.5k]")
zic="ATGTTCATATTCATGGGGCCGTACTGGTTGTGGAGTTGCGCCGCCGAGTA"; k=int(("G"+zic[1:]).translate(tbl),4); print("ZIC2 f249247 A0G verdict now:", "DENY" if k in deny else ("KEEP" if k in present else "absent"))
PY
common=(--runThreadN 32 --dynamicThreadInterface 1 --genomeDir "$REF/star_index" --soloType CB_UMI_Simple --soloCBstart 1 --soloUMIstart 17 --soloCBlen 16 --soloUMIlen 12
  --soloBarcodeReadLength 0 --soloCBwhitelist "$REF/737K-fixed-rna-profiling.txt" --flex yes --soloFlexExpectedCellsPerTag 3000
  --soloSampleWhitelist "$REF/sample_whitelist_full_16.tsv" --soloProbeList "$REF/star_index/flex_probe_artifacts/probe_list.txt"
  --soloSampleProbes "$REF/probe-barcodes-fixed-rna-profiling-rna.txt" --soloSampleProbeOffset 68 --soloFlexAllowedTags "$REF/sample_whitelist_full_16.tsv"
  --limitIObufferSize 50000000 50000000 --outSJtype None --outSAMtype None --outSAMattributes None --soloFeatures Gene --soloCellFilter None
  --soloMultiMappers Rescue --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloStrand Unstranded
  --clipAdapterType CellRanger4 --alignEndsType Local --chimSegmentMin 0 --soloKeysCompat cr --soloHashScreenFile "$NEWCACHE"
  --soloBucketMode ram --soloBucketCount 256 --flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0 --flexNoAlign 1 --crAssignConsumerThreads -1 --crAssignSearchThreads 1
  --readFilesType Binseq PE --readFilesCbqRangeMode range --readFilesIn "$(ls $CBQ/lane_00?.cbq | paste -sd,)")
o=$OUT/E_edgefix; rm -rf "$o"; mkdir -p "$o"; sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' || true; say "E_edgefix starting"
/usr/bin/time -v -o "$T/tagfix_E_edgefix.time.txt" "$STAR" "${common[@]}" --soloFlexOutputPrefix "$o/per_sample" --outTmpDir "$o/_STARtmp" --outFileNamePrefix "$o/" > "$T/tagfix_E_edgefix.stdout" 2>&1
say "E_edgefix exit=$? wall=$(grep -aoP 'Elapsed.*: \K.*' $T/tagfix_E_edgefix.time.txt) out=$(manifest $o/per_sample)  (A_fix_on was 23cd7d3c0cc9a00f)"
grep -aE "Hash screen: (KEEP|DENY|unmatched sample tag) %" "$o/Log.final.out"
python3 "$CONC" "$o/per_sample" star-edgefix > "$T/concordance_edgefix.txt" 2>&1; grep -vE '^#' "$T/concordance_edgefix.txt"
python3 - "$o/per_sample" <<'PY'
import sys,scipy.io,numpy as np,gzip; sys.path.insert(0,"/mnt/pikachu/STAR-suite/docs/benchmarks/jax_matrix_20260904"); import concordance_vs_cr as cc
for bc,s in cc.DEFAULT_GROUPS:
    cx,cb,cg=cc.read_mex(cc.DEFAULT_CRROOT/bc/"count"/"sample_filtered_feature_bc_matrix"); qx,qb,qg=cc.read_mex(f"{sys.argv[1]}/{s[0]}/Gene/filtered"); ax,ab,ag=cc.read_mex(f"/storage/jax_tagfix/A_fix_on/{'per_sample'}/{s[0]}/Gene/filtered")
    print(f"  {bc:18s} UMIs ours/CR = {qx.sum()/cx.sum():.4f}   (A_fix_on {ax.sum()/cx.sum():.4f}; pre-fix 0.9737/0.9757/0.9764/0.9731)")
    gi={g:i for i,g in enumerate(cg)}; qi={g:i for i,g in enumerate(qg)}; ai={g:i for i,g in enumerate(ag)}
    for gene,gid in (("ZIC2","ENSG00000043355"),("NNAT","ENSG00000053438"),("APOE","ENSG00000130203"),("N4BP2","ENSG00000078177")):
        print(f"      {gene:6s} CR {int(cx[:,gi[gid]].sum()):>8,}  A_fix_on {int(ax[:,ai[gid]].sum()):>8,}  edgefix {int(qx[:,qi[gid]].sum()):>8,}")
PY
echo "finished=$(date -u +%FT%TZ)" > "$SENT"
