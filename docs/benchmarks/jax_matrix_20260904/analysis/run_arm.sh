#!/bin/bash
# run_arm.sh ARM [extra STAR args]: fix-on JAX CBQ no-align with the edge-fix cache; report manifest, screen stats,
# concordance, UMI ratios and the hole genes against A_fix_on. Sentinel: $T/arm_$ARM.done
set -u; export PATH="$HOME/.cargo/bin:$HOME/.local/bin:$PATH"; ARM=$1; shift
T=/mnt/pikachu/nvme_jax_v184; W=/mnt/pikachu/STAR-suite-fix-sampletag-20260905; STAR=$W/core/legacy/source/STAR
REF=/home/lhhung/jax_stage_20260903/ref; CBQ=/home/lhhung/jax_stage_20260903/cbq; OUT=/storage/jax_tagfix; CACHE=${CACHE:-$REF/h01_cache_edgefix.bin}
CONC=/mnt/pikachu/STAR-suite/docs/benchmarks/jax_matrix_20260904/concordance_vs_cr.py; SENT=$T/arm_$ARM.done; rm -f "$SENT"
say(){ echo "[$(date -u +%FT%TZ)] $*"; }
manifest(){ (cd "$1" && find . -path '*/Gene/filtered/*' -type f ! -name gdna_metrics.json | sort | xargs sha256sum | sha256sum | cut -c1-16); }
common=(--runThreadN 32 --dynamicThreadInterface 1 --genomeDir "$REF/star_index" --soloType CB_UMI_Simple --soloCBstart 1 --soloUMIstart 17 --soloCBlen 16 --soloUMIlen 12
  --soloBarcodeReadLength 0 --soloCBwhitelist "$REF/737K-fixed-rna-profiling.txt" --flex yes --soloFlexExpectedCellsPerTag 3000
  --soloSampleWhitelist "$REF/sample_whitelist_full_16.tsv" --soloProbeList "$REF/star_index/flex_probe_artifacts/probe_list.txt"
  --soloSampleProbes "$REF/probe-barcodes-fixed-rna-profiling-rna.txt" --soloSampleProbeOffset 68 --soloFlexAllowedTags "${ALLOWED:-$REF/sample_whitelist_full_16.tsv}"
  --limitIObufferSize 50000000 50000000 --outSJtype None --outSAMtype None --outSAMattributes None --soloFeatures Gene --soloCellFilter None
  --soloMultiMappers Rescue --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloStrand Unstranded
  --clipAdapterType CellRanger4 --alignEndsType Local --chimSegmentMin 0 --soloKeysCompat cr --soloHashScreenFile "$CACHE"
  --soloBucketMode ram --soloBucketCount 256 --flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0 --flexNoAlign 1 --crAssignConsumerThreads -1 --crAssignSearchThreads 1
  --readFilesType Binseq PE --readFilesCbqRangeMode range --readFilesIn "$(ls $CBQ/lane_00?.cbq | paste -sd,)")
o=$OUT/$ARM; rm -rf "$o"; mkdir -p "$o"; sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' || true; say "$ARM starting (binary $(sha256sum $STAR | cut -c1-16), cache $(basename $CACHE))"
/usr/bin/time -v -o "$T/tagfix_$ARM.time.txt" "$STAR" "${common[@]}" --soloFlexOutputPrefix "$o/per_sample" --outTmpDir "$o/_STARtmp" --outFileNamePrefix "$o/" "$@" > "$T/tagfix_$ARM.stdout" 2>&1
say "$ARM exit=$? wall=$(grep -aoP 'Elapsed.*: \K.*' $T/tagfix_$ARM.time.txt) out=$(manifest $o/per_sample)"
grep -aE "Hash screen: (KEEP|DENY|unmatched sample tag) %" "$o/Log.final.out"
python3 "$CONC" "$o/per_sample" star-$ARM > "$T/concordance_$ARM.txt" 2>&1; grep -vE '^#' "$T/concordance_$ARM.txt" | tail -1
python3 - "$o/per_sample" <<'PY'
import sys; sys.path.insert(0,"/mnt/pikachu/STAR-suite/docs/benchmarks/jax_matrix_20260904"); import concordance_vs_cr as cc
for bc,s in cc.DEFAULT_GROUPS:
    cx,cb,cg=cc.read_mex(cc.DEFAULT_CRROOT/bc/"count"/"sample_filtered_feature_bc_matrix"); qx,qb,qg=cc.read_mex(f"{sys.argv[1]}/{s[0]}/Gene/filtered"); ex,eb,eg=cc.read_mex(f"/storage/jax_tagfix/E_edgefix/per_sample/{s[0]}/Gene/filtered")
    gi={g:i for i,g in enumerate(cg)}; qi={g:i for i,g in enumerate(qg)}; ei={g:i for i,g in enumerate(eg)}
    genes=" ".join(f"{n} {int(qx[:,qi[g]].sum())}/{int(ex[:,ei[g]].sum())}/{int(cx[:,gi[g]].sum())}" for n,g in (("ZIC2","ENSG00000043355"),("APOE","ENSG00000130203"),("N4BP2","ENSG00000078177"),("H2BC15","ENSG00000277224"),("SOX2","ENSG00000181449")) if g in qi and g in gi)
    print(f"  {bc:18s} UMIs ours/CR = {qx.sum()/cx.sum():.4f} (E_edgefix {ex.sum()/cx.sum():.4f})   this/E/CR: {genes}")
PY
echo "finished=$(date -u +%FT%TZ)" > "$SENT"
