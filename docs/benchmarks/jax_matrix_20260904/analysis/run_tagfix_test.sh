#!/bin/bash
# Validate the sample-tag tolerance fix on JAX CBQ (no-align, all NVMe):
#   A: fix on  (defaults: --soloSampleTagMismatch 1, --soloSampleSearchNearby yes)
#   B: fix off (--soloSampleTagMismatch 0 --soloSampleSearchNearby no) -> must hash 2d2fead4b43f2403
set -u; export PATH="$HOME/.cargo/bin:$HOME/.local/bin:$PATH"
T=/mnt/pikachu/nvme_jax_v184; STAR=/mnt/pikachu/STAR-suite-fix-sampletag-20260905/core/legacy/source/STAR; REF=/home/lhhung/jax_stage_20260903/ref; CBQ=/home/lhhung/jax_stage_20260903/cbq; OUT=/storage/jax_tagfix
CONC=/mnt/pikachu/STAR-suite/docs/benchmarks/jax_matrix_20260904/concordance_vs_cr.py; SENT="$T/tagfix.done"; rm -f "$SENT"; mkdir -p "$OUT"
say(){ echo "[$(date -u +%FT%TZ)] $*"; }
manifest(){ (cd "$1" && find . -path '*/Gene/filtered/*' -type f ! -name gdna_metrics.json | sort | xargs sha256sum | sha256sum | cut -c1-16); }
common=(--runThreadN 32 --dynamicThreadInterface 1 --genomeDir "$REF/star_index" --soloType CB_UMI_Simple --soloCBstart 1 --soloUMIstart 17 --soloCBlen 16 --soloUMIlen 12
  --soloBarcodeReadLength 0 --soloCBwhitelist "$REF/737K-fixed-rna-profiling.txt" --flex yes --soloFlexExpectedCellsPerTag 3000
  --soloSampleWhitelist "$REF/sample_whitelist_full_16.tsv" --soloProbeList "$REF/star_index/flex_probe_artifacts/probe_list.txt"
  --soloSampleProbes "$REF/probe-barcodes-fixed-rna-profiling-rna.txt" --soloSampleProbeOffset 68 --soloFlexAllowedTags "$REF/sample_whitelist_full_16.tsv"
  --limitIObufferSize 50000000 50000000 --outSJtype None --outSAMtype None --outSAMattributes None --soloFeatures Gene --soloCellFilter None
  --soloMultiMappers Rescue --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR --soloStrand Unstranded
  --clipAdapterType CellRanger4 --alignEndsType Local --chimSegmentMin 0 --soloKeysCompat cr --soloHashScreenFile "$REF/h01_cache.bin"
  --soloBucketMode ram --soloBucketCount 256 --flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0 --flexNoAlign 1 --crAssignConsumerThreads -1 --crAssignSearchThreads 1
  --readFilesType Binseq PE --readFilesCbqRangeMode range --readFilesIn "$(ls $CBQ/lane_00?.cbq | paste -sd,)")
run(){ local arm="$1"; shift; local o="$OUT/$arm"; rm -rf "$o"; mkdir -p "$o"; sudo -n sh -c 'sync; echo 3 > /proc/sys/vm/drop_caches' || true; say "$arm starting"
  /usr/bin/time -v -o "$T/tagfix_$arm.time.txt" "$STAR" "${common[@]}" --soloFlexOutputPrefix "$o/per_sample" --outTmpDir "$o/_STARtmp" --outFileNamePrefix "$o/" "$@" > "$T/tagfix_$arm.stdout" 2>&1
  local rc=$?; say "$arm exit=$rc wall=$(grep -aoP 'Elapsed.*: \K.*' $T/tagfix_$arm.time.txt) out=$(manifest $o/per_sample)"; grep -a "Flex sample tag detection" "$o/Log.out"; grep -aE "Hash screen: (KEEP|unmatched sample tag) " "$o/Log.final.out"; }
run A_fix_on
run B_fix_off --soloSampleTagMismatch 0 --soloSampleSearchNearby no
say "concordance (fix on) vs Cell Ranger"; python3 "$CONC" "$OUT/A_fix_on/per_sample" star-tagfix-on > "$T/concordance_tagfix_on.txt" 2>&1; grep -vE '^#' "$T/concordance_tagfix_on.txt"
python3 - "$OUT/A_fix_on/per_sample" <<'PY'
import sys,scipy.io; sys.path.insert(0,"/mnt/pikachu/STAR-suite/docs/benchmarks/jax_matrix_20260904"); import concordance_vs_cr as cc
for bc,s in cc.DEFAULT_GROUPS:
    cx,_,_=cc.read_mex(cc.DEFAULT_CRROOT/bc/"count"/"sample_filtered_feature_bc_matrix"); qx,_,_=cc.read_mex(f"{sys.argv[1]}/{s[0]}/Gene/filtered")
    print(f"  {bc:18s} UMIs ours/CR = {qx.sum()/cx.sum():.4f}   (before fix: 0.9737/0.9757/0.9764/0.9731)")
PY
echo "finished=$(date -u +%FT%TZ)" > "$SENT"
