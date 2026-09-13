#!/bin/bash
# STAR Suite 1.9.4: the half-probe (H1X2) route is the Flex default and every other
# route is LEGACY. Checks, on the JAX SC2300771 8-lane x 100K-pair fixture:
#   - plain --flex yes selects the fused no-alignment route and is byte-identical to
#     the explicit route flags (and, if STAR_REF_BIN is set, to that reference binary)
#   - each legacy route stops with an explanation unless --flexLegacy yes is given
#   - --flexLegacy yes and --runMode hashCacheGenerate are not blocked
# Override the data locations with the variables below.
set -u
unset STAR_FLEX_HASH_SCREEN_CACHE
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
S194="${STAR_BIN:-${REPO_ROOT}/core/legacy/source/STAR}"
S193="${STAR_REF_BIN:-}"
IFMT="${FLEX_HALF_CACHE_DIR:-/mnt/pikachu/storage_offload_20260912/inputfmt_20260910/stage}"
FQ="${FLEX_FIXTURE_DIR:-/mnt/pikachu/storage_symlinked_20260912/downsampled_100K/SC2300771}"
V="${TEST_WORKDIR:-/tmp/star_suite_flex_v194_default_route_$$}"
[[ -x "$S194" ]] || { echo "ERROR: STAR binary not executable: $S194" >&2; exit 2; }
[[ -f "$IFMT/model_h01x2_cache.half.khash" ]] || { echo "ERROR: half-probe cache missing under FLEX_HALF_CACHE_DIR=$IFMT" >&2; exit 2; }
[[ -d "$FQ" ]] || { echo "ERROR: fixture missing: FLEX_FIXTURE_DIR=$FQ" >&2; exit 2; }
mkdir -p "$V"
P=SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5
EMPTY=$V/empty_genome_index; rm -rf $EMPTY; mkdir -p $EMPTY
R2=(); R1=(); for L in L001 L002 L003 L004 L005 L006 L007 L008; do
  R2+=("$FQ/${P}_${L}_R2_001.fastq.gz"); R1+=("$FQ/${P}_${L}_R1_001.fastq.gz"); done
base=(--runThreadN 32 --genomeDir "$EMPTY" --soloType CB_UMI_Simple --soloCBstart 1 --soloUMIstart 17
  --soloCBlen 16 --soloUMIlen 12 --soloBarcodeReadLength 0
  --soloCBwhitelist /mnt/pikachu/hdd_320k_stage/737K-fixed-rna-profiling.txt --flex yes --soloRemoveDeprecated No
  --soloSampleWhitelist /mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv
  --soloProbeList "$IFMT/model_gene_ids.txt" --soloFlexFilteredGeneList "$IFMT/included_gene_ids.txt"
  --soloSampleProbes /mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt
  --soloSampleProbeOffset 68 --soloFlexAllowedTags /mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv
  --limitIObufferSize 50000000 50000000 --outSJtype None --outSAMattributes None
  --soloFeatures Gene --soloCellFilter None --soloMultiMappers Rescue
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR
  --soloStrand Unstranded --clipAdapterType CellRanger4 --alignEndsType Local --chimSegmentMin 0
  --soloKeysCompat cr --soloProbeMismatch 1 --soloSampleTagMismatch 1 --soloSampleSearchNearby no
  --soloBucketMode ram --soloBucketCount 256 --dynamicThreadInterface 1
  --crAssignConsumerThreads -1 --crAssignSearchThreads 1
  --soloRunFlexFilter yes --soloFlexCellCaller tag-aware
  --soloCellFilterMitochondrialGenes /mnt/pikachu/star_suite_paper/analysis/full320k_star_deprecated_20260909/completed/config/mitochondrial_gene_ids.tsv
  --soloCellFilterBootstrapThreads 32 --soloFlexFatalOnError no --soloFlexEdFdrThreshold 0.01
  --soloFlexEdNiters 0 --soloFlexKeepCBTag yes --soloFlexDebugTagLog no --soloFlexInvariantChecks no
  --soloFlexDecisionSidecar - --readFilesBgzfMode off
  --readFilesIn "$(IFS=,; echo "${R2[*]}")" "$(IFS=,; echo "${R1[*]}")")
CACHE=(--soloHashScreenFile "$IFMT/model_h01x2_cache.half.khash")
EXPLICIT=(--outSAMtype None --flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0 --flexNoAlign 1)
pass=0; fail=0
ok(){ echo "  PASS  $*"; pass=$((pass+1)); }
bad(){ echo "  FAIL  $*"; fail=$((fail+1)); }
run(){ # name binary args...
  local n=$1 b=$2; shift 2; local o=$V/$n; rm -rf $o; mkdir -p $o
  "$b" "${base[@]}" "$@" --soloFlexOutputPrefix "$o/per_sample" --soloFlexDebugOutputDir "$o/caller_diagnostics" \
      --outTmpDir "$o/_STARtmp" --outFileNamePrefix "$o/" > "$o/stdout.txt" 2>&1
  echo $? > "$o/rc"; rm -rf "$o/_STARtmp"
}
logs(){ cat "$V/$1/stdout.txt" "$V/$1/Log.out" 2>/dev/null; }
digest(){ (cd "$V/$1" && find . -type f ! -name 'Log.*' ! -name stdout.txt ! -name rc | sort | xargs -r sha256sum) ; }

echo "== positive runs =="
run A $S194 "${CACHE[@]}" "${EXPLICIT[@]}"
run B $S194 "${CACHE[@]}"
POS="A B L2"
if [[ -n "$S193" ]]; then run C $S193 "${CACHE[@]}" "${EXPLICIT[@]}"; POS="A B C L2"; fi
run L2 $S194 "${CACHE[@]}" "${EXPLICIT[@]}" --flexLegacy yes
for n in $POS; do [ "$(cat $V/$n/rc)" = 0 ] && ok "$n completed (rc 0)" || bad "$n rc=$(cat $V/$n/rc): $(tail -2 $V/$n/stdout.txt)"; done
nfiles=$(digest A | wc -l)
for n in ${POS#A }; do
  if diff <(digest A) <(digest $n) >/dev/null; then ok "$n byte-identical to A ($nfiles files)"
  else bad "$n differs from A: $(diff <(digest A) <(digest $n) | grep -c '^[<>]') lines"; fi
done
logs B | grep -q 'flexProbeRoute=half-probe H1X2 (default)' && ok "B logs half-probe default route" || bad "B missing route log"
logs B | grep -q 'flexNoAlign=1 (half-probe default)' && ok "B flexNoAlign defaulted to 1" || bad "B flexNoAlign not defaulted"
logs B | grep -q 'flexPipelineNTriage=0 (half-probe default)' && ok "B NTriage defaulted to 0" || bad "B NTriage not defaulted"
logs B | grep -q 'flexPipelineNSolo=0 (half-probe default)' && ok "B NSolo defaulted to 0" || bad "B NSolo not defaulted"
logs B | grep -q 'outSAMtype=None (half-probe default' && ok "B outSAMtype defaulted to None" || bad "B outSAMtype not defaulted"
logs B | grep -q 'Flex probe route: half-probe H1X2 (1.9.4 default)' && ok "B loaded H1X2 cache on default route" || bad "B route-after-load log missing"
logs B | grep -qi 'fused no-genome threads started' && ok "B ran the fused no-genome pipeline" || bad "B did not run fused no-genome"

echo "== runs that must stop with an explanation =="
neg(){ # name pattern args...
  local n=$1 pat=$2; shift 2; run $n $S194 "$@"
  if [ "$(cat $V/$n/rc)" != 0 ] && logs $n | grep -qF -- "$pat"; then ok "$n stopped: '$pat'"
  else bad "$n rc=$(cat $V/$n/rc), expected '$pat'; got: $(logs $n | grep -m1 -E 'EXITING' | cut -c1-120)"; fi
}
neg N1_nocache      "requires a half-probe (H1X2) cache, but no cache was given"
neg N2_h0h1cache    "requires a half-probe (H1X2) cache" --soloHashScreenFile /mnt/pikachu/hdd_320k_stage/h01_cache.bin
neg N3_bam          "--outSAMtype BAM (SAM/BAM output needs genomic alignment)" "${CACHE[@]}" --outSAMtype BAM Unsorted
neg N4_noalign0     "--flexNoAlign 0 (aligning hash-screen misses) is a LEGACY" "${CACHE[@]}" --flexNoAlign 0
neg N5_nohash       "--no-hash-screen yes is a LEGACY" "${CACHE[@]}" --no-hash-screen yes
run N6_spatial $S194 "${CACHE[@]}" --soloSpatialFlexIntegrated yes
if logs N6_spatial | grep -qF -- "spatial Flex is not yet migrated"; then ok "N6_spatial stopped: spatial legacy check"
elif logs N6_spatial | grep -qF -- "lacks an exact immutable source revision"; then
  echo "  DEFERRED  N6_spatial: stopped earlier by the spatial source-revision guard (uncommitted build); re-check on the release build"
else bad "N6_spatial unexpected: $(logs N6_spatial | grep -m1 EXITING | cut -c1-120)"; fi
neg N7_badlegacy    "unrecognized option in --flexLegacy=maybe" "${CACHE[@]}" --flexLegacy maybe

echo "== runs that must get PAST the check (may fail later for lack of a genome) =="
run L1_legacy_nohash $S194 --flexLegacy yes --no-hash-screen yes
if ! logs L1_legacy_nohash | grep -q 'is a LEGACY Flex route in STAR Suite 1.9.4'; then ok "L1 legacy + no-hash-screen not blocked by the check (rc $(cat $V/L1_legacy_nohash/rc))"
else bad "L1 blocked by the check despite --flexLegacy yes"; fi
run G1_generate $S194 --runMode hashCacheGenerate --hashCacheOutput $V/G1_generate/cache.bin
if ! logs G1_generate | grep -qE 'requires a half-probe|is a LEGACY Flex route'; then ok "G1 cache generation not blocked by the check (rc $(cat $V/G1_generate/rc))"
else bad "G1 cache generation blocked by the check"; fi

echo "== summary: $pass passed, $fail failed =="
exit $fail
