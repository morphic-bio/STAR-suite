#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
star="${STAR_BIN:-${repo}/core/legacy/source/STAR}"
caller="${CALL_FEATURES_BIN:-${repo}/core/features/process_features/call_features}"
genome="${MSK_MULTI_GENOME:-/storage/autoindex_110_44/bulk_index}"
tru=/storage/scRNAseq_output/whitelists/3M-february-2018_TRU.txt
nxt=/storage/scRNAseq_output/whitelists/3M-february-2018_NXT.txt
grna_ref=/mnt/pikachu/MSK-whitelists/ref_feature_geneBC_crispr.csv

for input in "$star" "$caller" "$genome" "$tru" "$nxt" "$grna_ref"; do
    [[ -e "$input" ]] || { echo "Missing test input: $input" >&2; exit 1; }
done

work=$(mktemp -d /tmp/larry_dominant_pf_multi.XXXXXX)
trap 'rm -rf "$work"' EXIT
bash "${repo}/tests/multi_feature/create_fixture_downsampled.sh" \
    "$work/fixture" 100000 5000 > "$work/fixture.log"
fixture="$work/fixture"
run="$work/run"
mkdir -p "$run"

printf '%s\n' \
    '[libraries]' \
    'fastqs,sample,library_type,feature_types,star_chemistry,star_whitelist,star_feature_ref,star_library_id,star_feature_caller' \
    "$fixture/mRNA,DE_30KO,Gene Expression,Gene Expression,TRU,,,gex_de," \
    "$fixture/PolyIII,DE_30KO,CRISPR Guide Capture,CRISPR Guide Capture,NXT,$nxt,$grna_ref,grna_de," \
    "$fixture/LARRY,DE_30KO,Custom,Custom,TRU,$tru,$fixture/ref_feature_larryBC.csv,larry_de,dominant" \
    > "$run/multi_config.csv"

"$star" --runThreadN 8 --genomeDir "$genome" \
    --readFilesIn "$fixture/mRNA/mRNA_L001_R2_001.fastq.gz" \
                  "$fixture/mRNA/mRNA_L001_R1_001.fastq.gz" \
    --readFilesCommand zcat --outFileNamePrefix "$run/" \
    --outSAMtype None --clipAdapterType CellRanger4 \
    --soloType CB_UMI_Simple --soloCBwhitelist "$tru" \
    --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
    --soloBarcodeReadLength 0 --soloFeatures GeneFull \
    --soloCellFilter EmptyDrops_CR --pfMultiConfig "$run/multi_config.csv" \
    --crMinUmi 2 --defaultCrCompat yes --dynamicThreadInterface 1 \
    --crAssignConsumerThreads -1 --crAssignSearchThreads 1 \
    > "$work/star_stdout.log" 2>&1 || {
        tail -60 "$work/star_stdout.log" >&2
        tail -60 "$run/Log.out" >&2
        exit 1
    }

calls="$run/outs/feature_analysis/larry_de/feature_calls.csv"
source_mex="$run/cr_assign/Custom/larry_de/LARRY/filtered"
[[ -s "$calls" && -s "$source_mex/matrix.mtx" ]] || {
    echo 'Integrated LARRY calls or filtered MEX missing' >&2; exit 1;
}
awk 'END { if (NR < 2) exit 1 }' "$calls" || {
    echo 'Integrated LARRY calls have no observed barcodes' >&2; exit 1;
}
grep -q $'feature_caller\tdominant' \
    "$run/cr_assign/Custom/larry_de/LARRY/pf_library_provenance.tsv"
"$caller" --guide-caller dominant --min_counts 1 --fraction 0 \
    --margin 1 --min-ratio 1 "$source_mex" "$work/control" \
    > "$work/control.log"
cmp "$calls" "$work/control/feature_calls.csv"

echo 'PASS: integrated STAR Custom/LARRY dominance matches process_features control'
