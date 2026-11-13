#!/usr/bin/env bash
# Run STARgrid, the Visium HD pipeline of STAR Suite, from FASTQ files to count
# matrices at 2, 8 and 16 micrometres under the strict, soft-expected, hard and
# gated-hard resolution policies.
#
#   --mode gex    3' gene expression: reads are aligned to --genome-dir.
#   --mode flex   probe-based (FFPE) Flex: reads are assigned to probes from a
#                 half-probe cache; no genome is loaded. The cache is built once
#                 from the probe set CSV and its GTF if --flex-reference lacks it.
#
# The barcode inputs are built from the Space Ranger slide-design file with
# scripts/make_visium_hd_contract.py when --contract does not hold them.
# The STAR options are the fixed recipes used for the STARgrid paper; STAR
# rejects most other values in these modes.
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

usage() {
    cat <<'EOF'
usage: run_visium_hd.sh --mode gex|flex --fastq-dir DIR --out-dir DIR
                        (--slide FILE | --contract DIR) [options]

Required:
  --mode gex|flex           3' gene expression or probe-based Flex
  --fastq-dir DIR           folder of *_R1_*.fastq.gz and matching *_R2_*.fastq.gz
  --out-dir DIR             output folder; matrices go to DIR/SpatialGex.out or
                            DIR/SpatialFlex.out
  --slide FILE              visium_hd_v1.slide from a Space Ranger installation
                            (lib/python/cellranger/barcodes/), used to build
                            the barcode inputs
  --contract DIR            folder written by scripts/make_visium_hd_contract.py;
                            built there from --slide if it is missing
                            (default with --slide: OUT_DIR/barcode_inputs)

GEX:
  --genome-dir DIR          STAR genome index of the reference

Flex:
  --flex-reference DIR      folder with flex_h01x2.khash and probe_gene_list.txt;
                            built there if missing (about 30 GB of memory and
                            30 GB of disk while building)
  --probe-set FILE          probe set CSV, to build the Flex reference
  --gtf FILE                GTF of the probe set's reference, to build it

Options:
  --sample NAME             use only FASTQ files whose names start with NAME
  --threads N               threads (default 16)
  --reads N                 read pairs in the FASTQ files (default: counted)
  --memory-fraction F       share of available memory for the spatial step
                            before it spills to disk (default 0.5)
  --star PATH               STAR binary (default: $STAR_BINARY, STAR on PATH,
                            then core/legacy/source/STAR in this repository)
  --dry-run                 print the commands without running them
EOF
}

die() {
    echo "error: $*" >&2
    exit 2
}

mode=
fastq_dir=
out_dir=
slide=
contract=
genome_dir=
flex_reference=
probe_set=
gtf=
sample=
threads=16
reads=
memory_fraction=0.5
star=
dry_run=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode) mode=$2; shift 2 ;;
        --fastq-dir) fastq_dir=$2; shift 2 ;;
        --out-dir) out_dir=$2; shift 2 ;;
        --slide) slide=$2; shift 2 ;;
        --contract) contract=$2; shift 2 ;;
        --genome-dir) genome_dir=$2; shift 2 ;;
        --flex-reference) flex_reference=$2; shift 2 ;;
        --probe-set) probe_set=$2; shift 2 ;;
        --gtf) gtf=$2; shift 2 ;;
        --sample) sample=$2; shift 2 ;;
        --threads) threads=$2; shift 2 ;;
        --reads) reads=$2; shift 2 ;;
        --memory-fraction) memory_fraction=$2; shift 2 ;;
        --star) star=$2; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; die "unknown option: $1" ;;
    esac
done

case "$mode" in
    gex|flex) ;;
    *) usage >&2; die "--mode must be gex or flex" ;;
esac
[[ -n "$fastq_dir" && -d "$fastq_dir" ]] || die "--fastq-dir must be an existing folder"
[[ -n "$out_dir" ]] || die "--out-dir is required"
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || die "--threads must be a positive integer"
[[ -z "$reads" || "$reads" =~ ^[1-9][0-9]*$ ]] || die "--reads must be a positive integer"
if [[ "$mode" == gex ]]; then
    [[ -n "$genome_dir" && -f "$genome_dir/SA" ]] || die "--mode gex needs --genome-dir with a STAR index"
else
    [[ -n "$flex_reference" ]] || die "--mode flex needs --flex-reference"
fi
if [[ -z "$contract" ]]; then
    [[ -n "$slide" ]] || die "give --slide or --contract"
    contract=$out_dir/barcode_inputs
fi

if [[ -z "$star" ]]; then
    if [[ -n "${STAR_BINARY:-}" ]]; then
        star=$STAR_BINARY
    elif command -v STAR >/dev/null 2>&1; then
        star=$(command -v STAR)
    else
        star=$repo_root/core/legacy/source/STAR
    fi
fi
[[ -x "$star" ]] || die "STAR binary not found or not executable: $star"

run() {
    printf '%q ' "$@"
    echo
    if [[ $dry_run -eq 0 ]]; then
        "$@"
    fi
}

# FASTQ pairs, lanes in name order.
r1_files=()
r2_files=()
while IFS= read -r r1; do
    r2=${r1/_R1_/_R2_}
    [[ -f "$r2" ]] || die "no R2 file for $r1"
    r1_files+=("$r1")
    r2_files+=("$r2")
done < <(find "$fastq_dir" -maxdepth 1 -type f -name "${sample}*_R1_*.fastq.gz" | LC_ALL=C sort)
[[ ${#r1_files[@]} -gt 0 ]] || die "no ${sample}*_R1_*.fastq.gz files in $fastq_dir"
r1_list=$(IFS=,; echo "${r1_files[*]}")
r2_list=$(IFS=,; echo "${r2_files[*]}")

if [[ $dry_run -eq 0 ]]; then
    mkdir -p "$out_dir"
    product=$out_dir/Spatial$([[ $mode == gex ]] && echo Gex || echo Flex).out
    [[ ! -e "$product" ]] || die "output exists; use a new --out-dir: $product"
    # Launchpad's Run in shell discards the terminal output; keep a copy.
    exec > >(tee -a "$out_dir/run_visium_hd.log") 2>&1
fi

# Barcode inputs.
contract_files=(barcode_coords.tsv bc1_full_oligos.txt bc2_full_oligos.txt bc1_whitelist.txt bc2_whitelist.txt)
missing=0
for name in "${contract_files[@]}"; do
    [[ -f "$contract/$name" ]] || missing=1
done
if [[ $missing -eq 1 ]]; then
    [[ -n "$slide" && -f "$slide" ]] || die "$contract lacks the barcode inputs; give --slide to build them"
    run python3 "$repo_root/scripts/make_visium_hd_contract.py" --slide "$slide" --out-dir "$contract"
fi

# Flex reference: the probe index, the half-probe cache and its gene list.
if [[ "$mode" == flex && ! ( -f "$flex_reference/flex_h01x2.khash" && -f "$flex_reference/probe_gene_list.txt" ) ]]; then
    [[ -n "$probe_set" && -f "$probe_set" && -n "$gtf" && -f "$gtf" ]] \
        || die "$flex_reference lacks flex_h01x2.khash; give --probe-set and --gtf to build it"
    command -v g++ >/dev/null 2>&1 || die "building the Flex reference needs g++"
    build=$flex_reference/build
    if [[ $dry_run -eq 0 ]]; then
        mkdir -p "$build/index" "$build/logs" "$build/cache_logs"
        rm -f "$build/h01x2.bin"
        # A one-contig placeholder genome: the index holds only the probes.
        printf '>B1_UNUSED_BASE\n%s\n' "$(printf 'A%.0s' {1..100})" > "$build/unused_base.fa"
        # Visium HD has no sample tags; cache generation needs one placeholder.
        printf 'BC001\tACGTACGT\n' > "$build/sample.tsv"
        printf 'ACGTACGT\tACGTACGT\tBC001\n' > "$build/sample_probes.tsv"
        printf 'ACGTACGTACGTACGT\n' > "$build/cell_barcode.txt"
    fi
    run "$star" --runMode genomeGenerate --runThreadN "$threads" --genomeDir "$build/index" \
        --genomeFastaFiles "$build/unused_base.fa" --sjdbGTFfile "$gtf" \
        --flexGeneProbeSet "$probe_set" --sjdbOverhang 49 --genomeSAindexNbases 8 \
        --genomeChrBinNbits 10 --outFileNamePrefix "$build/logs/"
    run "$star" --runMode hashCacheGenerate --runThreadN "$threads" --genomeDir "$build/index" \
        --soloType CB_UMI_Simple --soloCBlen 16 --soloUMIlen 12 --soloUMIstart 17 --soloCBstart 1 \
        --soloBarcodeReadLength 0 --soloCBwhitelist "$build/cell_barcode.txt" --flex yes \
        --soloProbeList "$build/index/probe_gene_list.txt" \
        --soloSampleWhitelist "$build/sample.tsv" --soloSampleProbes "$build/sample_probes.tsv" \
        --soloSampleProbeOffset 68 --soloFlexAllowedTags "$build/sample.tsv" \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR --soloCellFilter None --soloFeatures Gene --alignEndsType Local \
        --soloStrand Unstranded --soloKeysCompat cr --outSAMtype None --hashCacheTiers H0,H1X2 \
        --hashCacheOutput "$build/h01x2.bin" --outFileNamePrefix "$build/cache_logs/" \
        --soloRemoveDeprecated No
    core=$repo_root/core/legacy/source
    run g++ -O3 -std=c++11 -I"$core" -I"$repo_root/flex/source" -o "$build/flex_hash_cache_pack" \
        "$repo_root/flex/tools/hash_screen_replay/flex_hash_cache_pack.cpp" \
        "$core/FlexHashCacheStorage.cpp" "$core/FlexProbePairKhash.cpp"
    run "$build/flex_hash_cache_pack" "$build/h01x2.bin" "$flex_reference/flex_h01x2.khash"
    run cp "$build/index/probe_gene_list.txt" "$flex_reference/probe_gene_list.txt"
    run rm -f "$build/h01x2.bin"
fi

# Read capacities: STAR stops if the reads or their candidate locations exceed
# them. Candidate locations have stayed below 1.2 per read pair.
if [[ -z "$reads" ]]; then
    if [[ $dry_run -eq 1 ]]; then
        echo "# read pairs are counted before STAR runs; the capacities below are placeholders"
        reads=1
    else
        echo "counting read pairs in ${#r1_files[@]} R1 file(s)"
        decompress=(gzip -dc)
        command -v pigz >/dev/null 2>&1 && decompress=(pigz -dc -p 2)
        counts=$(mktemp -d)
        pids=()
        for i in "${!r1_files[@]}"; do
            ("${decompress[@]}" "${r1_files[$i]}" | wc -l > "$counts/$i") &
            pids+=($!)
        done
        for pid in "${pids[@]}"; do
            wait "$pid" || die "could not read the R1 files"
        done
        lines=$(cat "$counts"/* | awk '{ total += $1 } END { print total }')
        rm -rf "$counts"
        reads=$((lines / 4))
        echo "$reads" > "$out_dir/read_pairs.txt"
    fi
fi
max_capacity=4294967295
(( reads < max_capacity )) || die "STARgrid handles fewer than $max_capacity read pairs per run"
candidates=$(( (reads * 3 + 1) / 2 ))
(( candidates <= max_capacity )) || candidates=$max_capacity

read_command=(--readFilesCommand zcat)
command -v pigz >/dev/null 2>&1 && read_command=(--readFilesCommand pigz -dc -p 2)
common=(
    --runThreadN "$threads"
    --readFilesIn "$r2_list" "$r1_list"
    --readFilesBgzfMode auto --bgzfReaderThreads 0 "${read_command[@]}"
    --outFileNamePrefix "$out_dir/"
)
spatial=(
    --soloSpatialBarcodeContract "$contract"
    --soloSpatialBc1Oligos "$contract/bc1_full_oligos.txt"
    --soloSpatialBc2Oligos "$contract/bc2_full_oligos.txt"
    --soloSpatialAssignmentProducts all --soloSpatialBinSizes 2,8,16
    --soloSpatialExpectedReads "$reads" --soloSpatialExpectedCandidates "$candidates"
    --soloSpatialMemoryFraction "$memory_fraction" --soloSpatialOverflowPolicy Spill
)
if [[ "$mode" == gex ]]; then
    run "$star" "${common[@]}" --genomeDir "$genome_dir" \
        --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloType None \
        --soloFeatures GeneFull --soloCrGexFeature GeneFull --soloCrMultimapRescue yes \
        --soloCrMultimapRescueEvidence compatibility --soloCrMultimapRescueIntronic auto \
        --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR --soloMultiMappers Unique \
        --soloStrand Forward --soloCellFilter None --outSAMtype None \
        --soloSpatialGexIntegrated yes "${spatial[@]}"
else
    run "$star" "${common[@]}" \
        --outSAMtype None --outSJtype None --chimSegmentMin 0 --outSAMattributes None \
        --soloType CB_UMI_Complex --soloBarcodeMate 0 \
        --soloCBwhitelist "$contract/bc1_whitelist.txt" "$contract/bc2_whitelist.txt" \
        --soloStrand Unstranded --soloFeatures Gene --soloUMIdedup 1MM_CR \
        --soloCBmatchWLtype EditDist_2 --soloCBposition 0_11_0_24 0_25_0_38 \
        --soloUMIposition 0_0_0_8 --soloCellFilter None --soloUMIfiltering MultiGeneUMI_CR \
        --soloMultiMappers Rescue --soloSkipProcessing yes --soloKeysCompat cr \
        --soloProbeList "$flex_reference/probe_gene_list.txt" --soloInlineHashMode yes \
        --soloHashScreenFile "$flex_reference/flex_h01x2.khash" --flex yes --soloProbeMismatch 1 \
        --soloRemoveDeprecated No --flexPipeline yes --flexPipelineNTriage 0 \
        --flexPipelineNSolo 0 --flexNoAlign 1 --soloRunFlexFilter no \
        --soloFlexMinimalMemory no --soloSpatialFlexIntegrated yes "${spatial[@]}"
fi

if [[ $dry_run -eq 0 ]]; then
    echo "matrices: $product/<policy>/square_{002,008,016}um/"
fi
