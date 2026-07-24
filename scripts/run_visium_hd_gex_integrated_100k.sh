#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 FRESH_OUTPUT_ROOT" >&2
    exit 2
fi

output_root=$1
if [[ -e "$output_root" ]]; then
    echo "error: output root already exists: $output_root" >&2
    exit 2
fi

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
star_binary=${STAR_BINARY:-$repo_root/core/legacy/source/STAR}
fixture=${VISIUM_HD_GEX_FIXTURE:-${VISIUM_HD_GEX_100K_FIXTURE:-/mnt/pikachu/star-spatial/10x/visium_hd_3prime_human_ovarian_ff_min_depth/downsample_100k_v1}}
genome_dir=${VISIUM_HD_GEX_GENOME_DIR:-/mnt/pikachu/star-spatial/references/refdata-gex-GRCh38-2024-A_STAR-2.7.11a/star}
bc1_oligos=${VISIUM_HD_BC1_OLIGOS:-/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc1_full_oligos.txt}
bc2_oligos=${VISIUM_HD_BC2_OLIGOS:-/storage/star-spatial/runs/cleanroom_hd_mouse_brain/slide_oligos/bc2_full_oligos.txt}
barcode_contract=${VISIUM_HD_BARCODE_CONTRACT:-/storage/star-spatial/runs/cleanroom_hd_mouse_brain/barcode_contract}
threads=${STAR_THREADS:-16}
overflow_policy=${SPATIAL_OVERFLOW_POLICY:-Fail}
memory_fraction=${SPATIAL_MEMORY_FRACTION:-0.8}
spill_high_water=${SPATIAL_SPILL_HIGH_WATER_CANDIDATES:-0}
diagnostic_sidecar=${SPATIAL_DIAGNOSTIC_SIDECAR_PREFIX:-}
expected_reads=${SPATIAL_EXPECTED_READS:-100000}
expected_candidates=${SPATIAL_EXPECTED_CANDIDATES:-111744}

for path in "$star_binary" "$fixture/checksums.sha256" "$genome_dir/Genome" \
    "$genome_dir/SA" "$genome_dir/SAindex" "$genome_dir/geneInfo.tab" \
    "$bc1_oligos" "$bc2_oligos" "$barcode_contract/barcode_coords.tsv" \
    "$barcode_contract/raw28_barcode_coords.tsv"; do
    if [[ ! -e "$path" ]]; then
        echo "error: required source does not exist: $path" >&2
        exit 2
    fi
done

if pgrep -x STAR >/dev/null 2>&1; then
    echo "error: another STAR process is active; benchmarks must be serialized" >&2
    exit 2
fi

mkdir -p "$output_root"
(
    cd "$fixture"
    sha256sum -c checksums.sha256
) >"$output_root/fixture_checksums.log"

sha256sum "$star_binary" >"$output_root/star_binary.sha256"
sha256sum "$genome_dir/Genome" "$genome_dir/SA" "$genome_dir/SAindex" \
    "$genome_dir/geneInfo.tab" >"$output_root/index.sha256"
stat -c '%n\t%s\t%y' "$genome_dir/Genome" "$genome_dir/SA" \
    "$genome_dir/SAindex" "$genome_dir/geneInfo.tab" >"$output_root/index.stat.tsv"
git -C "$repo_root" status --short >"$output_root/git_status.txt"
git -C "$repo_root" rev-parse HEAD >"$output_root/git_head.txt"

r2_files=
r1_files=
for lane in 001 002 003 004; do
    stem="$fixture/fastqs/Visium_HD_3prime_Human_Ovarian_Cancer_FF_Min_Depth_S1_L${lane}"
    r2="$stem"_R2_001.fastq.gz
    r1="$stem"_R1_001.fastq.gz
    [[ -f "$r2" && -f "$r1" ]] || {
        echo "error: FASTQ pair is incomplete for lane $lane" >&2
        exit 2
    }
    r2_files=${r2_files:+$r2_files,}$r2
    r1_files=${r1_files:+$r1_files,}$r1
done

spatial_optional_args=()
if [[ "$spill_high_water" != 0 ]]; then
    spatial_optional_args+=(--soloSpatialSpillHighWaterCandidates "$spill_high_water")
fi
if [[ -n "$diagnostic_sidecar" && "$diagnostic_sidecar" != "-" && "$diagnostic_sidecar" != "None" ]]; then
    spatial_optional_args+=(--soloSpatialFeatureSidecar "$diagnostic_sidecar")
fi

star_pid=0
terminate_star() {
    if [[ $star_pid -gt 0 ]] && kill -0 "$star_pid" 2>/dev/null; then
        kill -TERM -- "-$star_pid" 2>/dev/null || true
        wait "$star_pid" 2>/dev/null || true
    fi
}
trap 'terminate_star; exit 130' INT
trap 'terminate_star; exit 143' TERM

start_epoch=$(date +%s)
setsid /usr/bin/time -v "$star_binary" \
    --runThreadN "$threads" \
    --genomeDir "$genome_dir" \
    --readFilesIn "$r2_files" "$r1_files" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$output_root/star" \
    --clipAdapterType CellRanger4 \
    --outFilterScoreMin 30 \
    --soloType None \
    --soloFeatures GeneFull \
    --soloCrGexFeature GeneFull \
    --soloCrMultimapRescue yes \
    --soloCrMultimapRescueEvidence annotated \
    --soloCrMultimapRescueIntronic auto \
    --soloUMIdedup 1MM_CR \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloMultiMappers Unique \
    --soloStrand Forward \
    --soloCellFilter None \
    --outSAMtype None \
    --soloSpatialGexIntegrated yes \
    --soloSpatialBarcodeContract "$barcode_contract" \
    --soloSpatialBc1Oligos "$bc1_oligos" \
    --soloSpatialBc2Oligos "$bc2_oligos" \
    --soloSpatialAssignmentProducts all \
    --soloSpatialBinSizes 2,8,16 \
    --soloSpatialExpectedReads "$expected_reads" \
    --soloSpatialExpectedCandidates "$expected_candidates" \
    --soloSpatialMemoryFraction "$memory_fraction" \
    --soloSpatialOverflowPolicy "$overflow_policy" \
    "${spatial_optional_args[@]}" \
    >"$output_root/benchmark.stdout.log" \
    2>"$output_root/benchmark.stderr.log" &
star_pid=$!

set +e
wait "$star_pid"
exit_code=$?
set -e
star_pid=0
finish_epoch=$(date +%s)

{
    echo "exit_code=$exit_code"
    echo "started_unix_seconds=$start_epoch"
    echo "finished_unix_seconds=$finish_epoch"
    echo "wall_seconds=$((finish_epoch - start_epoch))"
    echo "overflow_policy=$overflow_policy"
    echo "memory_fraction=$memory_fraction"
    echo "expected_reads=$expected_reads"
    echo "expected_candidates=$expected_candidates"
    echo "spill_high_water_candidates_per_thread=$spill_high_water"
    echo "diagnostic_sidecar=${diagnostic_sidecar:--}"
} >"$output_root/BENCHMARK_SUMMARY.txt"

if [[ $exit_code -ne 0 ]]; then
    exit "$exit_code"
fi
if [[ ! -f "$output_root/starSpatialGex.out/RUN_COMPLETE" ]]; then
    echo "error: STAR exited successfully without integrated RUN_COMPLETE" >&2
    exit 1
fi

touch "$output_root/RUN_COMPLETE"
