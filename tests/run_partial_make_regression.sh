#!/usr/bin/env bash
# Build each requested target from tracked source bytes, not prior build output.
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
out=${OUT_ROOT:-$(mktemp -d /tmp/star-partial-make.XXXXXX)}
jobs=${MAKE_JOBS:-8}
cases=("$@")
if [[ ${#cases[@]} -eq 0 ]]; then
    cases=(core core-static core-htslib release-companion-tools
           star-feature-call feature-barcodes-tools process-features-lib yremove-tools)
fi
mkdir -p "$out"
out=$(cd "$out" && pwd)
if [[ -e "$out/summary.tsv" ]]; then
    echo "ERROR: use a fresh OUT_ROOT; $out already contains a summary." >&2
    exit 2
fi
printf 'case\tstatus\tlog\n' > "$out/summary.tsv"

# Export working-tree versions of tracked files, including local Makefile fixes.
# Do not copy ignored objects/binaries or large local datasets. New test files
# need not be present in the export: this runner executes outside the test trees.
git -C "$root" ls-files -z | tar -C "$root" --null -T - -cf "$out/source.tar"

for case_name in "${cases[@]}"; do
    source_dir="$out/$case_name/source"
    log="$out/$case_name/build.log"
    mkdir -p "$source_dir"
    tar -xf "$out/source.tar" -C "$source_dir"
    src="$source_dir/core/legacy/source"
    args=()
    artifacts=()
    case "$case_name" in
        core|core-portable|core-static)
            args=("$case_name")
            artifacts=(core/legacy/source/STAR)
            ;;
        core-external-htslib)
            args=(core HTSLIB=external)
            artifacts=(core/legacy/source/STAR)
            ;;
        core-htslib)
            args=(core-htslib)
            artifacts=(core/legacy/source/htslib/libhts.a)
            ;;
        release-companion-tools)
            args=(release-companion-tools)
            artifacts=(core/legacy/source/transcriptvb_finalize
                       core/legacy/source/trim_qc_fastq core/legacy/source/trim_qc_merge)
            ;;
        star-feature-call)
            args=(star-feature-call)
            artifacts=(core/legacy/source/star_feature_call)
            ;;
        feature-barcodes-tools)
            args=(feature-barcodes-tools)
            artifacts=(core/features/process_features/assignBarcodes
                       core/features/process_features/demux_fastq
                       core/features/process_features/demux_bam
                       core/features/process_features/call_features)
            ;;
        process-features-lib)
            args=(process-features-lib)
            artifacts=(core/features/process_features/libprocess_features.a)
            ;;
        yremove-tools)
            args=(yremove-tools)
            artifacts=(core/features/yremove_fastq/tools/remove_y_reads/remove_y_reads)
            ;;
        *) echo "ERROR: unknown partial build case: $case_name" >&2; exit 2 ;;
    esac
    echo "Building $case_name from a fresh export; log: $log"
    if ! make -C "$source_dir" -j"$jobs" "${args[@]}" > "$log" 2>&1; then
        printf '%s\tFAIL\t%s\n' "$case_name" "$log" >> "$out/summary.tsv"
        tail -n 60 "$log" >&2
        exit 1
    fi
    for artifact in "${artifacts[@]}"; do
        test -s "$source_dir/$artifact"
        case "$artifact" in
            *.a) ar t "$source_dir/$artifact" >> "$log" ;;
            *) test -x "$source_dir/$artifact" ;;
        esac
    done
    if [[ -x "$src/STAR" ]]; then
        "$src/STAR" --version >> "$log"
    fi
    printf '%s\tPASS\t%s\n' "$case_name" "$log" >> "$out/summary.tsv"
done
echo "PASS: partial Make targets; results: $out/summary.tsv"
