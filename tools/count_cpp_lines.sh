#!/bin/bash
# Count C/C++ lines and files in STAR-suite vs upstream STAR.
#
# Layout support:
# - Modular tree: core/legacy/source + core/features
# - Upstream tree: source
#
# Exclusions:
# - htslib
# - opal
#
# Notes:
# - Default upstream ref is STAR tag 2.7.11b.
# - If the ref exists locally, comparison uses the local git object database.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)"
UPSTREAM_REF="2.7.11b"
SKIP_UPSTREAM="No"

usage() {
    cat << EOF
Usage: $0 [options]

Options:
  --upstream-ref REF    Git ref to compare against (default: ${UPSTREAM_REF})
  --skip-upstream       Skip upstream comparison (local counts only)
  -h, --help            Show help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --upstream-ref) UPSTREAM_REF="$2"; shift 2 ;;
        --skip-upstream) SKIP_UPSTREAM="Yes"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done

echo "=========================================="
echo "C/C++ Code Statistics: STAR-suite vs Upstream"
echo "=========================================="
echo ""

local_source_roots() {
    if [[ -d "$REPO_ROOT/core/legacy/source" ]]; then
        printf '%s\n' "$REPO_ROOT/core/legacy/source"
    elif [[ -d "$REPO_ROOT/source" ]]; then
        printf '%s\n' "$REPO_ROOT/source"
    fi

    if [[ -d "$REPO_ROOT/core/features" ]]; then
        printf '%s\n' "$REPO_ROOT/core/features"
    fi
}

ref_source_roots() {
    local ref="$1"

    if git -C "$REPO_ROOT" cat-file -e "${ref}:core/legacy/source" 2>/dev/null; then
        printf '%s\n' "core/legacy/source"
    elif git -C "$REPO_ROOT" cat-file -e "${ref}:source" 2>/dev/null; then
        printf '%s\n' "source"
    fi

    if git -C "$REPO_ROOT" cat-file -e "${ref}:core/features" 2>/dev/null; then
        printf '%s\n' "core/features"
    fi
}

primary_ref_source_prefix() {
    local ref="$1"
    if git -C "$REPO_ROOT" cat-file -e "${ref}:core/legacy/source" 2>/dev/null; then
        echo "core/legacy/source"
    elif git -C "$REPO_ROOT" cat-file -e "${ref}:source" 2>/dev/null; then
        echo "source"
    else
        echo ""
    fi
}

collect_source_files_local() {
    local roots=()
    local root
    mapfile -t roots < <(local_source_roots)

    for root in "${roots[@]}"; do
        find "$root" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" \) \
            ! -path "*/htslib/*" \
            ! -path "*/opal/*" \
            -print0
    done
}

collect_source_files_ref() {
    local ref="$1"
    local roots=()
    local root
    mapfile -t roots < <(ref_source_roots "$ref")

    for root in "${roots[@]}"; do
        git -C "$REPO_ROOT" ls-tree -r --name-only "$ref" -- "$root"
    done | grep -E '\.(c|cpp|h|hpp)$' | grep -vE '/(htslib|opal)/' || true
}

matches_any() {
    local file="$1"
    shift
    local pattern
    for pattern in "$@"; do
        if [[ "$file" == $pattern ]]; then
            return 0
        fi
    done
    return 1
}

# Category heuristics (edit patterns to tune classification).
FLEX_PATTERNS=(
    "*/libflex/*"
    "*/Flex*"
    "*/Flex*/*"
    "*/ProbeList*"
    "*/FlexProbe*"
)
BULK_PATTERNS=(
    "*/libtrim/*"
    "*/libem/*"
    "*/CellRangerFormatter*"
    "*/TranscriptQuant*"
    "*/Quantifications*"
    "*/ReadAlign_quantTranscriptome*"
    "*/Transcriptome_quantAlign*"
    "*/Transcriptome_geneCounts*"
    "*/TranscriptomeFasta*"
    "*/Genome_genomeGenerate*"
)

count_code_local() {
    local label="$1"

    local files=()
    mapfile -d '' -t files < <(collect_source_files_local)
    local file_count=${#files[@]}
    local line_count=0

    if [[ "$file_count" -gt 0 ]]; then
        line_count=$(printf '%s\0' "${files[@]}" | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    fi

    echo "$label:"
    echo "  Files: $file_count"
    echo "  Lines: $line_count"
    echo ""

    eval "${label//[^a-zA-Z]/_}_files=$file_count"
    eval "${label//[^a-zA-Z]/_}_lines=$line_count"
}

count_code_ref() {
    local ref="$1"
    local label="$2"

    local files=()
    mapfile -t files < <(collect_source_files_ref "$ref")
    local file_count=${#files[@]}
    local line_count=0

    if [[ "$file_count" -gt 0 ]]; then
        while IFS= read -r f; do
            local n
            n=$(git -C "$REPO_ROOT" show "${ref}:${f}" 2>/dev/null | wc -l | awk '{print $1}')
            line_count=$((line_count + n))
        done < <(printf '%s\n' "${files[@]}")
    fi

    echo "$label ($ref):"
    echo "  Files: $file_count"
    echo "  Lines: $line_count"
    echo ""

    eval "${label//[^a-zA-Z]/_}_files=$file_count"
    eval "${label//[^a-zA-Z]/_}_lines=$line_count"
}

echo "--- Current STAR-suite Repository ---"
count_code_local "STAR_Suite"

echo "--- Category Breakdown ---"
flex_files=()
bulk_files=()
shared_files=()
all_files=()
mapfile -d '' -t all_files < <(collect_source_files_local)

for f in "${all_files[@]}"; do
    is_flex=0
    is_bulk=0
    if matches_any "$f" "${FLEX_PATTERNS[@]}"; then
        is_flex=1
    fi
    if matches_any "$f" "${BULK_PATTERNS[@]}"; then
        is_bulk=1
    fi

    if [[ "$is_flex" -eq 1 && "$is_bulk" -eq 0 ]]; then
        flex_files+=("$f")
    elif [[ "$is_bulk" -eq 1 && "$is_flex" -eq 0 ]]; then
        bulk_files+=("$f")
    else
        shared_files+=("$f")
    fi
done

count_category() {
    local label="$1"
    shift
    local files=("$@")
    local file_count=${#files[@]}
    local line_count=0
    if [[ "$file_count" -gt 0 ]]; then
        line_count=$(printf '%s\0' "${files[@]}" | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    fi
    echo "  $label:"
    echo "    Files: $file_count"
    echo "    Lines: $line_count"
}

count_category "Flex-specific" "${flex_files[@]}"
count_category "Bulk/PE-specific" "${bulk_files[@]}"
count_category "Shared/core (includes overlap)" "${shared_files[@]}"
echo ""

upstream_primary_source=""
if [[ "$SKIP_UPSTREAM" == "Yes" ]]; then
    echo "--- Upstream comparison skipped (--skip-upstream) ---"
    echo ""
else
    if ! git -C "$REPO_ROOT" rev-parse --verify "${UPSTREAM_REF}^{commit}" >/dev/null 2>&1; then
        echo "--- Upstream ref not available: ${UPSTREAM_REF} ---" >&2
        echo "SOLUTION: fetch upstream or pass --upstream-ref to an existing ref, or use --skip-upstream." >&2
        exit 1
    fi

    upstream_primary_source="$(primary_ref_source_prefix "$UPSTREAM_REF")"

    echo "--- Upstream STAR (from git ref) ---"
    count_code_ref "$UPSTREAM_REF" "Upstream_STAR"
fi

echo "=========================================="
echo "Comparison"
echo "=========================================="

total_files=()
mapfile -d '' -t total_files < <(collect_source_files_local)
suite_files=${#total_files[@]}
suite_lines=0
if [[ "$suite_files" -gt 0 ]]; then
    suite_lines=$(printf '%s\0' "${total_files[@]}" | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')
fi

upstream_files=0
upstream_lines=0
if [[ "$SKIP_UPSTREAM" != "Yes" ]]; then
    upstream_files=$(collect_source_files_ref "$UPSTREAM_REF" | wc -l | awk '{print $1}')
    upstream_lines=${Upstream_STAR_lines:-0}
fi

file_diff=$((suite_files - upstream_files))
line_diff=$((suite_lines - upstream_lines))

echo ""
echo "STAR-suite:    $suite_files files,  $suite_lines lines"
if [[ "$SKIP_UPSTREAM" == "Yes" ]]; then
    echo "Upstream STAR: (skipped)"
else
    echo "Upstream STAR: $upstream_files files,  $upstream_lines lines"
fi
echo ""
if [[ "$SKIP_UPSTREAM" == "Yes" ]]; then
    echo "Difference:    (skipped)"
else
    printf "Difference:    %+d files, %+d lines\n" "$file_diff" "$line_diff"
fi
echo ""

primary_local_source=""
if [[ -d "$REPO_ROOT/core/legacy/source" ]]; then
    primary_local_source="$REPO_ROOT/core/legacy/source"
elif [[ -d "$REPO_ROOT/source" ]]; then
    primary_local_source="$REPO_ROOT/source"
fi

echo "=========================================="
echo "New/Modified Files in primary source root"
echo "=========================================="
echo ""

if [[ -n "$primary_local_source" ]]; then
    echo "New C++ files added:"
    for f in $(find "$primary_local_source" -maxdepth 1 -type f \( -name "*.cpp" -o -name "*.h" \) -printf "%f\n" | sort); do
        if [[ "$SKIP_UPSTREAM" == "Yes" ]]; then
            lines=$(wc -l < "$primary_local_source/$f")
            echo "  + $f ($lines lines) [upstream comparison skipped]"
            continue
        fi
        if [[ -z "$upstream_primary_source" ]] || ! git -C "$REPO_ROOT" cat-file -e "${UPSTREAM_REF}:${upstream_primary_source}/${f}" 2>/dev/null; then
            lines=$(wc -l < "$primary_local_source/$f")
            echo "  + $f ($lines lines)"
        fi
    done
else
    echo "No primary source root found."
fi

echo ""
echo "Files in libflex/:"
if [[ -n "$primary_local_source" ]] && [[ -d "$primary_local_source/libflex" ]]; then
    libflex_files=$(find "$primary_local_source/libflex" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | wc -l)
    libflex_lines=$(find "$primary_local_source/libflex" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | xargs wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    echo "  $libflex_files files, $libflex_lines lines"
fi

echo ""
echo "Files in libtrim/:"
if [[ -n "$primary_local_source" ]] && [[ -d "$primary_local_source/libtrim" ]]; then
    libtrim_files=$(find "$primary_local_source/libtrim" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | wc -l)
    libtrim_lines=$(find "$primary_local_source/libtrim" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) -print0 | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    echo "  $libtrim_files files, $libtrim_lines lines"
fi

echo ""
echo "Files in libem/:"
if [[ -n "$primary_local_source" ]] && [[ -d "$primary_local_source/libem" ]]; then
    libem_files=$(find "$primary_local_source/libem" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | wc -l)
    libem_lines=$(find "$primary_local_source/libem" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) -print0 | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    echo "  $libem_files files, $libem_lines lines"
fi

echo ""
echo "Files in solo/:"
if [[ -n "$primary_local_source" ]] && [[ -d "$primary_local_source/solo" ]]; then
    solo_files=$(find "$primary_local_source/solo" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | wc -l)
    solo_lines=$(find "$primary_local_source/solo" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | xargs wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    echo "  $solo_files files, $solo_lines lines"
fi

echo ""
echo "Done."
