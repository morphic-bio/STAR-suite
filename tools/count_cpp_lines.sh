#!/bin/bash
# Count C/C++ lines and files in STAR-Flex vs upstream STAR.
#
# Notes:
# - By default this compares against the upstream STAR tag `2.7.11b` using the local git object database
#   (no network required as long as the tag exists in this repo).
# - Counts include only compiled C/C++ sources under source/ and exclude htslib/opal and third-party folders.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
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

echo "=========================================="
echo "C/C++ Code Statistics: STAR-Flex vs Upstream"
echo "=========================================="
echo ""

# Collect source files (excluding third-party)
collect_source_files() {
    local dir="$1"
    find "$dir/source" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" \) \
        ! -path "*/htslib/*" \
        ! -path "*/opal/*" \
        ! -path "*third_party*" \
        -print0
}

# Collect source files from a git ref (excluding third-party)
collect_source_files_ref() {
    local ref="$1"
    git -C "$REPO_DIR" ls-tree -r --name-only "$ref" -- source \
        | grep -E '\.(c|cpp|h|hpp)$' \
        | grep -vE '^source/(htslib|opal)/' \
        | grep -v 'third_party' \
        || true
}

# Match file against any glob pattern in the list
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

# Count C/C++ lines and files by label
# Function to count lines and files
count_code() {
    local dir="$1"
    local label="$2"
    
    # Find all .c, .cpp, .h, .hpp files, excluding htslib and other third-party
    local files=()
    mapfile -d '' -t files < <(collect_source_files "$dir")
    local file_count=${#files[@]}
    local line_count=0
    
    if [[ "$file_count" -gt 0 ]]; then
        line_count=$(printf '%s\0' "${files[@]}" | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    fi
    
    echo "$label:"
    echo "  Files: $file_count"
    echo "  Lines: $line_count"
    echo ""
    
    # Return values for later use
    eval "${label//[^a-zA-Z]/_}_files=$file_count"
    eval "${label//[^a-zA-Z]/_}_lines=$line_count"
}

count_code_ref() {
    local ref="$1"
    local label="$2"

    mapfile -t files < <(collect_source_files_ref "$ref")
    local file_count=${#files[@]}
    local line_count=0

    if [[ "$file_count" -gt 0 ]]; then
        # Sum line counts without checking out a worktree.
        # (This is slower than wc -l on disk, but avoids network and extra temp directories.)
        while IFS= read -r f; do
            # git show prints the blob contents; wc -l counts newlines.
            local n
            n=$(git -C "$REPO_DIR" show "${ref}:${f}" 2>/dev/null | wc -l | awk '{print $1}')
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

# Count current repo
echo "--- Current STAR-Flex Repository ---"
count_code "$REPO_DIR" "STAR_Flex"

echo "--- STAR-Flex Category Breakdown ---"
flex_files=()
bulk_files=()
shared_files=()
all_files=()
mapfile -d '' -t all_files < <(collect_source_files "$REPO_DIR")

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

if [[ "$SKIP_UPSTREAM" == "Yes" ]]; then
    echo "--- Upstream comparison skipped (--skip-upstream) ---"
    echo ""
else
    if ! git -C "$REPO_DIR" rev-parse --verify "${UPSTREAM_REF}^{commit}" >/dev/null 2>&1; then
        echo "--- Upstream ref not available: ${UPSTREAM_REF} ---" >&2
        echo "SOLUTION: fetch upstream or pass --upstream-ref to an existing ref, or use --skip-upstream." >&2
        exit 1
    fi

    echo "--- Upstream STAR (from git ref) ---"
    count_code_ref "$UPSTREAM_REF" "Upstream_STAR"
fi

# Calculate differences
echo "=========================================="
echo "Comparison"
echo "=========================================="

flex_files=$(find "$REPO_DIR/source" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" \) \
    ! -path "*/htslib/*" ! -path "*/opal/*" ! -path "*third_party*" 2>/dev/null | wc -l | awk '{print $1}')
flex_lines=$(find "$REPO_DIR/source" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" \) \
    ! -path "*/htslib/*" ! -path "*/opal/*" ! -path "*third_party*" 2>/dev/null -print0 | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')

upstream_files=0
upstream_lines=0
if [[ "$SKIP_UPSTREAM" != "Yes" ]]; then
    upstream_files=$(collect_source_files_ref "$UPSTREAM_REF" | wc -l | awk '{print $1}')
    upstream_lines=${Upstream_STAR_lines:-0}
fi

file_diff=$((flex_files - upstream_files))
line_diff=$((flex_lines - upstream_lines))

echo ""
echo "STAR-Flex:     $flex_files files,  $flex_lines lines"
if [[ "$SKIP_UPSTREAM" == "Yes" ]]; then
    echo "Upstream STAR: (skipped)"
else
    echo "Upstream STAR: $upstream_files files,  $upstream_lines lines"
fi
echo ""
if [[ "$SKIP_UPSTREAM" == "Yes" ]]; then
    echo "Difference:    (skipped)"
else
    echo "Difference:    +$file_diff files, +$line_diff lines"
fi
echo ""

# Show new files in STAR-Flex
echo "=========================================="
echo "New/Modified Files in STAR-Flex (source/)"
echo "=========================================="
echo ""

# List files that are new or significantly different
echo "New C++ files added:"
for f in $(find "$REPO_DIR/source" -maxdepth 1 -type f \( -name "*.cpp" -o -name "*.h" \) -printf "%f\n" | sort); do
    if [[ "$SKIP_UPSTREAM" == "Yes" ]]; then
        lines=$(wc -l < "$REPO_DIR/source/$f")
        echo "  + $f ($lines lines) [upstream comparison skipped]"
        continue
    fi
    if ! git -C "$REPO_DIR" cat-file -e "${UPSTREAM_REF}:source/${f}" 2>/dev/null; then
        lines=$(wc -l < "$REPO_DIR/source/$f")
        echo "  + $f ($lines lines)"
    fi
done

echo ""
echo "Files in libflex/:"
if [[ -d "$REPO_DIR/source/libflex" ]]; then
    libflex_files=$(find "$REPO_DIR/source/libflex" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | wc -l)
    libflex_lines=$(find "$REPO_DIR/source/libflex" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | xargs wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    echo "  $libflex_files files, $libflex_lines lines"
fi

echo ""
echo "Files in libtrim/:"
if [[ -d "$REPO_DIR/source/libtrim" ]]; then
    libtrim_files=$(find "$REPO_DIR/source/libtrim" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | wc -l)
    libtrim_lines=$(find "$REPO_DIR/source/libtrim" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) -print0 | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    echo "  $libtrim_files files, $libtrim_lines lines"
fi

echo ""
echo "Files in libem/:"
if [[ -d "$REPO_DIR/source/libem" ]]; then
    libem_files=$(find "$REPO_DIR/source/libem" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | wc -l)
    libem_lines=$(find "$REPO_DIR/source/libem" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) -print0 | xargs -0 wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    echo "  $libem_files files, $libem_lines lines"
fi

echo ""
echo "Files in solo/:"
if [[ -d "$REPO_DIR/source/solo" ]]; then
    solo_files=$(find "$REPO_DIR/source/solo" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | wc -l)
    solo_lines=$(find "$REPO_DIR/source/solo" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.h" \) | xargs wc -l 2>/dev/null | tail -1 | awk '{print $1}')
    echo "  $solo_files files, $solo_lines lines"
fi

echo ""
echo "Done."
