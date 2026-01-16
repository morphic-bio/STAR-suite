#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

AUTO_DIR="${AUTO_DIR:-}"
BASE_DIR="${BASE_DIR:-}"
OUT_DIR="${OUT_DIR:-}"
RAW_R1="${RAW_R1:-}"
RAW_R2="${RAW_R2:-}"
TRIMMED_R1="${TRIMMED_R1:-}"
TRIMMED_R2="${TRIMMED_R2:-}"
TRIM_SAMPLE="${TRIM_SAMPLE:-100000}"
SKIP_TRIM="${SKIP_TRIM:-0}"
SKIP_BAM="${SKIP_BAM:-0}"
AUTO_LABEL="${AUTO_LABEL:-auto}"
BASE_LABEL="${BASE_LABEL:-base}"
TRIMVALIDATE_BIN="${TRIMVALIDATE_BIN:-$REPO_ROOT/tools/trimvalidate/trimvalidate}"
IDXSTATS_SORT_UNSORTED="${IDXSTATS_SORT_UNSORTED:-0}"
IDXSTATS_SORT_THREADS="${IDXSTATS_SORT_THREADS:-4}"
TXIMPORT_CLI="${TXIMPORT_CLI:-$REPO_ROOT/tools/tximport_compat/tximport_compat}"
TX2GENE="${TX2GENE:-}"
AUTO_SALMON_DIR="${AUTO_SALMON_DIR:-}"
BASE_SALMON_DIR="${BASE_SALMON_DIR:-}"
MORPHIC_SALMON_DIR="${MORPHIC_SALMON_DIR:-}"
MORPHIC_LABEL="${MORPHIC_LABEL:-morphic}"
SKIP_TXIMPORT="${SKIP_TXIMPORT:-0}"

usage() {
    cat <<'USAGE'
Usage: compare_run_parity.sh --auto-dir DIR --base-dir DIR --out-dir DIR [options]

Options:
  --auto-dir DIR          Run directory for autoindex/trimCutadapt run
  --base-dir DIR          Run directory for Trim Galore baseline run
  --out-dir DIR           Output report directory
  --raw-r1 PATH           Raw R1 FASTQ (gz or plain) for trim parity
  --raw-r2 PATH           Raw R2 FASTQ (gz or plain) for trim parity
  --trimmed-r1 PATH       Trim Galore R1 output (gz or plain) for trim parity
  --trimmed-r2 PATH       Trim Galore R2 output (gz or plain) for trim parity
  --trim-sample N         Number of read pairs to sample for trim parity (default: 100000)
                          Use 0 to compare full files (may be large)
  --skip-trim             Skip trim parity checks
  --skip-bam              Skip BAM-level stats comparisons
  --auto-label LABEL      Label for auto run (default: auto)
  --base-label LABEL      Label for base run (default: base)
  --trimvalidate PATH     Path to trimvalidate binary
  --idxstats-sort-unsorted  Sort unsorted BAMs to temp for idxstats
  --tx2gene PATH          tx2gene mapping for tximport_compat
  --tximport-cli PATH     Path to tximport_compat binary
  --auto-salmon-dir DIR   Salmon output dir for auto run
  --base-salmon-dir DIR   Salmon output dir for base run
  --morphic-salmon-dir DIR  Salmon output dir for morphic run
  --morphic-label LABEL   Label for morphic run (default: morphic)
  --skip-tximport         Skip tximport comparisons
  -h, --help              Show this help

Environment variables can be used instead of flags:
  AUTO_DIR, BASE_DIR, OUT_DIR, RAW_R1, RAW_R2, TRIMMED_R1, TRIMMED_R2,
  TRIM_SAMPLE, SKIP_TRIM, SKIP_BAM, AUTO_LABEL, BASE_LABEL, TRIMVALIDATE_BIN,
  IDXSTATS_SORT_UNSORTED, IDXSTATS_SORT_THREADS, TXIMPORT_CLI, TX2GENE,
  AUTO_SALMON_DIR, BASE_SALMON_DIR, MORPHIC_SALMON_DIR, MORPHIC_LABEL, SKIP_TXIMPORT
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --auto-dir) AUTO_DIR="$2"; shift 2 ;;
        --base-dir) BASE_DIR="$2"; shift 2 ;;
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --raw-r1) RAW_R1="$2"; shift 2 ;;
        --raw-r2) RAW_R2="$2"; shift 2 ;;
        --trimmed-r1) TRIMMED_R1="$2"; shift 2 ;;
        --trimmed-r2) TRIMMED_R2="$2"; shift 2 ;;
        --trim-sample) TRIM_SAMPLE="$2"; shift 2 ;;
        --skip-trim) SKIP_TRIM=1; shift ;;
        --skip-bam) SKIP_BAM=1; shift ;;
        --auto-label) AUTO_LABEL="$2"; shift 2 ;;
        --base-label) BASE_LABEL="$2"; shift 2 ;;
        --trimvalidate) TRIMVALIDATE_BIN="$2"; shift 2 ;;
        --idxstats-sort-unsorted) IDXSTATS_SORT_UNSORTED=1; shift ;;
        --tx2gene) TX2GENE="$2"; shift 2 ;;
        --tximport-cli) TXIMPORT_CLI="$2"; shift 2 ;;
        --auto-salmon-dir) AUTO_SALMON_DIR="$2"; shift 2 ;;
        --base-salmon-dir) BASE_SALMON_DIR="$2"; shift 2 ;;
        --morphic-salmon-dir) MORPHIC_SALMON_DIR="$2"; shift 2 ;;
        --morphic-label) MORPHIC_LABEL="$2"; shift 2 ;;
        --skip-tximport) SKIP_TXIMPORT=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$AUTO_DIR" || -z "$BASE_DIR" || -z "$OUT_DIR" ]]; then
    echo "Error: --auto-dir, --base-dir, and --out-dir are required." >&2
    usage
    exit 1
fi

mkdir -p "$OUT_DIR"
REPORT="$OUT_DIR/report.txt"
exec > >(tee "$REPORT") 2>&1

echo "== STAR-Flex run parity report =="
echo "Auto dir: $AUTO_DIR"
echo "Base dir: $BASE_DIR"
echo "Out dir : $OUT_DIR"
echo ""

extract_command() {
    local log="$1"
    awk '/^##### Command Line:/{getline; print; exit}' "$log"
}

extract_final_params() {
    local log="$1"
    awk '
        /^##### Final user re-defined parameters/ {in_block=1; next}
        in_block && /^-+$/ {exit}
        in_block {print}
    ' "$log"
}

save_log_sections() {
    local label="$1"
    local dir="$2"
    local log="$dir/Log.out"
    if [[ ! -f "$log" ]]; then
        echo "WARN: Missing Log.out in $dir"
        return
    fi
    extract_command "$log" > "$OUT_DIR/${label}.command.txt" || true
    extract_final_params "$log" > "$OUT_DIR/${label}.params.txt" || true
    head -n 30 "$log" > "$OUT_DIR/${label}.Log.head.txt" || true
    if [[ -f "$dir/Log.final.out" ]]; then
        cp "$dir/Log.final.out" "$OUT_DIR/${label}.Log.final.out"
    else
        echo "WARN: Missing Log.final.out in $dir"
    fi
}

compare_quant_corr() {
    local label="$1"
    local file_a="$2"
    local file_b="$3"
    python3 << PYEOF
import math

def read_quant(path):
    data = {}
    with open(path) as f:
        header = f.readline().strip().split("\\t")
        col_idx = {name: i for i, name in enumerate(header)}
        for line in f:
            fields = line.rstrip("\\n").split("\\t")
            name = fields[0]
            data[name] = {
                "TPM": float(fields[col_idx["TPM"]]),
                "NumReads": float(fields[col_idx["NumReads"]]),
            }
    return data

def pearson(x, y):
    n = len(x)
    if n < 2:
        return float("nan")
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    num = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    den_x = math.sqrt(sum((a - mean_x) ** 2 for a in x))
    den_y = math.sqrt(sum((b - mean_y) ** 2 for b in y))
    if den_x == 0 or den_y == 0:
        return float("nan")
    return num / (den_x * den_y)

def rankdata(vals):
    sorted_vals = sorted((v, i) for i, v in enumerate(vals))
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(sorted_vals):
        j = i
        while j < len(sorted_vals) and sorted_vals[j][0] == sorted_vals[i][0]:
            j += 1
        rank = (i + j + 1) / 2.0
        for k in range(i, j):
            ranks[sorted_vals[k][1]] = rank
        i = j
    return ranks

def spearman(x, y):
    if len(x) < 2:
        return float("nan")
    rx = rankdata(x)
    ry = rankdata(y)
    return pearson(rx, ry)

data_a = read_quant("$file_a")
data_b = read_quant("$file_b")
keys = sorted(set(data_a.keys()) & set(data_b.keys()))

if not keys:
    print("[$label] No shared genes to compare.")
    raise SystemExit(0)

xa = [data_a[k]["NumReads"] for k in keys]
xb = [data_b[k]["NumReads"] for k in keys]
ta = [data_a[k]["TPM"] for k in keys]
tb = [data_b[k]["TPM"] for k in keys]

print(f"[$label] Shared genes: {len(keys)}")
print(f"[$label] NumReads Spearman: {spearman(xa, xb):.6f}")
print(f"[$label] NumReads Pearson:  {pearson(xa, xb):.6f}")
print(f"[$label] TPM Spearman:      {spearman(ta, tb):.6f}")
print(f"[$label] TPM Pearson:       {pearson(ta, tb):.6f}")
PYEOF
}

run_tximport() {
    local quant="$1"
    local out="$2"
    if [ ! -x "$TXIMPORT_CLI" ]; then
        echo "WARN: tximport_compat not found: $TXIMPORT_CLI"
        return 1
    fi
    if [ -z "$TX2GENE" ]; then
        echo "WARN: TX2GENE not set; skipping tximport_compat for $quant"
        return 1
    fi
    if [ ! -f "$quant" ]; then
        echo "WARN: quant.sf not found: $quant"
        return 1
    fi
    "$TXIMPORT_CLI" --quant "$quant" --tx2gene "$TX2GENE" --output "$out" --stats
    return 0
}

echo "== Extracting Log.out sections =="
save_log_sections "$AUTO_LABEL" "$AUTO_DIR"
save_log_sections "$BASE_LABEL" "$BASE_DIR"

if [[ -f "$OUT_DIR/${AUTO_LABEL}.command.txt" && -f "$OUT_DIR/${BASE_LABEL}.command.txt" ]]; then
    diff -u "$OUT_DIR/${BASE_LABEL}.command.txt" "$OUT_DIR/${AUTO_LABEL}.command.txt" \
        > "$OUT_DIR/command.diff" || true
    echo "Command diff: $OUT_DIR/command.diff"
fi

if [[ -f "$OUT_DIR/${AUTO_LABEL}.params.txt" && -f "$OUT_DIR/${BASE_LABEL}.params.txt" ]]; then
    diff -u "$OUT_DIR/${BASE_LABEL}.params.txt" "$OUT_DIR/${AUTO_LABEL}.params.txt" \
        > "$OUT_DIR/params.diff" || true
    echo "Params diff:  $OUT_DIR/params.diff"
fi

if [[ -f "$OUT_DIR/${AUTO_LABEL}.Log.final.out" && -f "$OUT_DIR/${BASE_LABEL}.Log.final.out" ]]; then
    diff -u "$OUT_DIR/${BASE_LABEL}.Log.final.out" "$OUT_DIR/${AUTO_LABEL}.Log.final.out" \
        > "$OUT_DIR/Log.final.diff" || true
    echo "Log.final diff: $OUT_DIR/Log.final.diff"
fi

echo ""
echo "== BAM-level comparisons =="
if [[ "$SKIP_BAM" -eq 1 ]]; then
    echo "Skipping BAM comparisons (--skip-bam)."
else
    if ! command -v samtools >/dev/null 2>&1; then
        echo "WARN: samtools not found; skipping BAM comparisons."
    else
        pick_bam() {
            local dir="$1"
            if [[ -f "$dir/Aligned.sortedByCoord.out.bam" ]]; then
                echo "$dir/Aligned.sortedByCoord.out.bam"
            elif [[ -f "$dir/Aligned.out.bam" ]]; then
                echo "$dir/Aligned.out.bam"
            else
                echo ""
            fi
        }
        AUTO_BAM="$(pick_bam "$AUTO_DIR")"
        BASE_BAM="$(pick_bam "$BASE_DIR")"
        if [[ -z "$AUTO_BAM" || -z "$BASE_BAM" ]]; then
            echo "WARN: Could not locate BAMs in one or both run dirs."
        else
            echo "Auto BAM: $AUTO_BAM"
            echo "Base BAM: $BASE_BAM"
            samtools flagstat "$AUTO_BAM" > "$OUT_DIR/${AUTO_LABEL}.flagstat"
            samtools flagstat "$BASE_BAM" > "$OUT_DIR/${BASE_LABEL}.flagstat"
            diff -u "$OUT_DIR/${BASE_LABEL}.flagstat" "$OUT_DIR/${AUTO_LABEL}.flagstat" \
                > "$OUT_DIR/flagstat.diff" || true
            echo "flagstat diff: $OUT_DIR/flagstat.diff"

            bam_sort_order() {
                local bam="$1"
                samtools view -H "$bam" 2>/dev/null | awk -F'\t' '
                    /^@HD/ {
                        for (i=1;i<=NF;i++) {
                            if ($i ~ /^SO:/) { sub(/^SO:/,"",$i); print $i; exit }
                        }
                    }'
            }

            run_idxstats() {
                local label="$1"
                local bam="$2"
                local out="$OUT_DIR/${label}.idxstats"
                local so
                so="$(bam_sort_order "$bam")"
                if [[ "$so" == "coordinate" ]]; then
                    if samtools idxstats "$bam" > "$out"; then
                        echo "$out"
                    else
                        echo "WARN: idxstats failed for $bam"
                        rm -f "$out"
                    fi
                    return
                fi
                if [[ "$IDXSTATS_SORT_UNSORTED" -eq 1 ]]; then
                    local tmp_dir
                    tmp_dir="$(mktemp -d)"
                    local tmp_bam="$tmp_dir/coord.bam"
                    echo "Sorting unsorted BAM for idxstats ($label, SO=${so:-unknown})..."
                    if samtools sort -@ "$IDXSTATS_SORT_THREADS" -o "$tmp_bam" "$bam"; then
                        if samtools idxstats "$tmp_bam" > "$out"; then
                            rm -rf "$tmp_dir"
                            echo "$out"
                            return
                        fi
                    fi
                    echo "WARN: idxstats failed for sorted temp BAM ($label)"
                    rm -rf "$tmp_dir"
                else
                    echo "Skipping idxstats for unsorted BAM ($label, SO=${so:-unknown})."
                    echo "  Set --idxstats-sort-unsorted to enable sorting."
                fi
            }

            AUTO_IDXSTATS="$(run_idxstats "$AUTO_LABEL" "$AUTO_BAM" || true)"
            BASE_IDXSTATS="$(run_idxstats "$BASE_LABEL" "$BASE_BAM" || true)"
            if [[ -n "$AUTO_IDXSTATS" && -n "$BASE_IDXSTATS" ]]; then
                diff -u "$BASE_IDXSTATS" "$AUTO_IDXSTATS" > "$OUT_DIR/idxstats.diff" || true
                echo "idxstats diff: $OUT_DIR/idxstats.diff"
            else
                echo "idxstats diff: skipped (missing one or both idxstats outputs)"
            fi
        fi
    fi
fi

echo ""
echo "== Tximport comparisons =="
if [[ "$SKIP_TXIMPORT" -eq 1 ]]; then
    echo "Skipping tximport comparisons (--skip-tximport)."
else
    TXIMPORT_DIR="$OUT_DIR/tximport"
    mkdir -p "$TXIMPORT_DIR"

    auto_star_simple="$AUTO_DIR/quant.genes.sf"
    auto_star_tx="$AUTO_DIR/quant.genes.tximport.sf"
    base_star_simple="$BASE_DIR/quant.genes.sf"
    base_star_tx="$BASE_DIR/quant.genes.tximport.sf"

    if [[ -f "$auto_star_simple" && -f "$auto_star_tx" ]]; then
        compare_quant_corr "Auto STAR legacy vs tximport" "$auto_star_simple" "$auto_star_tx"
    else
        echo "WARN: Auto STAR gene files missing for legacy vs tximport comparison."
    fi

    if [[ -f "$base_star_simple" && -f "$base_star_tx" ]]; then
        compare_quant_corr "Base STAR legacy vs tximport" "$base_star_simple" "$base_star_tx"
    else
        echo "WARN: Base STAR gene files missing for legacy vs tximport comparison."
    fi

    if [[ -f "$auto_star_tx" && -f "$base_star_tx" ]]; then
        compare_quant_corr "STAR tximport auto vs base" "$auto_star_tx" "$base_star_tx"
    fi

    if [[ -n "$AUTO_SALMON_DIR" ]]; then
        auto_salmon_quant="$AUTO_SALMON_DIR/quant.sf"
        auto_salmon_simple="$AUTO_SALMON_DIR/quant.genes.sf"
        auto_salmon_tx="$TXIMPORT_DIR/${AUTO_LABEL}_salmon.tximport.sf"
        if run_tximport "$auto_salmon_quant" "$auto_salmon_tx"; then
            if [[ -f "$auto_salmon_simple" ]]; then
                compare_quant_corr "Auto Salmon simple vs tximport" "$auto_salmon_simple" "$auto_salmon_tx"
            else
                echo "WARN: Auto Salmon quant.genes.sf missing; simple vs tximport skipped."
            fi
            if [[ -f "$auto_star_tx" ]]; then
                compare_quant_corr "Auto STAR tximport vs Salmon tximport" "$auto_star_tx" "$auto_salmon_tx"
            fi
        fi
    fi

    if [[ -n "$BASE_SALMON_DIR" ]]; then
        base_salmon_quant="$BASE_SALMON_DIR/quant.sf"
        base_salmon_simple="$BASE_SALMON_DIR/quant.genes.sf"
        base_salmon_tx="$TXIMPORT_DIR/${BASE_LABEL}_salmon.tximport.sf"
        if run_tximport "$base_salmon_quant" "$base_salmon_tx"; then
            if [[ -f "$base_salmon_simple" ]]; then
                compare_quant_corr "Base Salmon simple vs tximport" "$base_salmon_simple" "$base_salmon_tx"
            else
                echo "WARN: Base Salmon quant.genes.sf missing; simple vs tximport skipped."
            fi
            if [[ -f "$base_star_tx" ]]; then
                compare_quant_corr "Base STAR tximport vs Salmon tximport" "$base_star_tx" "$base_salmon_tx"
            fi
        fi
    fi

    if [[ -n "$MORPHIC_SALMON_DIR" ]]; then
        morphic_salmon_quant="$MORPHIC_SALMON_DIR/quant.sf"
        morphic_salmon_simple="$MORPHIC_SALMON_DIR/quant.genes.sf"
        morphic_salmon_tx="$TXIMPORT_DIR/${MORPHIC_LABEL}_salmon.tximport.sf"
        if run_tximport "$morphic_salmon_quant" "$morphic_salmon_tx"; then
            if [[ -f "$morphic_salmon_simple" ]]; then
                compare_quant_corr "$MORPHIC_LABEL Salmon simple vs tximport" "$morphic_salmon_simple" "$morphic_salmon_tx"
            else
                echo "WARN: $MORPHIC_LABEL Salmon quant.genes.sf missing; simple vs tximport skipped."
            fi
            if [[ -f "$auto_star_tx" ]]; then
                compare_quant_corr "$MORPHIC_LABEL Salmon tximport vs Auto STAR tximport" "$morphic_salmon_tx" "$auto_star_tx"
            elif [[ -f "$base_star_tx" ]]; then
                compare_quant_corr "$MORPHIC_LABEL Salmon tximport vs Base STAR tximport" "$morphic_salmon_tx" "$base_star_tx"
            fi
        fi
    fi
fi

echo ""
echo "== Trim parity checks =="
if [[ "$SKIP_TRIM" -eq 1 ]]; then
    echo "Skipping trim parity checks (--skip-trim)."
else
    if [[ -z "$RAW_R1" || -z "$RAW_R2" || -z "$TRIMMED_R1" || -z "$TRIMMED_R2" ]]; then
        echo "WARN: Missing raw/trimmed inputs; skipping trim parity checks."
    elif [[ ! -x "$TRIMVALIDATE_BIN" ]]; then
        echo "WARN: trimvalidate not found or not executable: $TRIMVALIDATE_BIN"
    else
        TRIM_DIR="$OUT_DIR/trim_compare"
        mkdir -p "$TRIM_DIR"
        TMPDIR="$(mktemp -d)"
        trap "rm -rf \"$TMPDIR\"" EXIT

        materialize_fastq() {
            local src="$1"
            local dst="$2"
            local lines="${3:-}"
            if [[ "$src" == *.gz ]]; then
                if [[ -n "$lines" ]]; then
                    zcat "$src" | head -n "$lines" > "$dst"
                else
                    zcat "$src" > "$dst"
                fi
            else
                if [[ -n "$lines" ]]; then
                    head -n "$lines" "$src" > "$dst"
                else
                    cp "$src" "$dst"
                fi
            fi
        }

        if [[ "$TRIM_SAMPLE" -lt 0 ]]; then
            echo "Error: --trim-sample must be >= 0"
            exit 1
        fi

        if [[ "$TRIM_SAMPLE" -eq 0 ]]; then
            echo "Running full-file trim parity (may be large)..."
            RAW_R1_MAT="$TMPDIR/raw_R1.fastq"
            RAW_R2_MAT="$TMPDIR/raw_R2.fastq"
            materialize_fastq "$RAW_R1" "$RAW_R1_MAT"
            materialize_fastq "$RAW_R2" "$RAW_R2_MAT"
        else
            echo "Sampling first $TRIM_SAMPLE read pairs for trim parity..."
            LINES=$((TRIM_SAMPLE * 4))
            RAW_R1_MAT="$TMPDIR/raw_R1.sample.fastq"
            RAW_R2_MAT="$TMPDIR/raw_R2.sample.fastq"
            materialize_fastq "$RAW_R1" "$RAW_R1_MAT" "$LINES"
            materialize_fastq "$RAW_R2" "$RAW_R2_MAT" "$LINES"
        fi

        TRIM_OUT_R1="$TRIM_DIR/trimvalidate_R1.fastq"
        TRIM_OUT_R2="$TRIM_DIR/trimvalidate_R2.fastq"
        "$TRIMVALIDATE_BIN" -1 "$RAW_R1_MAT" -2 "$RAW_R2_MAT" \
            -o1 "$TRIM_OUT_R1" -o2 "$TRIM_OUT_R2" \
            --quality 20 --length 20 \
            > "$TRIM_DIR/trimvalidate.log" 2>&1

        OUT_LINES_R1=$(wc -l < "$TRIM_OUT_R1" || echo 0)
        OUT_READS=$((OUT_LINES_R1 / 4))
        if [[ "$OUT_READS" -eq 0 ]]; then
            echo "WARN: trimvalidate produced 0 reads; check $TRIM_DIR/trimvalidate.log"
            exit 1
        fi

        BASELINES_LINES=$((OUT_READS * 4))
        BASE_R1_SAMPLE="$TRIM_DIR/baseline_R1.fastq"
        BASE_R2_SAMPLE="$TRIM_DIR/baseline_R2.fastq"
        materialize_fastq "$TRIMMED_R1" "$BASE_R1_SAMPLE" "$BASELINES_LINES"
        materialize_fastq "$TRIMMED_R2" "$BASE_R2_SAMPLE" "$BASELINES_LINES"

        diff -u "$BASE_R1_SAMPLE" "$TRIM_OUT_R1" > "$TRIM_DIR/trim_R1.diff" || true
        diff -u "$BASE_R2_SAMPLE" "$TRIM_OUT_R2" > "$TRIM_DIR/trim_R2.diff" || true

        echo "Trim R1 diff: $TRIM_DIR/trim_R1.diff"
        echo "Trim R2 diff: $TRIM_DIR/trim_R2.diff"
    fi
fi

echo ""
echo "Done. Report: $REPORT"
