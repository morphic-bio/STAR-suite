#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BINSEQ_BIN="${BINSEQ_PROBE_HARNESS_BIN:-$ROOT_DIR/core/legacy/source/binseq_probe_harness}"
OUT_ROOT="${OUT_ROOT:-/tmp/star_suite_binseq_upstream_fixture_smoke}"
UPSTREAM_REPO="${UPSTREAM_REPO:-https://github.com/ArcInstitute/binseq.git}"

skip() {
    echo "SKIP: $*"
    exit 0
}

if [[ -n "${BQTOOLS:-}" ]]; then
    if [[ -x "$BQTOOLS" ]]; then
        BQTOOLS_BIN="$BQTOOLS"
    elif command -v "$BQTOOLS" >/dev/null 2>&1; then
        BQTOOLS_BIN="$(command -v "$BQTOOLS")"
    else
        skip "BQTOOLS is set but not executable or resolvable: $BQTOOLS"
    fi
else
    command -v bqtools >/dev/null 2>&1 || skip "bqtools not found in PATH; set BQTOOLS=/path/to/bqtools"
    BQTOOLS_BIN="$(command -v bqtools)"
fi

command -v git >/dev/null 2>&1 || skip "git not found in PATH"

if [[ ! -x "$BINSEQ_BIN" ]]; then
    make -C "$ROOT_DIR/core/legacy/source" binseq-probe-harness
fi

rm -rf "$OUT_ROOT"
mkdir -p "$OUT_ROOT"/{repo,info,dumps,probe}

if ! git clone --depth 1 "$UPSTREAM_REPO" "$OUT_ROOT/repo/binseq" \
    > "$OUT_ROOT/git_clone.stdout" 2> "$OUT_ROOT/git_clone.stderr"; then
    skip "could not clone $UPSTREAM_REPO; see $OUT_ROOT/git_clone.stderr"
fi

DATA_DIR="$OUT_ROOT/repo/binseq/data"
CBQ="$DATA_DIR/subset.cbq"
VBQ="$DATA_DIR/subset.vbq"
BQ="$DATA_DIR/subset.bq"

for path in "$CBQ" "$VBQ" "$BQ" "$DATA_DIR/subset_R1.fastq.gz" "$DATA_DIR/subset_R2.fastq.gz"; do
    [[ -s "$path" ]] || {
        echo "ERROR: expected upstream fixture missing or empty: $path" >&2
        exit 1
    }
done

"$BQTOOLS_BIN" info --json "$CBQ" > "$OUT_ROOT/info/subset.cbq.json"
"$BQTOOLS_BIN" info --json "$VBQ" > "$OUT_ROOT/info/subset.vbq.json"
"$BQTOOLS_BIN" info --json "$BQ" > "$OUT_ROOT/info/subset.bq.json"

records="$("$BQTOOLS_BIN" info --num "$CBQ" | awk 'NR == 1 { gsub("_", "", $1); print $1 }')"
if [[ ! "$records" =~ ^[0-9]+$ || "$records" -le 0 ]]; then
    echo "ERROR: could not determine positive record count from $CBQ" >&2
    exit 1
fi

"$BINSEQ_BIN" --readNameSeparator space \
    --mateCount 2 \
    --bqtools "$BQTOOLS_BIN" \
    --workDir "$OUT_ROOT/probe/subset_cbq" \
    --readFilesIn "$CBQ" \
    > "$OUT_ROOT/dumps/subset_cbq_probe.tsv"

line_count="$(wc -l < "$OUT_ROOT/dumps/subset_cbq_probe.tsv" | tr -d '[:space:]')"
expected_lines=$((records * 2))
if [[ "$line_count" -ne "$expected_lines" ]]; then
    echo "ERROR: expected $expected_lines probe rows for $records paired records, observed $line_count" >&2
    exit 1
fi

awk -F'\t' -v expected="$records" '
    $5 == 1 { mate1 += 1 }
    $5 == 2 { mate2 += 1 }
    END {
        if (mate1 != expected || mate2 != expected) {
            printf("ERROR: mate row counts mismatch: mate1=%d mate2=%d expected=%d\n", mate1, mate2, expected) > "/dev/stderr"
            exit 1
        }
    }
' "$OUT_ROOT/dumps/subset_cbq_probe.tsv"

python3 - "$DATA_DIR/subset_R1.fastq.gz" \
    "$DATA_DIR/subset_R2.fastq.gz" \
    "$OUT_ROOT/probe/subset_cbq/lane0_R1.fq" \
    "$OUT_ROOT/probe/subset_cbq/lane0_R2.fq" <<'PY'
import gzip
import sys

fastq_r1, fastq_r2, cbq_r1, cbq_r2 = sys.argv[1:5]

def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")

def read_fastq(path):
    records = []
    with open_text(path) as handle:
        while True:
            header = handle.readline().rstrip("\n")
            if not header:
                break
            sequence = handle.readline().rstrip("\n")
            plus = handle.readline().rstrip("\n")
            quality = handle.readline().rstrip("\n")
            if not plus.startswith("+"):
                raise SystemExit(f"bad FASTQ record in {path}: {header}")
            name = header[1:].split()[0]
            records.append((name, sequence, quality))
    return records

fq1 = read_fastq(fastq_r1)
fq2 = read_fastq(fastq_r2)
cb1 = read_fastq(cbq_r1)
cb2 = read_fastq(cbq_r2)

if len({len(fq1), len(fq2), len(cb1), len(cb2)}) != 1:
    raise SystemExit(
        "record count mismatch: "
        f"fastq_r1={len(fq1)} fastq_r2={len(fq2)} cbq_r1={len(cb1)} cbq_r2={len(cb2)}"
    )

for label, expected, observed in (("R1", fq1, cb1), ("R2", fq2, cb2)):
    expected_set = set(expected)
    observed_set = set(observed)
    if expected_set != observed_set:
        only_expected = sorted(expected_set - observed_set)[:3]
        only_observed = sorted(observed_set - expected_set)[:3]
        raise SystemExit(
            f"{label} order-independent FASTQ parity failed: "
            f"expected_only={only_expected} observed_only={only_observed}"
        )

expected_pairs = set(zip(fq1, fq2))
observed_pairs = set(zip(cb1, cb2))
if expected_pairs != observed_pairs:
    only_expected = sorted(expected_pairs - observed_pairs)[:3]
    only_observed = sorted(observed_pairs - expected_pairs)[:3]
    raise SystemExit(
        "paired order-independent FASTQ parity failed: "
        f"expected_only={only_expected} observed_only={only_observed}"
    )

print(f"order-independent FASTQ parity OK for {len(fq1)} paired records")
PY

cat > "$OUT_ROOT/NOTES.txt" <<EOF
Source: $UPSTREAM_REPO
Fixture: data/subset.cbq
Records: $records

The upstream repository also contains subset_R1.fastq.gz and subset_R2.fastq.gz.
Those files contain the same reads as subset.cbq but not in the same order, so
this smoke validates order-independent FASTQ parity after decoding subset.cbq
through the STAR-suite BINSEQ probe module.
EOF

echo "PASS: upstream BINSEQ fixture smoke completed at $OUT_ROOT"
