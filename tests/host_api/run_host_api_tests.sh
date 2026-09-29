#!/usr/bin/env bash
# Gate G-S2 for the STAR Suite host interface (docs/HOST_API.md):
#  1. A host with no callbacks (null hooks and empty hooks) reproduces STAR on
#     the synthetic scRNA fixture (tests/test_scrna_gex_counts.py, 8 UMI modes
#     plus the genome index); outputs are compared byte for byte except logs
#     and the command line in genomeParameters.txt.
#  2. A host with a dummy External domain, two threads taking and returning
#     permits while STAR maps, leaves Solo outputs identical to STAR's and the
#     permit exit invariant clean.
#  3. Unknown parameters are rejected without a host and recorded with one
#     (Log.out, effective command line, BAM @PG/@CO, --parametersFiles).
#  4. Failing lifecycle callbacks exit with STAR's exit codes and messages.
#  5. The permit unit tests, which use the public permit types and the
#     exported saturation controller, pass when linked with the library.
# Needs core/legacy/source/STAR and libstar_suite.a (make host-api-tests).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC="$ROOT/core/legacy/source"
STAR_BIN="${STAR_BIN:-$SRC/STAR}"
OUT="${HOST_API_TEST_OUT:-$(mktemp -d /tmp/star-host-api.XXXXXX)}"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
[[ -x "$STAR_BIN" ]] || { echo "ERROR: STAR not found at $STAR_BIN (make core)" >&2; exit 2; }
[[ -f "$SRC/libstar_suite.a" && -f "$SRC/libstar_suite.link" ]] || {
    echo "ERROR: build libstar_suite first (make star-host-lib)" >&2; exit 2; }

fail() { echo "FAIL: $*" >&2; exit 1; }

# genomeGenerate finds its expected-GC helper on PATH or relative to the
# executable; a host binary lives elsewhere, so put the helper on PATH for
# every binary compared here (docs/HOST_API.md).
GC_HELPER_DIR="$ROOT/core/features/vbem/tools/compute_expected_gc"
[[ -x "$GC_HELPER_DIR/compute_expected_gc" ]] || { echo "ERROR: build compute_expected_gc (make core)" >&2; exit 2; }
export PATH="$GC_HELPER_DIR:$PATH"

echo "Building the test host in $OUT"
# shellcheck disable=SC2046
g++ -std=c++11 -O2 -I"$SRC/host" "$ROOT/tests/host_api/star_host_test.cpp" \
    "$SRC/libstar_suite.a" $(cat "$SRC/libstar_suite.link") -o "$OUT/star_host_test"
HOST="$OUT/star_host_test"

# Byte comparison of two output trees, skipping logs, the gzip fixture inputs
# (their headers carry a timestamp) and the command-line header line of
# genomeParameters.txt; each tree's own path is replaced by a placeholder.
compare_trees() {
    python3 - "$1" "$2" <<'PY'
import sys
from pathlib import Path
a, b = Path(sys.argv[1]), Path(sys.argv[2])
skip = {"Log.out", "Log.progress.out", "Log.final.out", "Log.std.out",
        "run.log", "index.log", "host_events.tsv", "R1.gz", "R2.gz"}
def files(root):
    return {p.relative_to(root) for p in root.rglob("*") if p.is_file() and p.name not in skip}
fa, fb = files(a), files(b)
if fa != fb:
    sys.exit(f"file sets differ: only A {sorted(map(str, fa - fb))[:10]} only B {sorted(map(str, fb - fa))[:10]}")
diff = []
for rel in sorted(fa):
    x = (a / rel).read_bytes().replace(str(a).encode(), b"<RUN>")
    y = (b / rel).read_bytes().replace(str(b).encode(), b"<RUN>")
    if rel.name == "genomeParameters.txt":
        x, y = x.split(b"\n", 1)[1], y.split(b"\n", 1)[1]
    if x != y:
        diff.append(str(rel))
if diff:
    sys.exit(f"{len(diff)} files differ: {diff[:10]}")
print(f"identical: {len(fa)} files")
PY
}

echo "== 1. no-callback host reproduces STAR on the synthetic scRNA fixture"
for label in star null empty; do
    bin="$STAR_BIN"; mode=""
    [[ "$label" == star ]] || { bin="$HOST"; mode="$label"; }
    mkdir -p "$OUT/cwd_$label"
    ( cd "$OUT/cwd_$label" && STAR_HOST_TEST_MODE="$mode" \
        python3 "$ROOT/tests/test_scrna_gex_counts.py" --star "$bin" --outdir "$OUT/scrna_$label" \
        > "$OUT/scrna_$label.log" 2>&1 ) || { tail -20 "$OUT/scrna_$label.log"; fail "scRNA fixture with $label"; }
done
compare_trees "$OUT/scrna_star" "$OUT/scrna_null" || fail "null hooks differ from STAR"
compare_trees "$OUT/scrna_star" "$OUT/scrna_empty" || fail "empty hooks differ from STAR"

FIX="$OUT/scrna_star"
solo_args=(--genomeDir "$FIX/index" --readFilesIn "$FIX/R2.gz" "$FIX/R1.gz" --readFilesCommand zcat
    --soloType CB_UMI_Simple --soloCBwhitelist "$FIX/whitelist.txt" --soloCBlen 16 --soloUMIlen 12
    --soloUMIstart 17 --soloStrand Forward --soloFeatures Gene GeneFull
    --soloCellFilter CellRanger2.2 4 0.99 1 --clipAdapterType CellRanger4 --clip3pPolyG yes
    --outFilterScoreMinOverLread 0.8 --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR)
permit_args=(--runThreadN 2 --dynamicThreadInterface 1 --dynamicThreadTelemetry 1)

echo "== 2. dummy External domain on shared permits"
mkdir -p "$OUT/ext_star" "$OUT/ext_host"
"$STAR_BIN" "${solo_args[@]}" "${permit_args[@]}" --outSAMtype None \
    --outFileNamePrefix "$OUT/ext_star/" > "$OUT/ext_star/run.log" 2>&1 || fail "STAR permit run"
STAR_HOST_TEST_MODE=external "$HOST" "${solo_args[@]}" "${permit_args[@]}" --outSAMtype None \
    --outFileNamePrefix "$OUT/ext_host/" > "$OUT/ext_host/run.log" 2>&1 || {
    tail -20 "$OUT/ext_host/run.log"; fail "host external run"; }
compare_trees "$OUT/ext_star" "$OUT/ext_host" || fail "external-domain host outputs differ from STAR"
python3 - "$OUT" <<'PY'
import sys, re
from pathlib import Path
out = Path(sys.argv[1])
star_log = (out / "ext_star/Log.out").read_text()
host_log = (out / "ext_host/Log.out").read_text()
events = (out / "ext_host/host_events.tsv").read_text().splitlines()
kinds = [e.split("\t")[0] for e in events]
order = [k for k in kinds if k in ("preflight", "start", "extraPermitThreads", "initialFloors", "finish")]
assert order == ["preflight", "start", "extraPermitThreads", "initialFloors", "finish"], order
fin = dict(kv.split("=") for kv in next(e for e in events if e.startswith("finish\t")).split("\t")[1:])
assert fin["acquired"] == "400" and fin["externalAcquireCalls"] == "400", fin
assert fin["externalInUse"] == "0" and fin["available"] == fin["configured"] == "4", fin
assert "initialFloors\tconfigured=4\tin=0/0/0\tout=0/0/1" in events, events
assert "controllerLabels\tprobe-hosttest\texternal" in events, events
inv = re.search(r"Dynamic thread final invariant: .*", host_log).group(0)
assert "hosttestInUse=0" in inv and "inUse=0" in inv and "waiters=0" in inv, inv
assert "floors(map/feature/hosttest)=0/0/1" in host_log
assert "map permits=4 " in host_log and "[host test] preflight at " in host_log and "[host test] finish" in host_log
sinv = re.search(r"Dynamic thread final invariant: .*", star_log).group(0)
assert "externalInUse=0" in sinv, sinv
assert "floors(map/feature/external)=0/0/0" in star_log and "map permits=2 " in star_log
assert "FATAL" not in host_log and "FATAL" not in star_log
print("external domain: 400 host permits, invariant clean, labels and pool size as expected")
PY

echo "== 3. parameter pass-through"
mkdir -p "$OUT/par_star"
set +e
"$STAR_BIN" "${solo_args[@]}" --hostTestAlpha 7 --outFileNamePrefix "$OUT/par_star/" > "$OUT/par_star/run.log" 2>&1
rc=$?
set -e
[[ $rc -eq 102 ]] || fail "STAR accepted an unknown parameter (rc=$rc)"
grep -q 'unrecognized parameter name "hostTestAlpha" in input "Command-Line-Initial"' "$OUT/par_star/run.log" \
    || fail "STAR rejection message"
printf 'hostTestGamma 3\nhostTestDelta "a b"\n' > "$OUT/host_params.txt"
mkdir -p "$OUT/par_host"
STAR_HOST_TEST_MODE=params "$HOST" "${solo_args[@]}" --outSAMtype SAM --hostTestAlpha 7 \
    --hostTestBeta "two words" x --parametersFiles "$OUT/host_params.txt" --hostTestGamma 4 \
    --outFileNamePrefix "$OUT/par_host/" > "$OUT/par_host/run.log" 2>&1 || {
    tail -20 "$OUT/par_host/run.log"; fail "host parameter run"; }
python3 - "$OUT" <<'PY'
import sys
from pathlib import Path
out = Path(sys.argv[1]) / "par_host"
events = (out / "host_events.tsv").read_text().splitlines()
params = [e.split("\t")[1:] for e in events if e.startswith("parameter\t")]
names = [p[0] for p in params]
assert names == ["hostTestGamma", "hostTestDelta", "hostTestAlpha", "hostTestBeta", "hostTestGamma"], params
assert params[0][1:] == ["3", "5"] and params[1][1:] == ["a b", "5"], params
assert params[2][1:] == ["7", "2"] and params[3][1:] == ["two words|x", "2"] and params[4][1:] == ["4", "2"], params
pre = next(e for e in events if e.startswith("preflight\t"))
assert "commandLineFullHasHost=1" in pre, pre
log = (out / "Log.out").read_text()
full = log.split("##### Final effective command line:\n", 1)[1].splitlines()[0]
for token in ('--hostTestGamma 4', '--hostTestDelta "a b"', '--hostTestAlpha 7', '--hostTestBeta "two words"   x'):
    assert token in full, (token, full)
assert full.count("--hostTestGamma") == 1, full
sam = (out / "Aligned.out.sam").read_text().splitlines()
pg = next(l for l in sam if l.startswith("@PG"))
co = next(l for l in sam if l.startswith("@CO\tuser command line:"))
assert "--hostTestAlpha 7" in pg and "--hostTestGamma 4" in pg, pg
assert "--hostTestAlpha 7" in co and "--parametersFiles" in co, co
print("host parameters delivered in input order and recorded in Log.out, @PG and @CO")
PY
expect_fail() {  # rc, message, mode, extra args...
    local want_rc="$1" want_msg="$2" mode="$3"; shift 3
    local dir="$OUT/fail_${mode}_${want_rc}_$RANDOM"
    mkdir -p "$dir"
    set +e
    STAR_HOST_TEST_MODE="$mode" "$HOST" "${solo_args[@]}" "$@" --outFileNamePrefix "$dir/" > "$dir/run.log" 2>&1
    local rc=$?
    set -e
    [[ $rc -eq $want_rc ]] || fail "$mode $*: exit $rc, expected $want_rc"
    grep -qF "$want_msg" "$dir/run.log" || { cat "$dir/run.log"; fail "$mode $*: missing message: $want_msg"; }
}
expect_fail 102 'duplicate parameter "hostTestAlpha" in input "Command-Line"' params --hostTestAlpha 1 --hostTestAlpha 2
expect_fail 102 'invalid value for parameter "hostTestBad" in input "Command-Line": hostTestBad is always rejected' params --hostTestBad 1
expect_fail 102 'unrecognized parameter name "notAHostParameter" in input "Command-Line"' params --notAHostParameter 1

echo "== 4. failing lifecycle callbacks"
expect_fail 102 'EXITING because of fatal ERROR: host test preflight failure' fail-preflight
expect_fail 102 'SOLUTION: this failure is intentional' fail-preflight
expect_fail 103 'EXITING because of fatal ERROR: host test start failure' fail-start
expect_fail 103 'EXITING because of fatal ERROR: host test finish failure' fail-finish

echo "== 5. permit unit tests against libstar_suite"
# These tests include STAR-internal headers, so compile them with the HTSlib
# selection the library was built with (bundled libhts.a in the link file).
if grep -q '/core/legacy/source/htslib/libhts.a' "$SRC/libstar_suite.link"; then
    hts=(-DSTAR_EXTERNAL_HTSLIB=0 -I"$SRC/htslib")
else
    # shellcheck disable=SC2207
    hts=(-DSTAR_EXTERNAL_HTSLIB=1 $(pkg-config --cflags htslib 2>/dev/null))
fi
inc=("${hts[@]}" -I"$SRC" -I"$ROOT/flex/source" -I"$ROOT/flex/source/libflex"
     -I"$ROOT/core/features/libscrna/include" -I"$ROOT/core/features/vbem/source"
     -I"$ROOT/core/features/vbem/source/libem" -I"$ROOT/slam/source"
     -I"$ROOT/core/features/bamsort/source" -I"$ROOT/core/features/process_features/include")
for t in test_thread_control_permit_instrumentation test_bgzf_permit_hierarchy; do
    # shellcheck disable=SC2046
    g++ -std=c++11 -O2 -fopenmp "${inc[@]}" "$ROOT/tests/$t.cpp" "$SRC/libstar_suite.a" \
        $(cat "$SRC/libstar_suite.link") -o "$OUT/$t" || fail "build $t"
done
"$OUT/test_thread_control_permit_instrumentation" > "$OUT/permit_instrumentation.log" 2>&1 \
    || { tail -5 "$OUT/permit_instrumentation.log"; fail "test_thread_control_permit_instrumentation"; }
"$OUT/test_bgzf_permit_hierarchy" "$OUT/bgzf_hierarchy" > "$OUT/bgzf_hierarchy.log" 2>&1 \
    || { tail -5 "$OUT/bgzf_hierarchy.log"; fail "test_bgzf_permit_hierarchy"; }
echo "permit unit tests passed"

echo "PASS: STAR host API tests (G-S2); outputs in $OUT"
