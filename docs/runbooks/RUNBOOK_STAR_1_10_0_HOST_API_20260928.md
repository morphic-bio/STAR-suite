# Runbook: STAR Suite 1.10.0 host interface and Chromap removal (2026-09-28)

## Goal and scope

Carry out the STAR Suite side, steps S1-S4, of the single-binary design:
`/mnt/pikachu/multiomics-suite-single-binary-20260928/docs/design/SINGLE_BINARY_OWNERSHIP_20260928.md`
(latest copy, with the Step 0 corrections and the author decisions of
28 September). STAR Suite 1.10.0 is a minor release whose focus is removing
every dependency on other suites: STAR stops hosting Chromap and gains a small
generic host interface that knows nothing about Chromap or ATAC. The result is
a local `v1.10.0-rc1` tag on `dev-release-v1.10.0`; it stays an rc until
Multiomics Suite has exercised the host interface.

Author decisions that apply here: the third permit domain `ATAC` becomes a
neutral `EXTERNAL` with a host label; the STAR parser passes unknown
parameters to the host; `make core` uses bundled HTSlib with `HTSLIB=external`
optional; the CAT-ATAC guide layout stays in STAR; libscrna's multiome tools
move to Multiomics; `CbqChromapAdapter` with its harness and smoke test moves
(recorded for Multiomics, not deleted); `dev-release-v1.9.6` stays out; a new
feature is a minor release.

Worktree: `/mnt/pikachu/STAR-suite-v1100-20260928`, branch
`dev-release-v1.10.0` from `master` (`4824548`, v1.9.5.a).

The published `v1.9.5.a` tag is `f0d9f27`. Baseline `4824548` differs from
that tag only in `AGENTS.md`; it is code-equivalent, not the same Git commit.

## Rules

- **Clean room.** Never read, grep or summarize 10x Genomics code (Cell
  Ranger, Space Ranger, installed tarballs). Barcode whitelists used as test
  inputs are data. Stop and ask if a question needs vendor source.
- **Patent material.** Never read `/mnt/pikachu/libbfastq`,
  `/mnt/pikachu/fgqzip` or `/mnt/pikachu/zshard*`.
- **Git.** Commit on `dev-release-v1.10.0` only, plain messages, no AI or
  Claude attribution. Do not push anything, tags included. Do not merge into
  `master`. Do not touch other worktrees or branches
  (`feature/dogmaplex-gse309834`, `feat/atac-sidecar-without-bam`,
  `dev-release-v1.9.6`, `master`). The STAR paper cites 1.9.5.
- **Host sharing.** Timed runs and heavy regression batches hold
  `flock /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock <cmd>`. Before a
  timed run check `uptime` and `pgrep -af STAR`; time with
  `/mnt/pikachu/multiomics-suite-nm-refresh/scripts/run_timed.py` and report
  its verdict. Untimed builds: `nice -n 10`, at most `-j16`. Do not use
  `pkill -f` with a pattern that appears in your own command line.
- **G-S3 host-activity policy (owner clarification, 2026-09-29).** Preserve
  host-monitor flags, but do not reject a completed, output-matching workload
  when three-repeat median wall time and peak RSS meet the +3% / +2% limits.
  Keep every observation in the medians. If tolerance is exceeded, flag host
  activity as a possible contributor; do not automatically excuse the failure.
  This does not waive STAR errors, missing outputs or incomplete runs.
- **Design problems.** If a design assumption is wrong, stop and report.
- **Gate failures are not passes.** Matching failures on both binaries can
  establish that a problem predates this refactor, but do not satisfy coverage
  for an RC. Inspect per-case summaries and downstream failure markers, not
  just the batch driver's exit code or `finished_utc`.
- **Identical repeats require explicit owner approval.** Do not launch the
  three G-S3 repetitions just because they appear in this runbook. Preserve
  each existing run directory; the original scratch drivers remove outputs
  unconditionally and must not be used to restart an existing arm.

## Steps

**2026-09-29 candidate creation:** following the owner's request to start
1.10 version creation, create the local annotated `v1.10.0-rc1` tag from the
accepted branch. No remote push or stable merge is included. Record the source
archive and exact identities in `$D/releases/v1.10.0-rc1/LOCAL_RC.json`; see
the [handoff](../handoffs/HANDOFF_STAR_1_10_0_HOST_API_20260928.md) and
[distribution status](../Star-binary-distribution.md). Multiomics strict-build
blockers, its integration gates, canonical downstream integration and multiome
recipe migration remain pre-stable dependencies. Do not move rc1 to include
later fixes; create a new immutable candidate.

**2026-09-29 follow-up:** repeat approval is now recorded. The fresh candidate
passed the host API, expanded Tier A, SLAM SE/PE and seeded downstream checks.
The reused UCSF fixture configuration was repaired and matrix parity verified.
The expanded capture audit uncovered multi-thread TranscriptVB numerical
variation, also reproducible on the pristine baseline. The owner subsequently
directed single-thread, order-sensitive validation. All four 100K SE/PE cases
and both ordered sidecar payloads now match exactly; the kept-output audit is
closed. G-S1 is accepted and G-S2 passed. All eighteen G-S3 attempts match
outputs; median changes are below +3% wall / +2% RSS. **G-S3 is accepted**
under the subsequent owner clarification above, with all six host flags kept.
Wrapper errors were separately verified as harness-only, not failed STAR
workloads; their original exits and dispositions are preserved. Read-only
reassessment: `acceptance_20260928/gs3_report_v2.json` and
`G_S3_ACCEPTED.json` under the gate artifact root. No benchmarks were rerun.
Use immutable driver snapshots for any further authorized timing runs. See
[validation results](../VALIDATION_STAR_1_10_0_20260929.md).
Earlier execution counts and holds below are historical, not current approval
status. No tag, push or master merge occurred during validation; subsequent
local candidate creation is recorded above.

**2026-09-28 correction:** the modern Flex workload is
`tests/run_flex_half_probe_100k_smoke.sh`, using the established H0/H1X2 half-khash
and matching model/included-gene lists. It passed all 121 output comparisons to
saved v1.9.5.a results. The March H0/H1 hash-versus-alignment script is an obsolete
choice for the modern release gate; keep its unresolved drift as a separate
historical diagnostic. No STAR Flex classification change is included in this
follow-up. The historical standalone replay uses an explicit legacy policy.
See the handoff's **Latest follow-up** for exact paths, hashes and commands.

The PF stress suite and a separate five-epoch CUDA/layer smoke also passed.
The sparse UCSF 100K test now tests STAR/downstream H5ADs, with CUDA denoising
validated separately on raw A375 droplets. Seed/failure-handling fixes are in an
isolated canonical recipes worktree and one saved-MEX downstream execution passed.
Repeat proofs remain held for authorization. The G-S3 driver is now fail-closed,
preserves output and uses the modern Flex workload, but has not been launched.
Earlier batch totals below are historical, not retroactively green gates.

| Step | Work | Command / where | Gate |
|---|---|---|---|
| S1 | Branch from `master`; cherry-pick `7d9d71d` (paired hashtag demux, per-library ambient FDR) and `dd14f5b` (CellTag ambient-FDR test) | `git worktree add -b dev-release-v1.10.0 … master`; `git cherry-pick -x` | `make core-portable process-features-lib`; `tests/multi_feature/test_hash_pair_demux_pf_multi.sh`; `tests/test_ambient_fdr_feature_type.sh` |
| S3 (done before S2) | Remove the integration: orchestration and stub, `libchromap_contract`, 57 parameters, `WITH_CHROMAP` and Chromap/RapidMACS Makefile paths, Chromap tests/scripts/MCP/docs, libscrna ATAC parts, `CbqChromapAdapter`; `HTSLIB=bundled|external` replaces `WITH_CHROMAP` | commit `826e019` | clean build with no Chromap or RapidMACS on any path (`strace` of a fresh-export build) |
| S2 | Host interface: `star::host::runMain`, parameter pass-through, lifecycle and permit-plan callbacks, permit facade, `EXTERNAL` rename, `libstar_suite` target, HTSlib knob; `tests/host_api/` | commit `f8d2dae`; `make host-api-tests`; `docs/HOST_API.md` | G-S2: `tests/host_api/run_host_api_tests.sh` |
| Hand-over | `docs/HANDOVER_MULTIOMICS_1.10.md`: every removed file, block, symbol and parameter with its last STAR commit | commit `face24f` | review |
| Version | `core/legacy/source/VERSION` 1.10.0; `debian/changelog` 1.10.0-1 | commit `0b92ce2` | — |
| S4 G-S1 | `make core-portable` of v1.9.5.a vs `make core` of the rc in one container; 24 non-multiome manifest rows, Tier A (13), CI partial builds, release smokes; byte-identical except logs, BAM `@PG`/`@CO`, first line of `genomeParameters.txt` (and the Step 0 allowed items) | tools in `/mnt/pikachu/star_suite_v1100_gates_20260928/tools/` (section below) | new STAR-own tests pass; `test_parameters_default_generation.py` passes |
| S4 G-S3 | `tests/run_scrna_gex_100k_regression.sh`, one 10M-read scRNA run, the Flex 100k run; median of three; wall ≤ +3 %, peak RSS ≤ +2 % vs 1.9.5.a | `run_gs3.sh`; saved-evidence `reassess_gs3.py` | Completion, output parity and numerical limits; host verdicts diagnostic |
| Release | `docs/RELEASE_NOTES_v1.10.0.md` with gate results; `docs/Star-binary-distribution.md` entry; local annotated tag | `git tag -a v1.10.0-rc1 -m "STAR Suite v1.10.0-rc1"` | `python3 scripts/release/validate_official_snapshots.py` |

## Gate tooling

All tools are under `/mnt/pikachu/star_suite_v1100_gates_20260928/tools/`.
The commands below record the original execution layout: `$G` was
`/tmp/claude-1000/-mnt-pikachu-chromap-suite-paper/53e97281-48e6-45f8-90f8-a4913b843536/scratchpad/gate/`.
Both batches are complete; these are not instructions to restart them.
Verified durable copies are now under
`/mnt/pikachu/star_suite_v1100_gates_20260928/gate/`. Logs, fixtures and
selected external diagnostic directories are alongside it; see
`preservation.tsv` and `PRESERVATION_COMPLETE_UTC`. Original scratch paths
and symlinks are retained in the archive; no source outputs were removed.

1. Gate container (same toolchain for both builds; Ubuntu 22.04, g++-12
   12.3.0 as on the host, so binaries also run on the host):
   `docker build -t star-suite-gate-v1100:20260928 -f Dockerfile.gate .`
2. Source trees: `git clone --shared --no-checkout /mnt/pikachu/STAR-suite $G/src195 && git -C $G/src195 checkout --detach 4824548`; same for `$G/src110` at the rc code commit.
3. Build in the container (both trees, `WITH_CHROMAP=0` for 1.9.5.a):
   `docker run --rm --user 1000:1000 -e HOME=/tmp -v $G:$G -w $G/src195 star-suite-gate-v1100:20260928 nice -n 10 bash $G/build_in_container.sh $G/src195 core-portable <sha>`
   and `... src110 core <sha>`.
4. G-S1 runs: `nohup nice -n 10 bash $G/run_gs1.sh 195 &` then `... 110`.
   Each row and Tier A test runs under the lock; STAR is replaced by
   `star_capture_wrapper.sh` (real binary `STAR.real`), which copies each
   invocation's outputs to `$G/run<tree>/capture/<n>`. Row status:
   `$G/run<tree>/manifest_status.tsv`; Tier A: `$G/run<tree>/tierA/status.tsv`.
   Rows print `SKIP` and exit 0 when a prerequisite is missing: check stdout.
5. Compare using `env GATE_ROOT="$G" COMPARE_OUT=<fresh-report-directory>
   bash /mnt/pikachu/star_suite_v1100_gates_20260928/tools/compare_gs1.sh`
   under the shared host lock. Captures must be paired by unique normalized
   `argv.txt`, not by invocation number: held cases change subsequent
   numbering. This driver requires completion markers and propagates missing
   pairs, differing exit status and comparator failures. It does not implement
   every Step 0 semantic normalization; its differences require review.
6. CI partial builds: `run_partial_builds.sh` (both trees, under the lock).
7. G-S3: **held**, not launched. The driver has now been repaired to require
   explicit repeat approval, G-S1 acceptance and the accepted clean candidate
   binary. It refuses existing output directories and propagates failures.
   Review before use; completion of G-S1 alone is insufficient authorization.
   Default records: `$D/gs3_acceptance/<workload>/<tree>/rep<n>/`
   (`time.txt`, `timed/HOST_LOAD.json`). Modern Flex uses the half-probe wrapper.

## Gate audit and current hold

The read-only execution audit is available from this checkout. For a future
changed batch use its root and a new report path:

```bash
python3 tests/host_api/audit_gate_batch.py "$G/run195" \
  --manifest tests/production_module_regression_manifest.tsv \
  --report <fresh-report-path.json>
```

Use `run110` and a different report name for the candidate. Nonzero is expected
when any case fails, is skipped, lacks completion evidence, or has a
`CELLBENDER_FAILED.txt` marker. A successful execution audit still does not
establish output parity or performance. Unit coverage:

```bash
python3 -m unittest discover -s tests/host_api -p 'test_audit_gate_batch.py' -v
```

Completed on 2026-09-28 (12 auditor unit tests passed):

| Arm | Production PASS | FAIL | SKIP | Tier A |
|---|---:|---:|---:|---:|
| 1.9.5.a baseline | 18 | 6 | 0 | 13/13 |
| 1.10.0 candidate | 18 | 4 | 2 | 13/13 |

Reports in `/mnt/pikachu/star_suite_v1100_gates_20260928/`:
`gs1_baseline_execution_audit_v2.json`,
`gs1_candidate_execution_audit.json`, `compare/`,
`selected_output_differences.json`, and `BINARY_PROVENANCE.json`.
The four candidate failures also occur on the baseline. Two candidate SLAM
rows are held for explicit approval because they repeat identical FASTQ
executions internally; their gate-local manifest requires
`IDENTICAL_REPEATS_APPROVED` in that artifact root. Do not create this marker
without recording actual owner approval.

Selected comparisons establish equal PBMC matrices and keyed A375 feature
counts. UCSF `counts.h5ad` dataset values differ only in source-path provenance,
but downstream doublet identities and scores differ. The executed canonical
`morphic-recipes` R caller lacks the seed present in STAR's unused copy.
Resolve this external recipe reproducibility gap with its owner; do not
waive annotation parity. Both UCSF runs also fall back after CellBender
fails during prior estimation, so GPU denoising is not validated.
See the [handoff](../handoffs/HANDOFF_STAR_1_10_0_HOST_API_20260928.md)
for precise failures, comparison limits and remaining actions.

**Historical G-S3 diagnosis (superseded by the correction above).** The baseline batch exposed
an invalid legacy Flex invocation: `run_flex_hash_screen_internal_100k.sh`
sets `--flexLegacy yes` and legacy expected-cell tuning without selecting
`--soloFlexCellCaller legacy`. That flag is now fixed. Follow-up runs complete
STAR but still fail hash-on vs legacy matrix parity in both versions, with
identical discrepancies; cross-version hash-on counts match exactly. Resolve
this test contract and the remaining G-S1 failures, obtain repeat approval,
and make the timing driver preserve existing output directories and propagate
failures before proceeding. The follow-up also repairs the OCM/Velocyto
allocation guard, updates obsolete permit assertions, and isolates Flex replay
drift to a negative-cache policy change. See the handoff for individual
validation results; the historical batch totals above are not overwritten.

The official snapshot validator checks digest/count integrity, not whether
the six pinned multiome recipes have migrated to the Multiomics executable.
Keep that migration as a separate release dependency.
