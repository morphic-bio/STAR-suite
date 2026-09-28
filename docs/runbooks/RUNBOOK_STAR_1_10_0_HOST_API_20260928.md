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
- **Design problems.** If a design assumption is wrong, stop and report.

## Steps

| Step | Work | Command / where | Gate |
|---|---|---|---|
| S1 | Branch from `master`; cherry-pick `7d9d71d` (paired hashtag demux, per-library ambient FDR) and `dd14f5b` (CellTag ambient-FDR test) | `git worktree add -b dev-release-v1.10.0 … master`; `git cherry-pick -x` | `make core-portable process-features-lib`; `tests/multi_feature/test_hash_pair_demux_pf_multi.sh`; `tests/test_ambient_fdr_feature_type.sh` |
| S3 (done before S2) | Remove the integration: orchestration and stub, `libchromap_contract`, 57 parameters, `WITH_CHROMAP` and Chromap/RapidMACS Makefile paths, Chromap tests/scripts/MCP/docs, libscrna ATAC parts, `CbqChromapAdapter`; `HTSLIB=bundled|external` replaces `WITH_CHROMAP` | commit `826e019` | clean build with no Chromap or RapidMACS on any path (`strace` of a fresh-export build) |
| S2 | Host interface: `star::host::runMain`, parameter pass-through, lifecycle and permit-plan callbacks, permit facade, `EXTERNAL` rename, `libstar_suite` target, HTSlib knob; `tests/host_api/` | commit `f8d2dae`; `make host-api-tests`; `docs/HOST_API.md` | G-S2: `tests/host_api/run_host_api_tests.sh` |
| Hand-over | `docs/HANDOVER_MULTIOMICS_1.10.md`: every removed file, block, symbol and parameter with its last STAR commit | commit `face24f` | review |
| Version | `core/legacy/source/VERSION` 1.10.0; `debian/changelog` 1.10.0-1 | commit `0b92ce2` | — |
| S4 G-S1 | `make core-portable` of v1.9.5.a vs `make core` of the rc in one container; 24 non-multiome manifest rows, Tier A (13), CI partial builds, release smokes; byte-identical except logs, BAM `@PG`/`@CO`, first line of `genomeParameters.txt` (and the Step 0 allowed items) | tools in `/mnt/pikachu/star_suite_v1100_gates_20260928/tools/` (section below) | new STAR-own tests pass; `test_parameters_default_generation.py` passes |
| S4 G-S3 | `tests/run_scrna_gex_100k_regression.sh`, one 10M-read scRNA run, the Flex 100k run; median of three; wall ≤ +3 %, peak RSS ≤ +2 % vs 1.9.5.a | `run_gs3.sh` | `run_timed.py` verdicts |
| Release | `docs/RELEASE_NOTES_v1.10.0.md` with gate results; `docs/Star-binary-distribution.md` entry; local annotated tag | `git tag -a v1.10.0-rc1 -m "STAR Suite v1.10.0-rc1"` | `python3 scripts/release/validate_official_snapshots.py` |

## Gate tooling

All under `/mnt/pikachu/star_suite_v1100_gates_20260928/tools/`; working
copies and outputs are in the session scratchpad
`/tmp/claude-1000/-mnt-pikachu-chromap-suite-paper/53e97281-48e6-45f8-90f8-a4913b843536/scratchpad/gate/` (`$G`).

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
5. Compare: `python3 compare_outputs.py $G/run195/capture $G/run110/capture --roots-a '<T>=$G/src195' '<R>=$G/run195' --roots-b '<T>=$G/src110' '<R>=$G/run110' --report <json>`,
   and likewise the kept per-test output trees.
6. CI partial builds: `run_partial_builds.sh` (both trees, under the lock).
7. G-S3: `nohup bash $G/run_gs3.sh &` after G-S1 finishes (run_timed refuses
   while another STAR runs). Records in `$G/gs3/<workload>/<tree>/rep<n>/`
   (`time.txt`, `timed/HOST_LOAD.json`).
