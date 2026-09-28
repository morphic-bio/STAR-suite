# Handoff: STAR Suite 1.10.0 host interface and Chromap removal

Updated: 2026-09-28 17:50 UTC. Runbook:
`docs/runbooks/RUNBOOK_STAR_1_10_0_HOST_API_20260928.md`.

## State

- Worktree `/mnt/pikachu/STAR-suite-v1100-20260928`, branch
  `dev-release-v1.10.0`. Nothing pushed; no tag yet; `master` untouched.
- Commits on the branch (oldest first):
  - `10bb862` cherry-pick of `7d9d71d` (paired hashtag demux, per-library ambient FDR)
  - `90bb0e0` cherry-pick of `dd14f5b` (CellTag ambient-FDR test)
  - `826e019` removal of the Chromap integration; `HTSLIB=` knob (S3)
  - `f8d2dae` host interface, `tests/host_api/`, `docs/HOST_API.md` (S2)
  - `face24f` `docs/HANDOVER_MULTIOMICS_1.10.md`
  - `0b92ce2` version 1.10.0, `debian/changelog` 1.10.0-1
  - this runbook and handoff
- Uncommitted: none, except the draft `docs/RELEASE_NOTES_v1.10.0.md`
  (validation section pending), committed with this handoff as a draft.

### Done

| Item | Result | Evidence |
|---|---|---|
| S1 | both cherry-picks match their originals (ARTIFACTS.md conflict resolved; the ATAC-sampling bullet of `52c1ec8` left out); tests pass | `tests/multi_feature/test_hash_pair_demux_pf_multi.sh`, `tests/test_ambient_fdr_feature_type.sh` |
| S3 clean build | fresh `git archive` of `0b92ce2`, `make -j16 core star-host-lib` under `strace -f -e trace=%file,execve`: 0 paths containing chromap or rapidmacs (20.2M trace lines, scratch prefix masked); 0 such strings in `STAR` and `libstar_suite.a` | `/mnt/pikachu/star_suite_v1100_gates_20260928/S3_clean_build_summary.txt` |
| G-S2 | pass, bundled and `HTSLIB=external` builds; also as partial-build case `host-api` | `make host-api-tests`; runs in scratch `host_api_run1`, `ext_hostapi`, `partial/110/host-api/build.log` |
| Docker builder sequence | pass (`make core; flex; slam; feature-barcodes-tools; default; all`; unit tests) on a fresh export | scratch `logs/docker_seq.log` |
| `test_parameters_default_generation.py`, `test_htslib_build_discovery.py` (9 tests), MCP tests (589) | pass | — |
| CI partial builds, 1.10 | 11/11 PASS (core, core-static, core-htslib, host-api, release-companion-tools, star-feature-call, feature-barcodes-tools, process-features-lib, yremove-tools, core-external-htslib, core-portable) | `$G/partial/110/summary.tsv` |
| CI partial builds, v1.9.5.a baseline | 8/8 PASS (its default list) | `$G/partial/195/summary.tsv` |

`$G` = `/tmp/claude-1000/-mnt-pikachu-chromap-suite-paper/53e97281-48e6-45f8-90f8-a4913b843536/scratchpad/gate`.

### In progress (background jobs)

- `bash $G/run_gs1.sh 195` (PID 1061295), started 17:28 UTC: G-S1 baseline
  rows and Tier A with the v1.9.5.a gate binary. Log
  `$G/../logs/gs1_195.log`; status `$G/run195/manifest_status.tsv`,
  `$G/run195/tierA/status.tsv`; captures `$G/run195/capture/`. Rows 1-8
  exit 0. Row `cbq-ocm-composite-local` exits 1 on the v1.9.5.a baseline
  itself ("Velocyto's gene-like source did not populate per-read CB/UMI
  storage", with the current morphic-recipes OCM script): pre-existing; 1.10
  must fail the same way.
- Partial builds finished (both trees).
- Gate binaries (built in container `star-suite-gate-v1100:20260928`):
  `$G/src195/core/legacy/source/STAR.real` (1.9.5.a, `4824548`) and
  `$G/src110/core/legacy/source/STAR[.real]` (1.10.0, `0b92ce2`).
- 10M-read fixture for G-S3: `$G/../gs3/pbmc10k_10M/L00{1,2}_R{1,2}.fastq.gz`
  (first 5M read pairs of each PBMC 10k v3 lane).

## Next

1. When `run_gs1.sh 195` prints `done 195` in its log: `cd $G && nohup nice -n 10 bash run_gs1.sh 110 > $G/../logs/gs1_110.log 2>&1 &`.
2. Compare captures and kept outputs with `compare_outputs.py` (runbook,
   step 5). Allowed differences: logs, BAM `@PG`/`@CO`, first line of
   `genomeParameters.txt`, the Step 0 content-varying items (example read in
   `feature_sequences.txt`, EmptyDrops tie order), and fields that record the
   binary or version (report.json `binary_sha256`/`version`/`source_revision`,
   sidecar metadata). Anything else is a failure to investigate.
3. G-S3: `cd $G && nohup bash run_gs3.sh > $G/../logs/gs3.log 2>&1 &` after
   both G-S1 runs; summarize medians from `time.txt` (Elapsed, Maximum
   resident set size) and `timed/HOST_LOAD.json` verdicts.
4. Fill the validation section of `docs/RELEASE_NOTES_v1.10.0.md`, add the
   1.10.0 entry to `docs/Star-binary-distribution.md`, commit, then
   `git tag -a v1.10.0-rc1 -m "STAR Suite v1.10.0-rc1"`. No push.

## Deviations from the design

- **Order:** S3 (removal) was committed before S2 so every commit builds and
  behaves coherently; the end state is the design's.
- **57 parameters, not 58:** `chromapAtacBarcodeSampleLimit` (`52c1ec8`) was
  never on `master`; the hand-over lists it from `6e83853`.
- **Help and log text:** the `atacController=` log field is dropped (the host
  logs its own controller); the BGZF-hierarchy and floor error messages no
  longer mention ATAC.
- **Extras:** `SaturationPermitController` gained label-aware name overloads;
  `runMain` refuses a second call per process; `make host-api-tests` and a
  `host-api` partial-build case; the MuData half of the ambient-FDR test left
  with the MuData builder (`tests/test_ambient_fdr_anndata.py` keeps the AnnData checks).
- **G-S1 execution:** the gate container builds both trees; tests run on the
  host with those binaries (same glibc). The UCSF row hard-wires CPU
  CellBender, which AGENTS.md forbids, and defaults to the tracked
  `STAR.release`; the gate runs it with the gate STAR and `--cellbender-gpu`
  (identical gate-local edit on both sides).
- **G-S3 "Flex 100k replay"** is taken to be `tests/run_flex_hash_screen_internal_100k.sh`
  (STAR); `run_flex_hash_screen_replay_regression.sh` uses a standalone tool
  with no STAR code and unchanged sources.

## Open problems

- `share/star-suite/SNAPSHOTS.json` still pins STAR-suite-recipes `1b12325`
  with six multiome workflows; re-pinning needs a recipes revision that points
  at the Multiomics binary (after Multiomics M3).
- Parallel `make -j all` races at top level on v1.9.5.a as well (pre-existing);
  the Docker builder's serial sequence passes.
