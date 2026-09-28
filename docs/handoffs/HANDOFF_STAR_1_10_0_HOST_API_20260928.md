# Handoff: STAR Suite 1.10.0 host interface and Chromap removal

Updated: 2026-09-28 18:33 UTC. Runbook:
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
- Continuation changes: execution-coverage auditor and its unit tests,
  updated runbook warnings, and explicit unreleased status in the release
  notes. STAR's C++ implementation has not changed during this continuation.

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

### G-S1 continuation

- Baseline completed at 18:17 UTC (`finished_utc`; `done 195` in
  `$G/../logs/gs1_195.log`). All 24 production rows ran: 18 exited zero,
  six failed. All 13 Tier A tests passed. Do not restart this arm.
- Candidate started at 18:19:05 UTC with
  `nice -n 10 bash "$G/run_gs1.sh" 110`, under the driver's per-row shared
  host lock. At 18:31 UTC rows 1-16 were complete; the OCM failure reproduced.
  Inspect `$G/run110/manifest_status.tsv`, `$G/run110/tierA/status.tsv`, and
  `$G/run110/finished_utc`. Do not start a second candidate batch.
- First PBMC capture: 22 non-log files identical, no differing/missing files.
  The PBMC 100K `report.json` inputs and all three profile QC/matrix-checksum
  objects match (vanilla, modern, modern_bam). This is partial comparison
  evidence, not a full G-S1 pass.
- Partial builds finished (both trees).
- Gate binaries (built in container `star-suite-gate-v1100:20260928`):
  `$G/src195/core/legacy/source/STAR.real` (1.9.5.a, `4824548`) and
  `$G/src110/core/legacy/source/STAR[.real]` (1.10.0, `0b92ce2`).
- 10M-read fixture for G-S3: `$G/../gs3/pbmc10k_10M/L00{1,2}_R{1,2}.fastq.gz`
  (first 5M read pairs of each PBMC 10k v3 lane).

## Next

1. Let the existing candidate batch finish; do not restart either arm. Audit
   both completed batches with `tests/host_api/audit_gate_batch.py` (runbook).
   The new auditor's 12 unit tests pass. Exit zero from the original batch
   driver means only that it reached its end, not that all cases passed.
2. Compare captures and kept outputs with `compare_outputs.py` (runbook,
   step 5). Allowed differences: logs, BAM `@PG`/`@CO`, first line of
   `genomeParameters.txt`, the Step 0 content-varying items (example read in
   `feature_sequences.txt`, EmptyDrops tie order), and fields that record the
   binary or version (report.json `binary_sha256`/`version`/`source_revision`,
   sidecar metadata). Anything else is a failure to investigate.
3. G-S3 is on hold. Resolve the invalid Flex recipe below and obtain explicit
   owner approval for identical timing repeats. The question has been asked;
   no approval has been received. Do not run the existing driver unchanged:
   it overwrites runs and continues after failures. The 100K GEX workload has
   already run during G-S1, so timing it again also requires repeat approval.
4. Complete the output comparisons and disposition the failures before an RC
   tag. Keep release notes and distribution docs explicitly pending until
   gates are satisfied. Commit locally only; no tag, merge or push yet.

## Baseline failures requiring disposition

All findings below come from the clean, container-built baseline, not from
stale objects. Matching failures on 1.10 would establish that they predate
this refactor, not prove the affected workflow works.

| Row | Observed baseline failure | Required follow-up |
|---|---|---|
| `cbq-ocm-composite-local` | STAR exits 111: Velocyto's gene-like source did not populate per-read CB/UMI storage | Already reproduced on 1.10; investigate OCM/Velocyto coverage separately |
| `pf-dynamic-permit-100k` | Validator requires acquires == workUnits (7,264 vs 369,546); assignment now reports batched work, up to 64 records per permit | Correct the telemetry invariant with unit coverage; do not change the batching to satisfy an obsolete assertion |
| `flex-hash-screen-replay` | Flat and tiered agree, but both differ from pinned truth on 2,024/800,000 reads, all Pass -> Deny (negative codes 6: 1,873; 4: 136; 8: 15) | Establish the intended negative-cache contract before updating a golden fixture; do not overwrite the pinned truth |
| `flex-hash-screen-100k` | Exits 102 before mapping: legacy expected-cell tuning with default tag-aware caller | Explicit legacy caller is missing from the legacy test recipe; this also invalidates the planned G-S3 Flex timing command |
| `slam-cbq-divergence` | FASTQ vs CBQ SAM body order differs | Exact multiset is identical (113,516 records), and staged sorted-SAM, diagnostics and pre-NTR checks pass; decide whether record-order parity is required |
| `slam-cbq-divergence-pe` | Same exact-order failure | Staged sorted-SAM, junctions, diagnostics, transitions and pre-NTR checks pass |

Additional coverage notes:

- UCSF smoke exits zero by design when CellBender fails and writes a fallback
  H5AD (`tests/run_ucsf_corrected_production_100k_smoke.sh:112`). It did request
  `--cuda`; the 100K fixture has only 5,584 observed barcodes and CellBender did
  not produce its H5. This tests the fallback, not denoising. The failure
  marker is under `run195/ucsf/samples/EBs2_2/`
  `downstream_genefull_velocyto_cellbender/cellbender/`.
- The CBQ aggregate skips two network subtests with `RUN_NETWORK=0`; their
  separate manifest rows did run successfully. The auditor surfaces nested
  skips for review rather than silently ignoring them.
- `python3 scripts/release/validate_official_snapshots.py` passed: 11 recipes,
  10 public evidence records. It does not validate recipe executable routing.
- Baseline `4824548` differs from the published `v1.9.5.a` tag `f0d9f27` only
  in `AGENTS.md`; biological code is identical.

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
