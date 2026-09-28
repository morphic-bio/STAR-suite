# Handoff: STAR Suite 1.10.0 host interface and Chromap removal

Updated: 2026-09-28, after both G-S1 batches completed. Runbook:
`docs/runbooks/RUNBOOK_STAR_1_10_0_HOST_API_20260928.md`.

## State

### Latest follow-up: current Flex fixture corrected

This section supersedes the earlier interpretation of the legacy Flex failure
as a blocker to the current 100K route. The old batch records below are retained,
not rewritten into passes.

The failing `run_flex_hash_screen_internal_100k.sh` uses a March-2026 H0/H1
cache and the 2020 alignment reference. It is **not** the established 1.9.5.a
JAX half-probe regression. Selecting it as G-S3's Flex workload was a test-selection
error. Correct current inputs, under
`/mnt/pikachu/storage_offload_20260912/inputfmt_20260910/stage/`:

| Input | SHA256 |
|---|---|
| `model_h01x2_cache.half.khash` | `801fc6143a383b6a94c55307f816bf824c7a9c685d3049405b2ffcc4d72db296` |
| `model_gene_ids.txt` | `9f19693ee17820f068ce0dec1d27a905720cfe5912e2215c7d023d55feeb13aa` |
| `included_gene_ids.txt` | `5366fd7971227e5b9fab720affc0d93eab287c22254eaa01af35efd58f1ba431` |

Use `tests/run_flex_half_probe_100k_smoke.sh`, which invokes the established
`test_flex_v194_default_route.sh` in new single-run mode. It executes the current
binary once and compares to preserved successful 1.9.5.a output, without rerunning
the reference or the default/explicit equivalence arms. It uses all eight lanes,
100K pairs per lane, the matching model/filter gene axes, the tag-aware caller,
and the fused no-genome route. It is not a hash-versus-genomic-alignment test.
The reference is copied unchanged from
`/tmp/star-flex-removal-20260927/jax-parity/B` to the durable
`/mnt/pikachu/star_suite_v1100_gates_20260928/flex_modern_reference_v195a/`.
Its original logs retain their original paths.

Completed modern validation: **121/121 non-log, non-diagnostic outputs
byte-identical**, including raw/filtered matrices and caller outputs; seven
checks passed. Evidence: `$D/fixes_20260928/flex_modern_100k/`. This first run
used clean binary `ee0bfc87fd0d329247f42e23fbf04c8d6a9e83b1e8974a11d6420c080293b459`.
An inactive legacy helper was temporarily changed in that build; modern routing
was unchanged. All experimental STAR Flex source changes were subsequently
removed. Do not use the experiment binary for release.

After removing all production Flex edits, another clean build and the single-run
wrapper passed **121/121** again. Authoritative final evidence:
`$D/fixes_20260928/flex_modern_final/`; tested binary SHA256
`2746e2037d8632d8c516a195283cda830043ceb3e869ec86e834879d19704512`.
Command: `OUT_ROOT=$D/fixes_20260928/flex_modern_final bash tests/run_flex_half_probe_100k_smoke.sh`
under the shared host lock and nice 10. This was changed-source validation after
removing the experiment, not a timing repetition.

The historical replayer now has explicit `--legacy-negative-policy` for its
immutable March oracle. Its default certified-negative policy is preserved;
synthetic tests check both. This changes a standalone test tool, not STAR's
runtime classifier. Historical replay: 800,000/800,000 decisions agree, 20
synthetic checks pass. Evidence: `$D/fixes_20260928/flex_replay_versioned/`.
Do not use that historical success as a modern count-validation claim.

A temporary production-routing experiment did **not** solve historical count
drift (6,761 differing pooled coordinates, 0.9903%; previous 1.9299%). It was
removed, with the patch/binary and outputs retained under
`$D/fixes_20260928/{legacy_routing_experiment.patch,STAR.legacy-routing-experiment,flex_e2e/}`.
The old diagnostic remains available and failing; its exact cause is unresolved.
Never label it a benign difference merely because the baseline also failed.

### Other fixes and validation

- PF permit stress completed successfully with the earlier clean OCM-fix binary
  (`c74e9019...`): off/on MEX parity, 3->2->4 and 1->2->1 resizing,
  shadow/active controller, timeout and recovery. Evidence: `fixes_20260928/pf_stress/`.
- Canonical downstream edits are isolated in
  `/mnt/pikachu/morphic-recipes-v1100-20260928`, branch
  `fix/v1100-downstream-validation`, local commit `bbf8b54`; the main recipes worktree's cardiac edits
  are untouched. Set `MORPHIC_RECIPES_ROOT` to that worktree when testing STAR's
  compatibility launcher. The STAR R mirror matches the canonical helper.
- The R container now receives `SCDBLFINDER_SEED`, default 1; scDblFinder uses
  R's seed plus `SerialParam(RNGseed=...)`. Errors fail instead of marking all
  cells singlets. Saved-UCSF-MEX downstream validation completed with 3,943
  STAR cells and 281 doublets. This is one seeded run, not repeated-run proof.
  Evidence: `fixes_20260928/seeded_downstream/` and `seeded_downstream.log`.
- CellBender requires CUDA and successful nonempty output. Partial failures and
  stale/reused failed output cannot become successful fallback H5ADs. Ten
  orchestration tests pass in recipes. The 100K UCSF smoke now validates
  alignment/H5ADs without pretending its sparse input validates denoising.
- Separate CUDA smoke passed on 20,000 raw A375 droplets x 38,606 GEX features
  from `/storage/A375/paper_bench_20260326_134444/outs/raw_feature_bc_matrix`.
  Raw UMIs: 21,546,719; denoised: 20,993,095; denoised NNZ: 6,267,577.
  Raw X and barcode/gene order remain unchanged. `nvidia-smi` captured 1,104 MiB
  for CellBender. Five epochs test execution/layer integration, **not production
  convergence**. Evidence: `fixes_20260928/cellbender_cuda/{manifest.json,PASS.json,GPU_ACTIVE.txt}`.
- Downstream image ID: `sha256:c052c1568727f24a6e2c8ec793357ccfdf1ffdf5622542a7b2da40d8df47825b`;
  CellBender: `sha256:f31f1e993f3d87659c3dd541a22f505fec7ee4366f6d5da1324061c2c09dcb2a`.
- Five keyed-MEX comparator tests pass; automatically compare the emitted sample
  set without assuming the shallow fixture calls cells in all four tags. A
  missing sample on one side, missing matrix, or raw count drift still fails.
- `$D/tools/run_gs3.sh` now resolves the durable paths, uses the modern half-probe
  test, refuses output overwrites, propagates failures, records binary hashes,
  and checks repeat-approval/G-S1-acceptance records before any execution. Its
  no-approval dry check exits 2 without running a workload. Actual timings remain held.

### Remaining release work

Explicit repeat approval is still needed for G-S3, SE/PE SLAM determinism and the
second seeded doublet run. No approval markers were created. Finish output-pairing
normalization and the complete gate audit; the initial failed batches remain
failed historical records. The new 100K smoke scope and separate CUDA test are
not retroactive passes of the old combined wrapper. Multiomics integration and
its recipe executable/snapshot migration remain stable-release dependencies
after the local RC is available, not a reason to change the current Flex cache.
No push, tag or master merge is authorized by this follow-up.

Completed execution commands (records, not rerun authorization):

```bash
# Each data execution used flock $D/../e2e_bench_20260926/pikachu_timed.lock
# (actual lock: /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock), nice -n 10.
PF_DYNAMIC_100K_OUT_BASE=$D/fixes_20260928/pf_stress bash tests/run_pf_dynamic_permit_100k_smoke.sh
OUT_ROOT=$D/fixes_20260928/flex_replay_versioned bash tests/run_flex_hash_screen_replay_regression.sh
FLEX_SINGLE_RUN=1 FLEX_REFERENCE_OUTPUT=/tmp/star-flex-removal-20260927/jax-parity/B \
  TEST_WORKDIR=$D/fixes_20260928/flex_modern_100k bash tests/test_flex_v194_default_route.sh
MORPHIC_RECIPES_ROOT=/mnt/pikachu/morphic-recipes-v1100-20260928 \
  CELLBENDER_SMOKE_OUT=$D/fixes_20260928/cellbender_cuda bash tests/run_cellbender_cuda_smoke.sh
# In the isolated recipes worktree:
SCRNA_DOWNSTREAM_IMAGE=sha256:c052c1568727f24a6e2c8ec793357ccfdf1ffdf5622542a7b2da40d8df47825b \
  SCDBLFINDER_SEED=1 bash scripts/run_scrna_downstream_gene_full_velocyto.sh \
  --run-dir $D/gate/run110/ucsf/samples/EBs2_2/run \
  --output-dir $D/fixes_20260928/seeded_downstream --python-backend host
```

### Earlier branch history

- Worktree `/mnt/pikachu/STAR-suite-v1100-20260928`, branch
  `dev-release-v1.10.0`. Nothing pushed; no tag yet; `master` untouched.
- Commits on the branch (oldest first):
  - `10bb862` cherry-pick of `7d9d71d` (paired hashtag demux, per-library ambient FDR)
  - `90bb0e0` cherry-pick of `dd14f5b` (CellTag ambient-FDR test)
  - `826e019` removal of the Chromap integration; `HTSLIB=` knob (S3)
  - `f8d2dae` host interface, `tests/host_api/`, `docs/HOST_API.md` (S2)
  - `face24f` `docs/HANDOVER_MULTIOMICS_1.10.md`
  - `0b92ce2` version 1.10.0, `debian/changelog` 1.10.0-1
  - `a321d99` initial runbook, handoff and draft release notes
  - `c943825` partial-build checkpoint
  - `c762817` fail-closed execution audit, 12 unit tests and gate warnings
  - `35d3650` completed gate-batch results and comparison limitations
- Continuation changes: execution-coverage auditor and its unit tests,
  updated runbook warnings, and explicit unreleased status in the release
  notes. The initial gate continuation did not change STAR C++; the subsequent
  OCM storage fix described below does and requires new-binary validation.

### Owner-Requested Failure Follow-Up

Evidence root: `$D/followup_20260928/`; reproduction command:
`flock /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock nice -n 10 bash
$D/tools/run_followup_20260928.sh`. This script was run once; do not restart it.

- **Permit assertions updated.** Acquisitions count batches rather than
  records; GEX and feature assignment can overlap, so MAP activity in a PF
  snapshot is valid. The old equality and MAP-zero assertions were retired.
  Nonnegative integer counters, positive feature work/acquisitions, hook
  activation and the feature-wait ceiling remain checked. Do not add an
  exact cross-domain sum assertion: the legacy snapshot counters are sampled
  atomically but are not one coherent cross-domain snapshot.
  Nine new validator tests pass; saved baseline/candidate off/on records all
  pass via `--validate-api-run FILE 0|1`. No unchanged PF alignments were rerun.
  The original batch's unexecuted stress modes remain missing coverage.
- **Flex replay cause isolated, original truth retained.** The same
  checksum-pinned cache and 800,000-read March dump were replayed with fresh
  standalone builds immediately before (`ff8d53d`) and at (`da1b81d`) the
  2026-09-04 change. Before: zero mismatches. After: 2,024, with the diff TSV
  byte-identical to the current gate's. That commit makes unresolved/negative
  cache hits DENY instead of PASS-to-alignment. It is not a different cache or
  duplicate-handling change in this test. KEEP stays 670,450; DENY changes
  5,396 -> 7,420 and PASS changes 124,154 -> 122,130. A versioned current-policy
  oracle is still needed; no golden file or threshold was silently updated.
- **Legacy Flex invocation corrected.** The smoke now explicitly selects
  `--soloFlexCellCaller legacy`. Both old failed invocations are preserved.
  Both binaries complete mapping/counting, but their hash-on vs legacy
  comparator still fails: 13,178 pooled coordinates differ (1.92991%, limit
  0.2%), 250 barcodes differ (0.48096%, limit 0.1%), and maximum coordinate
  delta is 3 (limit 1). Baseline and candidate have the same discrepancies.
  Cross-version hash-on comparison is exact for both pooled counts (682,831
  coordinates; 689,934 UMIs) and the only emitted sample, BC006 (618
  coordinates; 672 UMIs). Cross-version legacy comparison is also exact
  (pooled: 682,742 coordinates / 689,861 UMIs; BC006: 620 / 674).
  This is a pre-existing within-version parity issue,
  not a passing smoke. No tolerances were relaxed.
- **OCM allocation fix.** `resetPackedStorage()` suppressed required
  storage for inline CB correction without BAM tags or the hash bridge.
  A narrow change honors `readInfoYes[featureType]`. A synthetic no-BAM
  control produces exact GeneFull and velocity counts, then reproduces exit
  111 when inline correction is enabled on the original clean binary.
  After a clean rebuild, plain, inline-corrected and native OCM modes pass
  exact GeneFull/spliced/unspliced/ambiguous checks, both normally and with
  `STAR_VELOCYTO_LOW_MEM=1` (six cases). Each barcode has 2 GeneFull molecules,
  1 spliced and 1 unspliced, and 0 ambiguous. The real CBQ arm also completes:
  the three raw GeneFull MEX files are byte-identical to the original 1.10
  candidate's pre-failure output; Velocyto gene/barcode axes match GeneFull.
  Pooled velocity UMIs are 367 spliced, 106 unspliced, 55 ambiguous. Per-sample
  splits (same order) are GCM1 103/25/11, GRHL1 84/25/17, OVOL1 104/30/13,
  WT-PrS-20pct 76/26/14, preserving all 528 molecules. CBQ completed at
  20:38:55 UTC and FASTQ at 20:49:06 UTC. Both have `STAR_COMPLETED.txt`.
  Their `Solo.out`, `outs`, and `samples` trees are byte-identical, including
  routed velocity matrices. The separate ordinary-GEX analytical regression
  also passes all eight UMI/counting configurations on the rebuilt binary.
  This 1,000-read fixture calls zero filtered cells; it validates raw counting
  and materialization, not production-scale cell-calling sensitivity.
  This is a no-BAM allocation correction, not a change to velocity counting
  rules or cell calling. The guard against missing readInfo stays in place.
  The new regression and permit-validator tests are added to Tier A (now
  15 cases). The original 13-case G-S1 capture audit remains frozen historical
  evidence, not the expanded CI list.
- **OCM wrapper follow-up, not an alignment failure.** After both OCM runs and
  the first two byte comparisons passed, the wrapper exited 1 on its final
  expectation that CBQ/Y-removal preparation must be rejected. The canonical
  recipe now supports this and correctly renders `--emitYNoY yes` plus
  `--emitYNoYFormat cbq`, without a BAM or FASTQ sidecar. Update the wrapper
  to validate that preparation contract and add per-sample output parity.
  Read-only `--validate-yremove-plan` passes the saved Y-removal script and
  rejects the saved no-Y-output script as expected. Native CBQ partitioning
  has its own `tests/run_cbq_ynoy_smoke.sh`; it was not run in this follow-up.
  The original wrapper failure/log is retained, and the driver's
  `finished_utc` is intentionally absent. No identical OCM alignments were
  repeated just to obtain a green wrapper exit. This is component validation,
  not a retrospectively passing full smoke or completed release gate.

Follow-up commands and provenance:

```bash
D=/mnt/pikachu/star_suite_v1100_gates_20260928
python3 -m unittest discover -s tests -p 'test_pf_dynamic_permit_validation.py' -v
flock /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock nice -n 10 \
  bash "$D/tools/run_followup_20260928.sh"
flock /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock nice -n 10 \
  bash "$D/tools/run_ocm_storage_fix_validation.sh" --resume-build
python3 tests/compare_flex_hash_screen_mex.py \
  "$D/followup_20260928/flex_195/hash_on" \
  "$D/followup_20260928/flex_110/hash_on" --samples BC006 \
  --max-mismatch-fraction 0 --max-count-delta-fraction 0 \
  --max-barcode-difference-fraction 0 --max-coordinate-delta 0
python3 tests/compare_flex_hash_screen_mex.py \
  "$D/followup_20260928/flex_195/legacy" \
  "$D/followup_20260928/flex_110/legacy" --samples BC006 \
  --max-mismatch-fraction 0 --max-count-delta-fraction 0 \
  --max-barcode-difference-fraction 0 --max-coordinate-delta 0
bash tests/run_cbq_ocm_composite_smoke.sh --validate-yremove-plan \
  "$D/ocm_storage_fix_validation/host_ocm/cbq_yremove_gate/RUN_STAR_COMPOSITE.sh"
# Negative control: expected exit 1, because Y/noY output is not enabled.
bash tests/run_cbq_ocm_composite_smoke.sh --validate-yremove-plan \
  "$D/ocm_storage_fix_validation/host_ocm/cbq/RUN_STAR_COMPOSITE.sh"
diff -qr "$D/ocm_storage_fix_validation/host_ocm/cbq/star_composite/samples" \
  "$D/ocm_storage_fix_validation/host_ocm/fastq/star_composite/samples"
# Executed from $D/ocm_storage_fix_validation/gex_control_work (isolated logs).
flock /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock nice -n 10 \
  python3 /mnt/pikachu/STAR-suite-v1100-20260928/tests/test_scrna_gex_counts.py \
  --star /mnt/pikachu/STAR-suite-v1100-20260928/core/legacy/source/STAR \
  --outdir "$D/ocm_storage_fix_validation/gex_control"
```

These commands record completed executions, not permission to
repeat them. The OCM driver first cleaned and built with a mistyped encoder
target (`cbq_ordered_encoder` instead of `cbq-ordered-encoder`); no STAR test
started in that failed build attempt. Its explicit resume verified unchanged
source and continued the same clean objects with the corrected target. Both
logs are preserved. The resulting host-built STAR SHA256 is
`c74e9019a46216080a61bb836e1d6ff7c24856cede56b7aa983d991edeb03ab0`;
source is `35d3650` plus the saved patch in
`$D/ocm_storage_fix_validation/source.patch`. This binary is distinct from the
initial container-built G-S1 candidate. The complete rebuild commands, fixture
commands, logs and binary checksum are under `ocm_storage_fix_validation/`.
An independent copy of the checked executable is preserved as `STAR.tested`
in that directory, with the same SHA256.

### Done

| Item | Result | Evidence |
|---|---|---|
| S1 | both cherry-picks match their originals (ARTIFACTS.md conflict resolved; the ATAC-sampling bullet of `52c1ec8` left out); tests pass | `tests/multi_feature/test_hash_pair_demux_pf_multi.sh`, `tests/test_ambient_fdr_feature_type.sh` |
| S3 clean build | fresh `git archive` of `0b92ce2`, `make -j16 core star-host-lib` under `strace -f -e trace=%file,execve`: 0 paths containing chromap or rapidmacs (20.2M trace lines, scratch prefix masked); 0 such strings in `STAR` and `libstar_suite.a` | `/mnt/pikachu/star_suite_v1100_gates_20260928/S3_clean_build_summary.txt` |
| G-S2 | pass, bundled and `HTSLIB=external` builds; also as partial-build case `host-api` | `$D/ext_hostapi`, `$D/logs/ext_hostapi.log`, `$G/partial/110/host-api/build.log` |
| Docker builder sequence | pass (`make core; flex; slam; feature-barcodes-tools; default; all`; unit tests) on a fresh export | `$D/logs/docker_seq.log` |
| `test_parameters_default_generation.py`, `test_htslib_build_discovery.py` (9 tests), MCP tests (589) | pass | — |
| CI partial builds, 1.10 | 11/11 PASS (core, core-static, core-htslib, host-api, release-companion-tools, star-feature-call, feature-barcodes-tools, process-features-lib, yremove-tools, core-external-htslib, core-portable) | `$G/partial/110/summary.tsv` |
| CI partial builds, v1.9.5.a baseline | 8/8 PASS (its default list) | `$G/partial/195/summary.tsv` |

`$D` = `/mnt/pikachu/star_suite_v1100_gates_20260928`; `$G` = `$D/gate`.
The original execution root was
`/tmp/claude-1000/-mnt-pikachu-chromap-suite-paper/53e97281-48e6-45f8-90f8-a4913b843536/scratchpad/gate`.
Logs retain that original path. The evidence copy preserves symlinks, not
their external targets; it is an archive, not a relocatable execution tree.
See `$D/preservation.tsv` and `$D/PRESERVATION_COMPLETE_UTC` for checksum-copy
verification. Scratch sources were not removed.

### G-S1 continuation

- Baseline completed at 18:17:21 UTC (`finished_utc`; `done 195` in
  `$D/logs/gs1_195.log`). All 24 production rows ran: 18 exited zero,
  six failed. All 13 Tier A tests passed. Do not restart this arm.
- Candidate started at 18:19:05 UTC with
  `nice -n 10 bash "$G/run_gs1.sh" 110`, under the driver's per-row shared
  host lock, and completed at 18:57:26 UTC. Production: 18 PASS, four FAIL,
  two SKIP; Tier A: 13/13 PASS. Inspect `$G/run110/manifest_status.tsv`,
  `$G/run110/tierA/status.tsv`, and `$G/run110/finished_utc`. The four failures
  also occur on the baseline. Do not restart either batch.
- The two SLAM determinism rows were held before launch because they repeat
  identical FASTQ executions internally. Only the gate-local candidate
  manifest gained a requirement for `$D/IDENTICAL_REPEATS_APPROVED`; that
  file does not exist. Explicit owner approval has been requested, not
  received. These are missing coverage, not successful tests.
- First PBMC capture: 22 non-log files identical, no differing/missing files.
  The PBMC 100K `report.json` inputs and all three profile QC/matrix-checksum
  objects match (vanilla, modern, modern_bam). This is selected comparison
  evidence, not a full G-S1 pass.
- Partial builds finished (both trees).
- Gate binaries (built in container `star-suite-gate-v1100:20260928`):
  `$G/src195/core/legacy/source/STAR.real` (1.9.5.a, `4824548`) and
  `$G/src110/core/legacy/source/STAR.real` (1.10.0, `0b92ce2`). Full source and
  executable SHA256s, image ID and compiler are in `$D/BINARY_PROVENANCE.json`.
  `STAR` is a capture wrapper; PBMC `report.json` hashes that wrapper, not
  the compiled executable. Use the `.real` hashes for binary provenance.
- 10M-read fixture for G-S3: `$D/gs3/pbmc10k_10M/L00{1,2}_R{1,2}.fastq.gz`
  (first 5M read pairs of each PBMC 10k v3 lane).

### Output comparisons (not a G-S1 pass)

- Read-only execution reports: `$D/gs1_baseline_execution_audit_v2.json` and
  `$D/gs1_candidate_execution_audit.json`. Both fail the audit; a zero exit
  from the batch driver is not a gate pass.
- `$D/tools/compare_gs1.sh` now requires finished markers, refuses report
  overwrites, pairs captures by unique normalized arguments rather than
  invocation index, and propagates failures. It returned 1. Of 108 baseline
  and 100 candidate captures, 83 paired uniquely (31 pairs fully identical
  under the strict comparator); 25/17 remain unmatched, including held
  calls and temporary-path differences. Across paired captures, 968 files
  match and 107 differ. The kept-output scan reports 2,042 matches, 196
  differences, 22 baseline-only paths and no candidate-only paths.
- These are triage counts, not counts of biological regressions. The strict
  comparator has incomplete Step 0 normalization, scans intermediates and Git
  metadata, and does not follow symlinks. Some log names, SAM headers and
  telemetry are not normalized. Baseline-only Chromap-adapter outputs belong
  to the deliberate Multiomics hand-over.
- `$D/selected_output_differences.json`: all four A375 feature MEX comparisons
  (raw/filtered, permits off/on) match by feature ID and barcode despite
  row/column ordering. All compared SAM bodies match after removing
  `@PG`/`@CO`. The upstream BINSEQ decoded R1 and R2 each contain the same
  24,893-record multiset; their order differs.
- UCSF `counts.h5ad` is **not byte-identical**: three `uns` source-path
  datasets differ. Its expression, Velocyto and observation dataset values
  match. The diagnostic compares HDF5 datasets, not attributes/layout.
- Downstream `final_counts.h5ad` and `unfiltered_counts.h5ad` have genuine
  differences in `obs/doublet`, `obs/doublet_scores`, and `obs/singlet`.
  Both arms call 243 doublets, but only 143 identities overlap (100 unique to
  each arm). Adaptive-MT summary fractions also differ. Other inspected
  H5AD dataset values match, apart from provenance paths.
- The smoke delegates to `morphic-recipes` at
  `e27f1745bf3d59aaf71cc0301b2f31845c759180`. Its
  `scripts/run_star_cell_doublets.R` calls scDblFinder without setting a seed.
  STAR's copy sets `SCDBLFINDER_SEED=1`, but is not the script this wrapper
  executes. This is a reproducibility gap consistent with the observed
  annotation differences, not evidence that STAR changed the count matrix.
  The external repository was not modified; no annotation difference is
  waived by this diagnosis.

## Earlier checkpoint (superseded by Latest follow-up above)

1. Finish disposition of the four original failures using the follow-up above;
   retain existing outputs and validate any changed code against a clean build.
   The execution auditor's 12 unit tests pass. Correct the affected test
   contracts with explicit coverage, not by accepting all baseline failures.
2. Resolve the canonical downstream recipe's seed gap in its own repository
   with owner coordination; pin recipe/image versions and validate changed
   downstream settings against the already-produced MEX. Complete capture
   pairing and Step 0 comparison normalization without relaxing biological
   comparisons. Do not silently ignore doublet differences.
3. G-S3 is on hold. The Flex caller flag is fixed, but its within-version
   hash-on/legacy matrix parity is still failing. Resolve that and obtain explicit
   owner approval for identical timing repeats. The question has been asked;
   no approval has been received. Do not run the existing driver unchanged:
   it overwrites runs and continues after failures. The 100K GEX workload has
   already run during G-S1, so timing it again also requires repeat approval.
4. Obtain approval before the held SLAM determinism repeats. Record approval
   provenance rather than creating the gate marker speculatively. Run only
   the missing cases, with preserved prior evidence.
5. Complete the output comparisons and disposition the failures before an RC
   tag. Keep release notes and distribution docs explicitly pending until
   gates are satisfied. Commit locally only; no tag, merge or push yet.

## Baseline failures requiring disposition

All findings below come from the clean, container-built baseline, not from
stale objects. Matching failures on 1.10 would establish that they predate
this refactor, not prove the affected workflow works.

| Row | Observed baseline failure | Required follow-up |
|---|---|---|
| `cbq-ocm-composite-local` | STAR exits 111: Velocyto's gene-like source did not populate per-read CB/UMI storage | Allocation guard fixed; follow-up validation above supersedes the initial diagnosis |
| `pf-dynamic-permit-100k` | Validator requires acquires == workUnits: baseline 7,264 vs 369,546; candidate 5,692 vs 369,546. Consumer permit acquisition is batched, not per record | Obsolete assertions replaced, nine validator tests and all four saved records pass. Later stress modes still lack execution coverage |
| `flex-hash-screen-replay` | Flat and tiered agree, but both differ from pinned truth on 2,024/800,000 reads, all Pass -> Deny (negative codes 6: 1,873; 4: 136; 8: 15) | Isolated to negative-cache policy change at da1b81d, not duplicate handling. Version the oracle without overwriting historical truth |
| `flex-hash-screen-100k` | Exits 102 before mapping: legacy expected-cell tuning with default tag-aware caller | Caller flag fixed; follow-up reveals the same hash-on/legacy matrix-parity failure in both versions |
| `slam-cbq-divergence` | Baseline only: FASTQ vs CBQ SAM body order differs; candidate held | Exact multiset is identical (113,516 records), and staged sorted-SAM, diagnostics and pre-NTR checks pass; decide whether record-order parity is required |
| `slam-cbq-divergence-pe` | Baseline only: same exact-order failure; candidate held | Staged sorted-SAM, junctions, diagnostics, transitions and pre-NTR checks pass |

Additional coverage notes:

- UCSF smoke exits zero by design when CellBender fails and writes a fallback
  H5AD (`tests/run_ucsf_corrected_production_100k_smoke.sh:112`). Both arms
  requested `--cuda`; the 100K fixture has only 5,584 observed barcodes.
  CellBender 0.3.2 fails during prior estimation with an IndexError and never
  reaches inference. This tests the fallback, not GPU denoising. The failure
  marker is under each arm's `ucsf/samples/EBs2_2/`
  `downstream_genefull_velocyto_cellbender/cellbender/`.
- Full baseline SLAM diagnostics are retained in
  `$D/slam_cbq_divergence_harness_20260928T181334Z` (SE) and
  `$D/slam_cbq_divergence_harness_20260928T181455Z` (PE).
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
