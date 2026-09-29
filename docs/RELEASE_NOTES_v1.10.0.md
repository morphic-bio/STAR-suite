# STAR Suite v1.10.0 Release Notes

Local release candidate: 2026-09-29 (`v1.10.0-rc1`; not published)

STAR Suite 1.10.0 removes every dependency on other suites. STAR no longer
links Chromap Suite or RapidMACS, and no longer hosts the joint RNA + ATAC
(multiome) run. That integration moves to Multiomics Suite, which builds the
multiome binary from released STAR Suite, Chromap Suite and RapidMACS
versions. In its place STAR gains a small, generic host interface that knows
nothing about Chromap or ATAC. The release also adds paired hashtag
demultiplexing and per-library ambient-FDR feature calling.

`STAR --version` reports `1.10.0`. Debian source packaging uses `1.10.0-1`;
Ubuntu packages use `1.10.0-1~ubuntu22.04.1` and `1.10.0-1~ubuntu24.04.1`.
Upstream STAR remains `2.7.11b`, genome-index compatibility remains `2.7.4a`,
and legacy compatibility remains `2.7.1a`. Existing indexes do not need
rebuilding.

This is an unreleased candidate. The final `v1.10.0` tag follows once Multiomics
Suite has built its binary against the host interface; a change the interface
needs would come as `v1.10.0-rc2`.

## Candidate limitations

- The rc1 source contains date/time macros in the PCG helper header and
  libem's hard-coded `-march=native`. Post-rc1 branch fixes replace the unused
  arbitrary seed with a fixed literal and remove host-CPU detection, preserving
  explicit scientific seeds. A new regression enforces both requirements.
  Clean strict STAR/host-library builds, EmptyDrops tests, host API checks and
  the single-thread TranscriptVB scatter/gather smoke pass. The same fix and
  compiler regression cover the unlinked old Flex PCG copy. The Multiomics
  integration gates and a new immutable candidate are still required before
  stable release; rc1 itself has not changed.
- The bundled official recipe catalog remains the pinned 1.9.5 snapshot.
  Its multiome recipes have not yet migrated to the Multiomics executable
  and must not be treated as supported standalone STAR 1.10 workflows.
  The validator confirms snapshot integrity, not runtime compatibility.
- The validated canonical downstream fixes are on
  `morphic-recipes` branch `fix/v1100-downstream-validation`, commit `bbf8b54`;
  their default-branch integration remains outstanding.
- The local tag and source handoff are not published binaries, installer
  bundles, Debian packages or container images. Those artifacts still need
  their release-pipeline build and runtime checks.

## Removed: the Chromap integration

- **The multiome run moves to Multiomics Suite.** Removed from STAR: the
  Chromap orchestration (asynchronous ATAC worker, permit-telemetry sampler,
  drain-time and saturation controllers), the libchromap contract library and
  its runners, the inline peak/matrix step, libscrna's ATAC evidence reader
  and its five multiome cell-calling tools, the CBQ-to-Chromap FASTQ adapter
  with its harness and smoke test, the multiome tests and scripts, the
  `morphic_multiome` MCP workflow and the integration runbooks.
- **57 parameters are gone:** the 42 `--chromapAtac*`, the 11
  `--multiomeAtac*`, and the ATAC-only permit parameters
  `--dynamicThreadAtacFloor`, `--dynamicThreadAtacController`,
  `--dynamicThreadAtacWorkEstimate` and `--dynamicThreadTelemetryIntervalSec`.
  STAR now rejects them as unknown parameters. Multiomics Suite accepts the
  same names, types and defaults.
- **Nothing is lost.** Every removed file, code block, symbol and parameter is
  listed with its last STAR commit in
  [the hand-over document](https://github.com/morphic-bio/STAR-suite/blob/v1.10.0-rc1/docs/HANDOVER_MULTIOMICS_1.10.md),
  so that Multiomics Suite can copy it unchanged. The CBQ adapter moves rather
  than being deleted.
- **What stays in STAR:** the CAT-ATAC guide arm (a feature-barcode arm that
  reads its barcode in the ATAC whitelist namespace), the shared thread-permit
  pool and its saturation controller, and libscrna's EmptyDrops, OrdMag and
  occupancy code.

## Build

- **`make core` no longer needs Chromap Suite.** It builds STAR with its
  bundled HTSlib, which is what `make core-portable` built before;
  `core-portable` remains as an alias. `WITH_CHROMAP` and `CHROMAP_SUITE_DIR`
  are no longer used.
- **Published artifacts were already Chromap-free.** Release tarballs were
  built with `core-static`, Debian packages with `core-portable`, and the
  Docker image with `STAR_WITH_CHROMAP=0`. Only a local `make core` linked
  Chromap. Users who built that way and write BAM will see different
  compressed BAM bytes with identical records, because STAR's bundled HTSlib
  replaces the system HTSlib.
- **`HTSLIB=external`** compiles and links STAR against an installed HTSlib
  (`pkg-config htslib`, or `HTSLIB_CFLAGS`/`HTSLIB_LIBS`), for programs that
  embed STAR alongside other HTSlib users.

## Host interface (new)

A program can now link STAR Suite as a library and run its own work beside
STAR in the same process, sharing STAR's thread permits. See the
[host interface reference](https://github.com/morphic-bio/STAR-suite/blob/v1.10.0-rc1/docs/HOST_API.md).

- `make star-host-lib` builds `libstar_suite.a` and `libstar_suite.link`
  (the libraries to link after it). The headers are in
  `core/legacy/source/host/`. Host API version 1.
- `star::host::runMain(argc, argv, hooks)` runs STAR. The `STAR` executable is
  `runMain(argc, argv, nullptr)`, so standalone behaviour is unchanged.
- Parameters STAR does not know go to the host, which accepts or rejects them.
  STAR records accepted ones in `Log.out` and in the command lines it writes to
  BAM headers, including those given in `--parametersFiles`.
- Callbacks run before read mapping (`preflight`, `start`) and after all of
  STAR's own work (`finish`), where STAR used to call Chromap. Further
  callbacks set the extra permit threads and the initial floors, refuse routes
  that cannot share the pool, and require a full pool at exit.
- The third permit domain, formerly `ATAC`, is now a neutral `EXTERNAL`
  domain that STAR lends to the host. Standalone permit logs name it
  `external`; a host supplies its own label (Multiomics Suite uses `atac`).
  The `atacController=` field of the "Dynamic thread interface enabled" line is
  gone.
- The saturation permit controller is now a public header,
  `host/SaturationPermitController.h`, in namespace `star::permits`.

## Feature calling

- **Paired hashtag demultiplexing.** `star_hash_demux_method=pair` in a
  pf-multi library assigns each cell to a sample from a
  `sample, hash_a, hash_b` table (`star_hash_sample_table`), with
  `star_hash_min_pair_ratio` (default 2.0). Cells whose top two tags form an
  unknown pair are reported as `unknown_pair`. The standalone
  `assignBarcodes` accepts `--hash-demux-method pair`. The default `ratio`
  method and its output are unchanged.
- **Per-library ambient-FDR calls.** `star_feature_caller=ambient-fdr` now
  calls a named non-GEX library on its own, with `star_feature_call_fdr`
  (default 0.01) and `star_feature_call_min_umi` (default 1), within the
  EmptyDrops cell set. Output is written to
  `outs/feature_analysis/<library_id>/ambient_fdr/`. CellTag libraries are
  selected by `feature_type=CellTag` and are not mixed with guides. The
  automatic CRISPR path and its `guide_*` outputs are unchanged.

## Correctness and regression coverage

- **Inline CB correction with Velocyto, without BAM.** A memory-saving
  allocation guard incorrectly skipped the per-read CB/UMI storage requested
  by Velocyto's gene-like source. Honor that existing requirement, including
  native OCM composite barcodes. No cell-calling or velocity-counting rule is
  changed. A new synthetic regression checks exact GeneFull and velocity
  layers in plain, inline-corrected and OCM modes.
- Update the feature-permit smoke to validate batched acquisitions and
  concurrent GEX work instead of expecting one permit per record. Add nine
  validator unit tests, including invalid and missing telemetry.
- Explicitly select the legacy caller in the legacy Flex smoke. Its remaining
  hash-on/legacy parity failure remains a historical diagnostic, not the modern
  release gate. The release gate now uses the established H0/H1X2 half-probe
  fixture with matching model/filtered gene axes and preserved 1.9.5.a output.
- Version the standalone March decision replay with an explicit legacy negative
  policy, leaving the default policy and STAR runtime classification unchanged.
- Seed scDblFinder in the canonical downstream recipe's R container and match
  the STAR mirror. Fail on doublet errors rather than silently assigning singlets.
  Require CUDA and successful CellBender output; reject stale/partial failures.

## Validation

The local STAR gates G-S1, G-S2 and G-S3 are accepted. Canonical recipe
integration and Multiomics integration/snapshot migration remain stable-release
dependencies. The owner subsequently authorized starting version creation;
rc1 is local only. Publication and the stable master merge remain separate.

**September 29 update:** the owner authorized repeats. Fresh-build G-S2 and
15/15 Tier A checks pass; SLAM SE/PE repeat and cross-version results match
when unsorted SAM record order is treated separately. Seeded downstream
repeats and cross-version H5ADs match, including every dataset and attribute
after explicit provenance-path normalization. The UCSF reused-fixture test
now regenerates its configuration with current paths and rejects zero-match
GEX filtering; its repaired MEX outputs match the previous run exactly.

The expanded audit found multi-thread TranscriptVB quantification variation
on the pristine baseline too. The owner selected exact single-thread validation:
all treated/no-4sU SE/PE tables and ordered evidence payloads match across
versions. No numerical tolerance or record sorting was used. The kept-output
audit is closed and **G-S1 is accepted**. All eighteen G-S3 attempts match
outputs and median timing/RSS changes are within limits. **G-S3 is accepted**
under the subsequent owner clarification that host activity is diagnostic
within +3% wall / +2% peak RSS tolerances. Maximum observed median increases
are 0.151% wall time and 0.035% peak RSS across the three workloads. All six
host flags remain; separately verified wrapper errors are documented rather
than rewritten as successful exits. No benchmarks were rerun for this
reassessment. Deterministic multi-thread TranscriptVB is not claimed.
Details:
[validation follow-up](VALIDATION_STAR_1_10_0_20260929.md).
The earlier evidence below is retained as history, not the current hold reason.

- Previously completed: clean no-Chromap build, bundled/external HTSlib host
  API tests, and partial builds (baseline 8/8; candidate 11/11).
- Both G-S1 arms completed. Tier A passed 13/13 on each. Production results:
  baseline 18 PASS / 6 FAIL; candidate 18 PASS / 4 FAIL / 2 SKIP. The four
  candidate failures also occur on the baseline; two SLAM determinism tests
  remain held for repeat approval. These results do not satisfy G-S1.
- Selected PBMC matrices and keyed A375 feature counts match. UCSF
  `counts.h5ad` dataset values differ only in provenance paths. Downstream
  doublet identities/scores differ; the external canonical recipe does not
  set the seed used in STAR's local R script. That gap is now fixed in the
  isolated recipe branch, with one successful saved-MEX downstream run. A
  repeated-run reproducibility check is still held for approval.
- Both UCSF rows continued after CellBender failed during prior estimation;
  the successful smoke wrapper exits validate fallback, not GPU denoising.
  The updated 100K smoke separates this sparse-input test from CUDA validation.
  A real raw-droplet CUDA smoke has since passed inference and H5AD layer
  integration on 20,000 droplets. It used five epochs, not production convergence.
- Follow-up: the clean-built OCM allocation fix passes six synthetic cases
  (three modes, normal and low-memory). Both real 1,000-read OCM executions
  complete with byte-identical CBQ/FASTQ pooled and per-sample outputs; all
  eight ordinary-GEX count configurations pass. The OCM wrapper's subsequent
  stale CBQ/Y-removal rejection check was corrected and validated on the
  saved plans, not by repeating alignment. Its original exit 1 is retained.
  Permit-validator units pass, and all
  four saved off/on telemetry records pass the corrected assertions; the
  resizing, shadow/active-controller and forced-exit/recovery modes have since
  passed. The expanded
  Tier A suite has 15 cases; the initial batch's 13/13 is historical evidence.
- Historical Flex replay differences are reproduced by the September 4 negative-cache
  policy change using exactly the same cache and March dump, not a change of
  reference or duplicate handling. The historical oracle is retained and passes
  all 800,000 decisions under its explicit historical policy.
- Current Flex 100K: **121/121 non-log, non-diagnostic outputs byte-identical**
  to preserved 1.9.5.a output with the correct H0/H1X2 half-probe fixture.
  The earlier March H0/H1 test was the wrong modern gate; its matrix drift is
  not evidence of a regression in the current half-probe route.
- G-S3 was initially held for repeat approval and the wider output audit.
  Both subsequently completed; current acceptance is summarized above. Its
  driver refuses overwrites, propagates failures, and selects the modern
  half-probe workload.
- Official snapshot digest/count validation passed (11 recipes, 10 evidence
  records). Migration of the pinned multiome recipes is still a separate
  dependency.

See the [runbook](runbooks/RUNBOOK_STAR_1_10_0_HOST_API_20260928.md) and
[handoff](handoffs/HANDOFF_STAR_1_10_0_HOST_API_20260928.md) for gate evidence
and remaining work. The observed local medians meet the regression limits;
they are not a claim of noise-free measurements or a speedup.
