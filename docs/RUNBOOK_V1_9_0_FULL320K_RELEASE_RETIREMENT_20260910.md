# STAR Suite 1.9.0: full 320K benchmarks, release, and instance retirement

**Status: COMPLETE on 2026-09-10. Release published; both full benchmarks validated; archive verified; instance stopped, snapshotted, and terminated.**

Authorization: “Go ahead with the runbook and the archive job.” The review gate is open.

The user requested this reviewable runbook before committing, merging, pushing,
delegating the release, running benchmarks, or retiring the instance. The user subsequently opened that gate with the approval recorded above.
Carry the approved sequence through without requesting the same authorization
again.

## Required outcome and order

1. Audit, commit, merge, and push the completed STAR changes.
2. Assign a subagent to prepare and deliver release **1.9.0** in an isolated
   checkout. Freeze the release's compiled source and version metadata.
3. In parallel with the subagent's release packaging and checks, benchmark that
   frozen source on the **entire 320K dataset**, first CBQ, then FASTQ-BGZF.
4. Validate both runs and save the important new and existing SSD files to S3.
   Verify the completed archive. Release publication continues independently
   on GitHub/local packaging infrastructure, with no EC2 dependency.
5. Stop the instance and wait until EC2 reports `stopped`.
6. Snapshot its EBS volume(s) and wait until every snapshot reports `completed`.
7. Terminate the instance and verify EC2 reports `terminated`.
8. Verify the published release assets, mirror them and final local records to
   S3, and close the overall runbook only when both release and retirement finish.

The release can be prepared while the benchmarks run. Publish the final
immutable `v1.9.0` tag only after its checks and the full-set validation pass.
The two cloud benchmarks must run sequentially, with no archive compression,
checksum sweep, upload, build, or other benchmark competing on that instance.

## Facts and resources

| Item | Value |
| --- | --- |
| Repository | `/mnt/pikachu/STAR-suite` |
| Original integration base commit | `4c33c0145c7bd76b2c6f5a915da0bcf0c2fc4998` |
| Public remote | `origin`: `git@github.com:morphic-bio/STAR-suite.git` |
| Integration destination | `master`, using `git merge --no-ff` |
| Release candidate branch | `dev-release-v1.9.0`, subject to checking its remote state |
| Release tag / binary version | `v1.9.0` / `1.9.0` |
| Debian source version | `1.9.0-1` |
| AWS profile / region | `uw` / `us-west-2` |
| Instance | `i-06de289faa5d78117` |
| Instance type | `m6id.12xlarge`, 48 vCPUs |
| Original root EBS volume (deleted on termination) | `vol-059ae86ad14651d0b`, 30 GiB, `/dev/sda1` |
| Local SSD storage | `/scratch`, RAID0 across two instance-store NVMe devices |
| Full dataset | 4 lanes, L001–L004, **7,303,142,230 read pairs** |
| Grouping | 8 samples, 16 probe tags; retain the existing fused-tag mapping |
| Work/results S3 bucket | `star-suite-320k-benchmark-alt-171440768238-us-west-2-20260904` |
| Original backup S3 bucket | `star-suite-320k-scratch-backup-171440768238-us-west-2-20260903` |

Refresh instance, volume, branch, and tag state before acting. Never substitute
another instance, overwrite an existing release tag, or force-push a shared
branch to make these names fit.

### Existing results are controls, not the requested new benchmarks

| Dataset / build | CBQ | FASTQ-BGZF |
| --- | ---: | ---: |
| Full 320K, September 9 all-feature caller | **1,217.150 s (20m17s)** | **1,452.071 s (24m12s)** |
| L004 only, latest optimized caller | **201.901 s (3m22s)** | Not measured with that build |
| L004 only, earlier half-khash build | 522.278 s | 1,030.833 s (17m11s) |

The 3m22s result processed only **1,823,648,323** read pairs. It cannot satisfy
the full-set requirement or be substituted into the full-set table.

Saved Cyto 0.4.7 controls on this same full dataset and 48-thread instance are
758.46 s (12m38s) for CBQ and 1,191.57 s (19m52s) for **ordinary gzip FASTQ**.
Retain them as measured controls; the requested new executions are the two
STAR full-set arms. Do not label Cyto's gzip timing as BGZF or rerun Cyto
implicitly.

## Historical review pause and archive state

At the initial review pause, no commits, merges, pushes, release delegation,
new full-set benchmarks, instance stops, snapshots, or termination had been
performed under this runbook. The execution ledger below records subsequent
authorized work.

An archive job launched under the earlier instruction was still `InProgress`
when this document was prepared:

- SSM command: `253d1a2e-6b11-4d99-8528-5bb53b7e5dca`.
- Local controller:
  `/mnt/pikachu/star_suite_paper/analysis/instance_retirement_20260910/archive_remote.py`.
- Cloud state: `/scratch/instance_retirement_20260910/`.
- S3 prefix: `instance-retirement/20260910/` in the work bucket.

**Review pause confirmed at 2026-09-10 07:54:09 UTC.** Archive Python PID
`552001` and its children `554486`, `554487`, and `554488` were suspended;
each reported process state `T`. The instance remains running. The pause
record, including process start ticks to guard against PID reuse, is at
`/mnt/pikachu/star_suite_paper/analysis/instance_retirement_20260910/review_pause.json`.
This suspension honors the user's review pause; it does not start an execution
phase of this runbook.

After user approval, the validated process identities were resumed at
07:55:47 UTC. The archive completed and all 90 bundles were independently
verified before full-set timing began. Those historical PIDs must not be
signalled again. Shutdown remains gated on the supplemental archive and
remaining release/validation phases below.

The previous archive script inventories the old SSD tree once. Even a valid
completion marker from it cannot cover the new full-set runs. A final inventory
and supplemental archive after the new runs are mandatory.

## 1. Commit, merge, and push the implementation

Owner: primary agent.

1. Save the current Git status, base commit, scoped diff, untracked source
   inventory, and source SHA-256 manifest outside the repository. Preserve
   unrelated local edits and artifacts. Use an isolated checkout for
   integration; do not use `git add -A`, destructive resets, or broad cleanup
   in the dirty primary checkout.
2. Audit the final changes against the validated source bundles and reports.
   Include the cache loading/storage and khash implementation, paired 25-base
   half-probe lookup, caller parallelism, sort-once OrdMag, shared thread
   permits, required build dependencies, meaningful tests, and their docs.
   Audit the reader/caller integration dependencies already present in the
   base commit. Include only relevant hunks in mixed-purpose files.
3. Confirm obsolete experiments, large matrices, generated binaries, plots,
   copied repositories, and source-derived Cell Ranger investigation material
   are excluded from the commit. Preserve the existing TODO for ambiguous
   same-sample tag matches without implementing that separate experiment.
4. Create a feature branch from the current remote integration base, commit
   the scoped implementation, and merge it without squashing into the release
   integration branch. Run the relevant clean-build tests and required CI.
5. Merge the accepted integration branch into `master` with `--no-ff`, push to
   `origin`, and verify the remote commit and CI outcome. Preserve the DAG and
   record the commit IDs. Keep the primary checkout on `master`; reconcile
   only integrated changes while preserving unrelated dirty work.

Follow `AGENTS.md` and `AGENTS.local.md`. All Cell Ranger work remains clean
room: documentation, help, logs, matrices, BAMs, and HDF5 results are allowed;
Cell Ranger code, modules, source archives, mirrors, and source-derived notes
are not evidence for implementation.

## 2. Delegate release 1.9.0 and freeze the benchmark source

Owner: a release subagent, explicitly requested by the user. Do not spawn it
until the review gate opens. The primary agent coordinates merges and the
benchmark handoff so the two agents do not edit the same checkout concurrently.

The subagent's assignment is to complete the release, including:

- Work in an isolated release checkout based on the pushed implementation.
- Check whether `v1.9.0` or its candidate branch already exists remotely.
- Update `core/legacy/source/VERSION`, `debian/changelog`, release notes, and
  other canonical packaging metadata as required by the existing release
  tooling. Preserve upstream and genome-index compatibility versions unless
  an actual compatibility change requires otherwise.
- Prepare accurate release notes covering the shipped changes since v1.8.4,
  with measured results labelled by dataset, source, and input format.
- Commit and hand off the versioned source before the primary agent builds
  the benchmark binary. Record its commit and compiled-source manifest.
- Run the existing packaging and release checks in
  `.github/workflows/release.yml`, `docs/Star-binary-distribution.md`, and
  `scripts/release/`. Validate supported tarballs, installer, Debian packages,
  checksums, expected version, and smoke tests as the workflow requires.
- Prepare publication concurrently with the primary agent's full-set runs.
  After the shared acceptance gate, coordinate the final `--no-ff` merge to
  `master`, push the immutable `v1.9.0` tag, verify the release workflow and
  published assets, and report the release URL and asset checksums.

The primary agent clean-builds the Chromap-enabled production binary from the
frozen versioned source. Both full-set runs use the **same binary SHA-256**.
Packaging can produce platform-specific binaries, but compiled source must
match the accepted release. Record the benchmark build commit and final tag
relationship explicitly; a documentation-only follow-up does not silently
change the reported benchmark build identity.

If the subagent must change compiled source after the freeze, notify the
primary agent before doing so. Rebuild and validate that changed source; do
not describe earlier timings as measurements of the changed binary.

## 3. Run the two full-set benchmarks

Owner: primary agent, using the cloud instance above.

### Inputs and configuration

- CBQ: `/scratch/cbq/lane_000.cbq` through `lane_003.cbq`, all four together.
- FASTQ-BGZF: all eight files
  `/scratch/fastqs/16-plex_GEM-X_FLEX_S1_L00{1,2,3,4}_R{1,2}_001.fastq.gz`.
  Verify actual BGZF member structure; the `.gz` suffix alone proves nothing.
  Preserve the existing STAR R2/R1 argument ordering and mate pairing.
- Compact half cache:
  `/scratch/flex_half_production_20260910_v1/model.half.khash`.
  Expected SHA-256:
  `801fc6143a383b6a94c55307f816bf824c7a9c685d3049405b2ffcc4d72db296`.
- Paired model/export config:
  `/scratch/full320k_star_deprecated_20260909_v1/config/`.
- Empty genome directory:
  `/scratch/full320k_star_deprecated_20260909_v1/empty_genome_index`.

Use the completed September 9 full-set `argv.json` files as the argument
baseline and audit every intentional difference: the new binary, compact
stored half cache, new output directories, and the newly integrated execution
optimizations. Preserve the statistical and counting settings.

| Setting | Required value / behavior |
| --- | --- |
| STAR threads / total shared permits | 48 |
| Model features / exported features | 19,068 / 18,129 |
| Deprecated features | Retained for modeling, excluded by the filtered-export allowlist |
| Caller | Tag-aware, grouped sample tags, joint occupancy after group calling |
| Simulations / BH FDR | 100,000 / 0.01, verify resolved settings in logs |
| Rank/tie policy and floor | Existing deterministic quality ties and model-count floor |
| Alignment | `--flexNoAlign 1`; no alignment |
| Reference index | No loading; empty genome directory; verify no-genome activation |
| BAM / SAM | `--outSAMtype None` |
| Per-read decision sidecar | `--soloFlexDecisionSidecar -` |
| Caller diagnostics | Enabled, matching the baseline |
| CBQ ingestion | `Binseq PE`, range mode, all four lanes |
| BGZF ingestion | Range mode, 48 reader-thread limit, CRC checks enabled |
| Bucket configuration | Preserve full-set baseline auto mode and 32 GiB setting |
| Input conversion | Outside timed runtime; use already staged inputs |

The full baseline sets `--soloFlexEdNiters 0` and resolves to 100,000
simulations. Confirm this in the new binary's log rather than treating `0`
as the simulation count.

### Execution and timing

1. Verify the release source/binary identity, input manifest, model/cache
   pairing, free disk and memory, fixtures, and absence of competing work.
   Complete required changed-build tests; do not repeat identical fixtures
   already executed against the exact same binary and environment.
2. Allocate a fresh cloud run root, proposed
   `/scratch/full320k_v190_20260910_v1/`, with separate `full_cbq/` and
   `full_bgzf/` directories. Never reuse an existing timed output directory.
3. Save explicit argument arrays, shell-readable commands, environment,
   binary/source/input hashes, host configuration, and start timestamps.
4. Run **CBQ once**, then **FASTQ-BGZF once**. Use the same frozen binary,
   48-thread host, and statistical settings. Clear the OS page cache before
   each timed run. Time the complete STAR process with monotonic wall time
   and GNU time, including cache loading, ingestion, counting, cell calling,
   occupancy, matrix output, and process cleanup.
5. Write a wrapper completion record only after successful process exit.
   Assert **7,303,142,230 input read pairs in each arm**. A completed L004 run,
   a partial progress log, or an old completion file cannot pass this gate.
6. Record wall time, CPU time, peak RSS, reader and caller phases, permit
   telemetry, model/export feature counts, read classifications, called cells,
   output size, and any spill activity. Do not compress/archive concurrently
   with the second timed arm; defer heavy result processing until both finish.

This authorizes one new execution of each arm against the selected build.
An absolutely identical rerun still requires explicit repeat authorization.
A source fix followed by a clean rebuild is a different execution. Retain
failed attempts and their status instead of overwriting them.

## 4. Validate and report before release/retirement

- Compare all raw matrix coordinates and all eight filtered matrices between
  formats and against the September 9 full-set controls, using canonical
  axes where needed. Check exact called-cell identities, not just totals.
- Expected old full-set reference values: 333,439 final STAR cells,
  7,924,268 raw barcodes, 725,480,265 nonzero coordinates, and
  1,310,330,528 raw UMIs. Investigate any difference; these are controls,
  not values to force into the new results.
- Compare OrdMag outputs, candidates, ambient profiles, p-values, rescue
  decisions, occupancy removals, and sample summaries. Execution metadata
  such as worker counts may differ; statistical differences need explanation.
- Verify no reference loading, BAM, SAM, or per-read sidecar, and that the
  shared caller permit budget does not exceed 48 and is fully returned.
- Reuse the existing CR/Cyto concordance only if exact count and callset parity
  establishes it. Otherwise recompute the relevant metrics from saved results
  and isolate the cause before claiming unchanged concordance.
- Publish a table comparing the new full-set times with **1,217.150 s CBQ**
  and **1,452.071 s BGZF**, with time ratios, memory, binary identity, and
  validation status. Keep L004 and full-set tables separate. Mark the two
  executions as single-run measurements, not repeated-trial estimates.

Once benchmark validation and release checks pass, complete the release
merge/tag/publication and verify the published artifacts. If a correctness
regression or missing SSD artifact remains unresolved, retain the instance.
Independent hosted publication can continue after verified SSD archival and
retirement; do not declare the overall sequence complete until it also passes.

## 5. Finish and verify S3 preservation

Resume the suspended archive safely after timing, or use its recorded partial
state to finish preservation. Verify its actual status before resuming: SSM
may time out during a long review pause. Never submit the old archive script
blindly or mistake a stale completion marker for coverage of later files.
Before sending a continuation signal, match each PID's saved start ticks and
parent identity; do not signal a reused PID. Resume surviving descendants and
then the controlling archive process. Record the result.

Create a supplemental archive for the new run root and any files created since
the earlier inventory. Proposed prefix in the work bucket:
`analysis-tools/full320k_v190_20260910_v1/`.

Preserve:

- Full raw and filtered MEX outputs, reusable binary/HDF5 count caches where
  already available, and caller diagnostics for both new runs.
- Arguments, input/config/model manifests, runtime and validation records,
  stdout/stderr, completion markers, parity results, and the final report.
- The exact benchmark binary, source commit/manifest, relevant source bundle,
  compact khash and paired gene lists; include release provenance and checksums.
- Existing important STAR and CR matrices, BAM/HDF5 oracle files, Cyto results,
  read audits, and scripts needed for future analysis. Preserve existing
  durable copies by verified manifest references rather than duplicating all
  objects. Keep excluded CR code/source-derived material outside clean-room
  evidence.
- This runbook, updated execution ledger, restore notes, original instance and
  volume descriptions, and final lifecycle API records.

Verify S3 object existence and sizes; record SHA-256 for new archive bundles
and compare their saved metadata. Existing objects referenced by the retirement
manifest must retain their verified content identities. Check a representative
restore of the new small metadata/diagnostics bundle before shutdown. Record
all dependencies on the two original buckets; delta archives alone do not
constitute a complete restore.

Write **`CLOUD_ARCHIVE_COMPLETE.json` only after the new full-set results
are included**. Save it locally and to S3 before stopping, because instance-store
data can be lost on stop. `FINAL_ARCHIVE_COMPLETE.json` adds the independently
published release assets and final local records once publication completes.

## 6. Stop → snapshot → terminate

Owner: primary agent, after benchmark acceptance and verified cloud archival.

1. Verify no required cloud process or upload remains, the cloud archive is complete,
   and the exact instance ID still matches this runbook. Save its current
   attached-volume list locally and to S3.
2. Request a normal stop of `i-06de289faa5d78117`. Record the API response and
   poll until its state is `stopped`. Do not force-stop by default.
3. Create a snapshot of every attached EBS volume identified at that point.
   Currently the only volume is `vol-059ae86ad14651d0b`. Tag snapshots with
   project, source instance, date, and archive prefix. Record every snapshot
   ID and wait for each to reach **`completed`**. A `pending` snapshot does
   not pass this gate.
4. Verify completed snapshot IDs and their source-volume identities. Store
   this evidence locally and in S3 before terminating the instance.
5. Terminate the exact instance. Poll until EC2 reports **`terminated`**.
   Confirm the snapshots remain available and record remaining volumes, if any;
   do not delete snapshots or unrelated resources.
6. Write and upload `instance_retirement_complete.json` containing the cloud state.
   After release verification, write combined `retirement_complete.json` with
   snapshot IDs, archive locations, release URL/tag/commit, benchmark results,
   and UTC timestamps. Update this runbook's ledger and give the user the
   final timing table, release link, archive location, and snapshot ID(s).

**The EBS snapshot covers the root disk, not `/scratch`.** The SSD tree is
restored from S3. Never stop the instance first and expect a later EBS snapshot
to preserve its NVMe data.

## Execution ledger

Update after each phase; do not mark queued work complete.

| Gate / phase | Status | Evidence to fill in |
| --- | --- | --- |
| User review of this runbook | **APPROVED** | User: “Go ahead with the runbook and the archive job.” |
| Earlier archive suspended or finished | **COMPLETE** | 90 bundles / 126,765,224,602 bytes independently verified; existing-object references also retained. Both new benchmark arms are covered by the verified supplemental archive. |
| Scoped implementation committed | **DONE** | `60fbe8cf6425596f10c1008ea5db7bd6c44ccf2e`; 47 scoped files matched tested manifests |
| Integration merged and pushed | **DONE** | Release implementation and notes merged/pushed; tag at `fc59b1b20b9fd34c1176acebd7e1688d4a82c1ec`; compiled-ancestor CI and release checks passed |
| Release subagent assigned | **COMPLETE** | `release_v190`; published release, verified downloads, and both container platforms |
| Versioned source frozen and clean-built | **DONE** | Build `a6c1c546cbb53403c542bb9819706631801c935b`; STAR 1.9.0 SHA `54311642da84160e2dd75321d1ea6306ef0a5724b32727b402216748a5707d67`; unit/cache tests and both input fixtures passed |
| New full CBQ benchmark | **COMPLETE; exact parity passed** | Exit 0; 7,303,142,230 pairs; 628.033 s; 103.968 GiB; 333,439 cells |
| New full FASTQ-BGZF benchmark | **COMPLETE; exact parity passed** | Exit 0; 7,303,142,230 pairs; 867.680 s; 93.794 GiB; 333,439 cells |
| Count/caller validation | **PASSED** | All 18 raw/filtered matrices, 114 caller diagnostic files, sample summaries and read classifications exact; 333,439 cells in both arms |
| v1.9.0 published and verified | **COMPLETE** | [v1.9.0](https://github.com/morphic-bio/STAR-suite/releases/tag/v1.9.0); tag `fc59b1b20b9fd34c1176acebd7e1688d4a82c1ec`; published checksums/versions/source commits and both container platforms verified |
| Final S3 archive includes new runs | **COMPLETE** | Original 90 + supplemental 3 bundles; restore test passed; 5,520 existing object references reverified; local records and 20 published assets also saved |
| Instance stopped | **COMPLETE** | Normal EC2 stop confirmed; `lifecycle/instance_stopped.json` |
| EBS snapshot(s) completed | **COMPLETE** | `snap-0428a020b3a0fbbd4`; completed before termination |
| Instance terminated | **COMPLETE** | `i-06de289faa5d78117` verified `terminated`; completed snapshots retained |
| Final report and restore records saved | **COMPLETE** | Local retirement records plus `s3://star-suite-320k-benchmark-alt-171440768238-us-west-2-20260904/instance-retirement/20260910/` |

## Supporting records

- [OrdMag and permits L004 benchmark](benchmarks/FLEX_ORDMAG_SORT_PERMITS_L004_20260910.md)
- [Half-khash L004 benchmark](benchmarks/FLEX_HALF_KHASH_L004_20260910.md)
- [September 9 full-set caller integration](HANDOFF_FULL320K_STAR_DEPRECATED_20260909.md)
- [Release packaging](Star-binary-distribution.md)
- [GitHub Actions](Github-actions.md)
- Local runtime/retirement records:
  `/mnt/pikachu/star_suite_paper/analysis/instance_retirement_20260910/`.
- September 9 full-set wrappers and completed evidence:
  `/mnt/pikachu/star_suite_paper/analysis/full320k_star_deprecated_20260909/`.

### Execution note: release wrapper compatibility

The prior master and first integration candidate failed the Flex tiny-public
smoke because the recipe wrapper passed a legacy fixed expected-cell override
to the default tag-aware caller. Release commit `1b7aa4e08c1da5e8d3a975893a06afd0983705fd`
selects the tag-aware caller explicitly and removes that incompatible override.
Its focused wrapper test and official snapshot validation passed. Compiled
source stayed identical to the benchmark freeze.

The local grouped unit test initially lacked its required artifact-directory
argument and stopped before executing its fixture; the corrected invocation
passed. The CBQ fixture comparison initially assumed nine matrices, but the
small fixture emits three. Existing CBQ output was validated against those
three control paths, then BGZF executed once and matched. Neither STAR fixture
was repeated. Corrections and original failure records are preserved in
`analysis/full320k_v190_20260910/`.

### Execution note: generated parameter defaults

Source-archive packaging exposed a stale tracked `parametersDefault.xxd`
whose timestamp matched its newer authoritative input. Commit
`91fdfd52cf5f629f4b61990291b433ff97d5f87a` refreshes that generated file and
content-checks it atomically during builds. The refreshed bytes have SHA-256
`033bba5177e252ebfd8bcf6a4067bd6a4eb14e455670af71d11eb6b86a2886d4`, exactly
matching the clean full-benchmark build. Effective compiled defaults are
unchanged; failed package artifacts and the corrected build checks are retained.

Master CI `34453932378` passed build/Tier A and secret checks. Its image
publication encountered Docker Hub HTTP 502 during login, before any image
build. The packaging-fix candidate passed its checks; final artifact and container
publication were verified separately as recorded below.

### Full-set benchmark acceptance

CBQ completed in 628.033 seconds (103.968 GiB peak RSS), versus the
1,217.150-second control. FASTQ-BGZF completed in 867.680 seconds (93.794 GiB),
versus 1,452.071 seconds. Each processed all 7,303,142,230 read pairs. Every raw
and filtered matrix, caller statistical diagnostic, sample summary and read
classification matched exactly. See the [full report](benchmarks/FULL320K_V190_20260910.md).
The validated 1.9.0 binary was atomically installed locally; the previous
binary is backed up in `analysis/full320k_v190_20260910/`.

### Execution scheduling update: independent release publication

After full-set validation and verified SSD archival, the remaining release
pipeline consists of hosted package/image builds and local download verification.
The release agent confirmed these have no dependency on the EC2 instance.
To avoid leaving the instance idle through those builds, instance retirement
now runs alongside publication. No release acceptance check is waived.
`CLOUD_ARCHIVE_COMPLETE.json` gates the stop and preserves all instance data;
`instance_retirement_complete.json` records cloud completion. Final published
assets and local records are mirrored afterward, then `FINAL_ARCHIVE_COMPLETE.json`
and `retirement_complete.json` close the combined task. The required cloud
order remains normal stop, completed EBS snapshot(s), then termination.

### Release asset filename compatibility

GitHub stored the six Debian asset filenames with `.ubuntu` where the CI
files used `~ubuntu`. Package contents and internal Debian versions stayed
unchanged. The public `SHA256SUMS` was corrected to the actual public names;
all 20 downloaded release assets were verified, and the 19 payload files
matched their accepted CI artifact hashes. The original checksum manifest
and `public_filename_mapping.json` are preserved in the release records.

Follow-up for the release workflow: generate checksum entries using GitHub's
normalized public filenames before upload. This operational correction did
not change the immutable v1.9.0 source tag.

### Container publication

The final multiarchitecture image was built with the existing local
`rmbuilder` from a fresh archive of the immutable release commit
`fc59b1b20b9fd34c1176acebd7e1688d4a82c1ec`. The image build and its checks
succeeded for `linux/amd64` and `linux/arm64`. The accepted image index is
`sha256:9c744f7a83ca96d76dc72d44c76ce3c4e1869cd752aaacaa6da07080e1ef0b99`.

After the staging image passed registry verification, the redundant hosted
image workflows `34455225261`, `34456618393`, and `34456680980` were canceled
and reached terminal states before promotion to `biodepot/star-suite:v1.9.0`
and `latest`. Their completed build/Tier A or artifact checks passed; the
canceled image stages are not reported as successful hosted executions.
The local build, guarded promotion, registry platforms, and binary identities
are recorded in `local_image_build.json`, `image_publication_selection.json`,
and `image_verification.json` under the archived release provenance.

### Retirement completion

Verified terminal state: `terminated`. Completed retained snapshot(s): `snap-0428a020b3a0fbbd4`.
Root EBS was stopped before snapshotting. `/scratch` is preserved through the
verified original and supplemental S3 manifests, not by the EBS snapshots.
The published release assets, source/build records, exact matrices, and
restore/lifecycle evidence are retained. No further instance work remains.
