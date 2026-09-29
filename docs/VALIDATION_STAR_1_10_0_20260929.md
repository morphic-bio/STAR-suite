# STAR Suite 1.10.0 validation follow-up

Date: 2026-09-29 UTC. **Local STAR gates G-S1, G-S2 and G-S3 accepted.** The owner
authorized the remaining repeats with "Ok finish the validation gates" and then
directed: "Do it on a single thread with a smaller set if needed - the order
matters unfortunately". The exact ordered controls below close the TranscriptVB
hold. G-S3 is accepted under the subsequent owner instruction to treat host
activity as diagnostic when measured runtime and memory stay within tolerance.
No benchmarks were rerun for this reassessment. This is not stable-release
acceptance. No tag, push, or master merge was performed.

## Build and provenance

Artifact root `D=/mnt/pikachu/star_suite_v1100_gates_20260928`;
new evidence `A=$D/acceptance_20260928`.

- Fresh candidate checkout/build: `a9c16362fcbc14119c4318a65ee031dda03068cf`,
  STAR 1.10.0, `make core` plus companion tools.
- Candidate executable: `$A/src110/core/legacy/source/STAR`, SHA256
  `bfb53b3f3227a37f5abc38afe1b8efd5ac7e617a8dec5b9e2b8f565be4c41e4c`.
- Baseline: preserved clean-container STAR 1.9.5.a build at `4824548`,
  code-equivalent to published `f0d9f27`; SHA256
  `066734e7c100757f28101567faf41e964823def6a7d599d2d2bf243949cdcdcf`.
- Both built with gate image
  `sha256:598e763f8241939b051732858428a4c42618135ace91186e5750d5806b8218aa`
  (Ubuntu 22.04 / g++ 12.3.0). New host-library build did not change the
  candidate executable hash.
- Canonical downstream recipe: isolated branch
  `fix/v1100-downstream-validation`, commit
  `bbf8b5408c0f0c706ff4a06ac3025de4f0ebce62`, selected explicitly with
  `MORPHIC_RECIPES_ROOT=/mnt/pikachu/morphic-recipes-v1100-20260928`.
  It has not been merged into the canonical default branch.
- Downstream R image pinned to
  `sha256:c052c1568727f24a6e2c8ec793357ccfdf1ffdf5622542a7b2da40d8df47825b`;
  `SCDBLFINDER_SEED=1`. Python backend `host` for saved-MEX comparisons.
- Heavy executions were serialized with
  `/mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock`, nice 10.

## Completed checks

| Check | Result | Evidence under A |
|---|---|---|
| Fresh core and companion builds | PASS | `build.log`, `BUILD_COMPLETE`, `STAR.sha256` |
| G-S2 host API | PASS: 180 files for each no-callback host; 20 external-host files; 400 external permits; parameter/lifecycle and permit tests | `host_api.log`, `host_api/`, `G_S2_COMPLETE` |
| Expanded Tier A | 15/15 PASS after correcting companion executable paths in the driver | `followup/status.tsv`, `companion_checks/status.tsv` |
| Paired-hashtag demux, CellTag ambient FDR, default-parameter generation | PASS | `companion_checks/` |
| SLAM SE and PE repeat/thread/CBQ checks | Each: 21 PASS, one ORDER_ONLY, zero FAIL | `followup/slam_{se,pe}/stage_status.tsv` |
| SLAM cross-version comparison | 40/40 exact: quantification, diagnostics, transitions, junctions and SAM record multisets across four arms per mode | `slam_cross_version.json` |
| Seeded downstream reproducibility | All six H5ADs: all 396 datasets and 974 attributes match; seed, doublets, scores and normalized QC JSON match | `followup_output_parity_v2.json` |
| Seeded downstream cross-version | All six H5ADs match, normalizing only explicit input/output provenance roots | `seeded_baseline_downstream/`, `followup_output_parity_v2.json` |
| Repaired UCSF fixture | 3,943 cells; all 35 checked MEX matrix/axis files identical to prior candidate | `ucsf_fixture_repaired/`, `followup_output_parity_v2.json` |
| Comparator/auditor unit tests | 18 host-audit/HDF5 tests; five SAM-multiset tests | Test commands below |

The new downstream aggregate report passes **51/51 checks**: six repeated-run
H5AD comparisons, four text/JSON comparisons, six cross-version H5AD
comparisons and 35 UCSF MEX/axis comparisons. All HDF5 datasets, attributes,
node sets, dataset shapes and dtypes are checked; categorical codes and axes
are not reordered or dropped. HDF5 physical layout is not compared.

Previously verified current Flex (121/121 outputs), OCM/Velocyto, PF dynamic
stress, separate CUDA CellBender/layer smoke, and partial-build gates remain
documented in the handoff. They were not all repeated in this follow-up.

## Test-harness corrections

### Unsorted SAM order

STAR's default `outSAMorder=Paired` does not guarantee global input order.
The SLAM harness now records a differing order as `ORDER_ONLY` only after
proving the complete SAM record multiset identical, preserving duplicates and
every field. Changed or missing alignments still fail. The independent sorted
SAM, junction, diagnostic and quantification checks remain required. Input
FASTQ-vs-CBQ payload comparison remains ordered and byte-exact.

The test refuses an existing output directory. Five unit tests cover exact
identity, permutation, duplicate loss, changed fields and empty inputs.

### Reused UCSF fixture configuration

The initial new smoke's reused FASTQs were linked under a new fixture path,
but its copied `pf_multi_config.csv` still named the old paths. STAR reported
`pf-multi GEX input filtering matched 0/20 FASTQ files; leaving input unchanged`.
That run completed but incorrectly included guide inputs in GEX counting,
producing 5,748 rather than 3,943 cells. It is **not accepted parity evidence**.

The smoke now lets the workflow generate its configuration from the actual
fixture paths and fails on the zero-match warning before downstream work.
The repaired run returns 3,943 cells and exact matrix/axis parity. No production
datasets were altered, and no STAR counting algorithm was changed here.

Preserved non-acceptance records:

- `followup/ucsf/`: stale fixture configuration, not a production regression.
- `followup/tier_molecule.log`: missing default `/usr/local/bin` helper.
- `followup/{paired_hash,ambient_fdr}.log`: tests were initially pointed at
  the working tree rather than the built libraries in the fresh checkout.
  Corrected companion checks use the fresh build and all pass.
- Early downstream comparison reports that flagged provenance paths remain;
  the final `followup_output_parity_v2.json` includes explicit root normalization.

## Original G-S1 issue: TranscriptVB

The expanded read-only capture audit paired all 100 original candidate STAR
invocations without ambiguity. **95 pairs pass**. Eight additional baseline
invocations are the previously held SLAM arms, now covered by the completed
follow-up. The five differing paired captures are:

- `00065`: no-4sU SE TranscriptVB quantification.
- `00066`: no-4sU PE TranscriptVB quantification.
- `00067`: one gene length differs by 0.001 in the treated SE output.
- Baseline `00103/00104` versus candidate `00095/00096`: TranscriptVB
  evidence sidecars.

The no-4sU R1 and R2 decoded FASTQs are identical between original arms
(hashes recorded in `transcriptvb_diagnostic/inputs.json`). Controlled
diagnostic runs use those same inputs and settings, changing only output
paths, executable, and the stated thread count:

| Comparison | Transcript NumReads differences | Maximum absolute difference |
|---|---:|---:|
| Original vs repeated **1.9.5.a**, PE, 16 threads | 1,863 transcripts | 34.017 |
| 1.9.5.a vs 1.10.0, PE, 16 threads | 1,875 transcripts | 44.376 |
| 1.9.5.a vs 1.10.0, SE, 16 threads, new controlled runs | 0 | 0 |
| 1.9.5.a vs 1.10.0, SE and PE, one thread | 0 | 0 |

In the one-thread controls, the entire parsed transcript, gene and tximport
tables match exactly, including TPM and effective lengths. The PE 16-thread
cross-version comparison changes gene NumReads in 38 rows, maximum 1.0;
the repeated baseline changes the same number of gene rows, maximum 1.0.
These are not merely output-order differences. The earlier SE discrepancy
was not reproduced by the new controlled SE comparison.

The unchanged `TranscriptQuantEC.cpp` seeds the learning RNG by thread ID,
uses stochastic alignment acceptance for online FLD/error-model updates,
and merges thread-local learned models. This is consistent with pre-existing
scheduling-dependent output; it does **not** justify silently accepting a
strict byte-identity gate. No TranscriptVB production algorithm was changed.

Both sidecars have valid checksums, identical equivalence-class multisets
and GC counts. One has an FLD-bin difference of `8.881784197001252e-16`;
the other has byte-identical evidence payload. Their provenance embeds the
different source revisions, temporary reference paths and path/mtime-based
FASTA fingerprints. These fields and the rounding difference are reported
explicitly, not relabeled as byte-identical.

Full reports: `capture_audit_v1.json`,
`transcriptvb_diagnostic/{commands.json,report.json}` and
`sidecar_content_diagnostics.json`.

## Ordered single-thread acceptance

The owner chose single-thread validation, not a numerical tolerance or a
sorted-record comparison. The 100K inputs were retained without further
downsampling. Reused the completed S50 no-4sU one-thread controls and executed
S43 treated SE/PE on both frozen binaries. All **12 quantification tables**
(transcript, gene and tximport for four cases) match byte-for-byte, including
row order. No production TranscriptVB code or defaults changed.

The synthetic 4,000-pair scatter/gather smoke ran on both versions with
`THREADS=1` (mapping, reference generation and finalize). The two sidecar
evidence payloads match byte-for-byte in order, including ECs, weights, GC and
FLD bins. Checksums and fixed metadata validate. Only the explicitly verified
source revision, reference directory and path/mtime-based FASTA fingerprint
differ in their provenance headers; reference/FASTA/read files match exactly.
In-process and gathered quantification outputs also match across versions.
The audit passes **27/27 checks** in
`single_thread_acceptance/ordered_parity.json`. The smoke's existing within-arm
scatter-versus-in-process tolerances are separate from this exact cross-version
audit; they were not used to accept cross-version differences.

The legacy external BINSEQ probe's auto-thread decoder also varied record
order. With the same pinned CBQ input and `bqtools decode --threads 1`, both
versions now produce byte-identical R1, R2 and probe TSV outputs, preserving
all headers and order (`binseq_single_thread/PASS`). No runtime reader code
was changed.

The final kept-output audit accounts for all **218 initially different or
missing paths**, plus the original 2,042 identical outputs. There are no
unexplained paths. Details are in `kept_output_disposition_v2.json`: logs,
fixture Git/font caches, verified provenance/timing fields, approved keyed PF
matrices, unordered transfer inventory, intentionally moved Chromap adapter
outputs, and explicitly superseded OCM/downstream attempts. The old reports
and failed outputs remain intact.

The seeded downstream comparison additionally passes all six summary, QC,
doublet and plot checks. The diagnostic `inspect_anndata.py` samples values
without a fixed seed when describing inferred semantic types: this caused one
different description of the identical ambiguous layer. Re-rendering both
summaries with NumPy seed 1 produces identical descriptions. This is a
diagnostic-only control, not a change to the H5AD data or production inspector.
PNG bytes match, and HTML matches after replacing its generated Plotly ID.

Acceptance record: `$D/G_S1_ACCEPTED.json`, generated only after these checks
and the frozen binary SHA256 checks succeeded. The historical multi-thread
differences remain documented; deterministic multi-thread TranscriptVB is not
claimed.

## Three-repeat timing results

All eighteen planned workload attempts ran, with matching scientific outputs
in all eighteen. The 100K scRNA regression includes vanilla, modern and
modern-BAM modes; 10M scRNA uses Gene/GeneFull without BAM; Flex uses the
established eight-lane H0/H1X2 fixture. Threads remain 8, 16 and 32 respectively;
these are separate from the single-thread TranscriptVB correctness controls.

**Accepted local regression-gate medians, using all three observations per
version.** These are not noise-free headline performance measurements.

| Workload | 1.9.5.a wall (s) | 1.10.0 wall (s) | Wall change | Peak RSS change | Host-flagged runs, old/new (of 3 each) | Gate |
|---|---:|---:|---:|---:|---|---|
| scRNA 100K, three modes | 55.65 | 55.67 | +0.036% | +0.0023% | 0 / 0 | PASS |
| scRNA 10M | 39.63 | 39.69 | +0.151% | -0.00005% | 0 / 1 | PASS |
| Flex 8 x 100K | 15.36 | 15.38 | +0.130% | +0.0348% | 2 / 3 | PASS |

All observed changes are below the +3% wall / +2% RSS limits. **Six attempts
retain their original host-monitor flags**. The first candidate 10M run started
at load 4.33 (limit 4). The first baseline Flex run started at load 6.38 with
about 496 MB of unattributed array I/O. Four later Flex runs had low start
load but 238-373 KB of unattributed I/O; this exceeds the monitor's 5% fraction
because these short, cached runs issue very little disk I/O themselves. These
flags remain in the original logs and both reports. No observation was dropped
from the medians.

The owner's subsequent instruction was: "Don't worry about host activity as
long as it falls within tolerances. If it doesn't then you can flag that as
a reason". Accordingly, host activity is diagnostic, not an independent veto
when median wall time is within +3% and median peak RSS within +2%. An
out-of-tolerance result still fails; host activity may be investigated as a
possible contributor, not assumed to excuse a regression. Completed workloads
and output parity are still required. This changes the acceptance policy, not
the host monitor's verdicts or the measurements.

Execution/disposition records:

- `$D/gs3_acceptance/`: initial five attempts. The baseline Flex wrapper exited
  2 on absolute build-header path checks in the archived copy; STAR itself
  completed and all 121 output files matched. The dependency paths name the
  original build directory, which still exists. Subsequent baseline executions
  use that original binary path, with the same verified SHA256. Reassessment
  independently verifies `out/A/rc=0`, the finished 800,000-read final log,
  matching output signatures, and exactly those two build-path guard failures.
  This observation is accepted as a completed STAR workload with a documented
  harness failure, not reclassified as a successful wrapper execution.
- `$D/gs3_remaining/`: the thirteen previously unexecuted attempts, each with
  zero exit status and `WORKLOAD_COMPLETE`. Untimed preflights drain preceding
  writes and wait for low load. Sampling changed from 5 s to 0.5 s to observe
  short-lived children; host-verdict thresholds did not change. The later idle
  preflight permits up to 1 MiB background writes per five seconds instead of
  demanding absolute zero; it does not decide timing acceptance.
- The resumed batch shell still exited 2 **after** its last completed workload:
  the agent edited the live driver, causing Bash to resume reading at a stale
  file offset. This execution mistake is documented in
  `gs3_remaining/BATCH_EXIT_DISPOSITION.md`. The per-run evidence is audited
  independently; no successful batch marker was fabricated. Future executions
  must use `tools/run_gs3_frozen.sh`, which snapshots the driver and idle helper.
  Its existing-output guard was tested without launching a data workload.
- `$D/gs3_round3/`: cancelled replacement driver; no timer or STAR workload
  launched. No data execution was duplicated. Original driver revisions and
  the first idle helper are retained under `$D/tools/`.
- `$A/gs3_report.json`: 18/18 matching output signatures, all attempt timings,
  peak RSS, host verdicts, the medians above and the original strict-policy
  `accepted: false` result. Preserved unchanged.
- `$A/gs3_report_v2.json`: read-only reassessment under the revised policy,
  **G-S3 accepted**. Binary hashes, raw timing/completion evidence and all
  eighteen output signatures were rechecked; no workload was launched.
  `tools/test_reassess_gs3.py` passes ten synthetic checks, including rejection
  of failed/missing workloads, duplicate repetitions, changed outputs and
  out-of-tolerance measurements despite host activity.

The prior request for additional clean-host measurements is superseded by this
policy clarification. No further identical repetitions were launched or are
needed solely for these host flags. Acceptance is recorded in
`$D/G_S3_ACCEPTED.json`; current combined status is
`$A/VALIDATION_STATUS_v3.json`. Earlier status files remain historical records.

## Release disposition

- **G-S1 accepted** under the owner-selected exact, ordered single-thread
  TranscriptVB contract and the explicit kept-output dispositions above.
- **G-S2 passed** on the new clean build.
- **G-S3 accepted** under the owner-clarified host-activity policy. All eighteen
  outputs match; all three workloads meet the unchanged numerical limits.
  Host flags and separately audited wrapper errors remain preserved.
- Recipe default-branch integration remains outstanding. Multiomics host
  integration and recipe/snapshot migration remain stable-release dependencies
  after the local RC. No release was pushed or tagged.

## Commands

These are records, not instructions to overwrite existing output directories.
One-off drivers and detailed commands stay in the artifact root:

```bash
flock /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock \
  bash "$D/tools/build_acceptance_candidate.sh"
bash "$D/tools/run_acceptance_followup.sh"
bash "$D/tools/run_ucsf_fixture_repair_validation.sh"
bash "$D/tools/run_host_acceptance.sh"
bash "$D/tools/run_companion_acceptance.sh"
bash "$D/tools/run_seeded_baseline_downstream.sh"
python3 "$D/tools/check_transcriptvb_repeatability.py"
python3 "$D/tools/complete_capture_audit.py" --report "$A/capture_audit_v1.json"
python3 "$D/tools/validate_followup_outputs.py" --report "$A/followup_output_parity_v2.json"
python3 "$D/tools/run_single_thread_acceptance.py"
python3 "$D/tools/verify_single_thread_acceptance.py"
bash "$D/tools/check_binseq_ordered.sh"
python3 "$D/tools/seed_inspection_summaries.py"
python3 "$D/tools/close_kept_output_audit.py"
python3 "$D/tools/accept_gs1.py"
STAR_CANDIDATE_BIN="$A/src110/core/legacy/source/STAR" bash "$D/tools/run_gs3.sh"
# Resume of only previously unexecuted attempts, with original build directory:
OLD=/tmp/claude-1000/-mnt-pikachu-chromap-suite-paper/53e97281-48e6-45f8-90f8-a4913b843536/scratchpad/gate
STAR_CANDIDATE_BIN="$A/src110/core/legacy/source/STAR" \
  STAR_BASELINE_BIN="$OLD/src195/core/legacy/source/STAR.real" \
  GS3_PRIOR_OUT="$D/gs3_acceptance" GS3_OUT="$D/gs3_remaining" \
  bash "$D/tools/run_gs3.sh"
python3 "$D/tools/report_gs3.py" \
  --roots "$D/gs3_acceptance" "$D/gs3_remaining" --report "$A/gs3_report.json"
# Saved-evidence reassessment only; does not execute STAR or overwrite reports:
python3 -m unittest discover -s "$D/tools" -p test_reassess_gs3.py -v
python3 "$D/tools/reassess_gs3.py"
# In the STAR candidate worktree:
python3 -m unittest discover -s tests/host_api -p 'test_*.py' -v
python3 -m unittest discover -s tests/slam -p test_compare_sam_records.py -v
```

The drivers take the host lock per workload. Original nonzero driver exits
are retained; replacement evidence is explicitly listed above rather than
rewriting historical status tables. The initial timing driver is preserved
as `tools/run_gs3.initial.sh`, the first resume as `run_gs3.resume_v1.sh`.
Do not edit an executing driver or helper; prepare a new immutable snapshot.
