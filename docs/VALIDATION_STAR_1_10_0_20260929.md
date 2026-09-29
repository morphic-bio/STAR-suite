# STAR Suite 1.10.0 validation follow-up

Date: 2026-09-29 UTC. **Not release acceptance.** The owner authorized the
remaining repeats with "Ok finish the validation gates" on September 28.
The previous authorization hold is resolved; the newly identified output-parity
issue below is not. No tag, push, or master merge was performed.

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

## Outstanding G-S1 issue: TranscriptVB

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

## Release disposition

- **G-S1 is not accepted.** Owner review requested: use the exact single-thread
  TranscriptVB control and document multi-thread variation, or retain the
  strict multi-thread gate and address deterministic model learning first.
  The broader kept-output audit still needs final acceptance/disposition.
- **G-S2 passed** on the new clean build.
- **G-S3 is not run.** Repeat authorization now exists, but the fail-closed
  driver also requires G-S1 acceptance. No acceptance marker was fabricated,
  and no performance claim is made.
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
# In the STAR candidate worktree:
python3 -m unittest discover -s tests/host_api -p 'test_*.py' -v
python3 -m unittest discover -s tests/slam -p test_compare_sam_records.py -v
```

The drivers take the host lock per workload. The follow-up driver's original
exit 1 is retained; replacement evidence is explicitly listed above rather
than rewriting its historical status table.
