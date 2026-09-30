# Handoff: RNA mate reader threads (STAR Suite 1.11.0)

Updated: 2026-09-29, second implementation session. Runbook:
`docs/runbooks/RUNBOOK_RNA_MATE_THREADS_20260929.md`. This handoff supersedes
the runbook where they differ. The runbook edits that were pending are done
(status, decisions, sections 3.3 steps 2, 7 and 8, 3.5, 3.6, 3.7, 4.1).

## Current state (read first)

**M0-M5 are done; stopped before G-S1 as instructed.** No output difference
from v1.10.0 outside the B0 variance classes was found in any gate.

- C is `4bc3b8f` (blank-line rule per the author's decisions of 29 and 30
  Sep), built from clean in `$W/src_c` with `libstar_suite.a`; G-R0 passes
  there (272 identity runs, behaviour, blank-line part: 0 failures).
- G-R1/G-R2 (M3, 21 variants; M4, F3 20M), G-R3 (Flex), G-S2 and G-M1 (5
  fixtures) pass; details below.
- D6 was carried out and cleaned up: the Multiomics development build used a
  local, unpushed tag `v1.11.0-dev.rna-mate-threads` in a scratch STAR Suite
  clone and a scratch manifest via `--manifest`; the scratch area
  (`$W/scratch/d6`, clones, tag and build) is deleted. The tag never existed
  outside that clone (checked: worktree, `/mnt/pikachu/STAR-suite`, GitHub).
  Build record kept in `$W/m5/d6_record/` (receipt, inputs, sources lock,
  scratch manifest, build log, binary sha256).
- Lock policy: queue behind the lock with no cutoff; stop only if a single
  wait passes about 3 hours. `run_star.sh` supports batch holding
  (`flock <lock> env RUN_STAR_LOCK_HELD=1 run_queue.sh <queue>`).

Earlier stop (resolved): the first F2 dry run was denied ("[Safety Bypass
Flag]") because of `--skip-active-check`; the coordinator pointed out that
the lane runner's active-run check runs only without `--dry-run`
(`run_dogmaplex_lane.py:335-343`), and the dry run without the flag was
allowed. Real runs never use `--skip-active-check`.

## State

- The author approved the runbook. Decisions:
  - D1 `auto`;
  - D2 truncate per file with a WARNING;
  - D3 changed: reader threads take MAP-domain decode permits, only while
    parsing bytes already in memory;
  - D4 keep the helper and FIFO;
  - D5 decide after M4;
  - D6 approved with conditions (below);
  - D7 separate change after 1.11.0;
  - D8 resolved.
- Branch `design/rna-mate-threads-20260929`, on top of `origin/master`
  `0f9701a`. There are no `core/` changes between `9090fb4` and `0f9701a`.
- Commits (on the branch only; not pushed, merged or tagged):
  - `130af40` M1: optional per-mate CRC32 columns in the input chunk trace
    (`STAR_INPUT_CHUNK_TRACE_DIGEST=1`).
  - `75608c3` M2: the mate reader threads, the option, the gate, open, close
    and mapping changes, the harness, the smoke script, the manifest row,
    `.gitignore` and the regenerated `parametersDefault.xxd`.
- The runbook and this handoff are uncommitted. Earlier release handoffs and
  runbooks were committed, so commit these later with the work.

## Verified in this session

- **G-R0 passes.** `tests/run_fastx_mate_reader_harness_smoke.sh`: at
  `4329a6a`, 272 identity runs (17 cases without blank lines), 0 failures, no
  unexpected WARNING; behaviour 0 failures; blank-line part 0 failures.
  (At `75608c3`: 304 identity runs including the two old blank-line cases.)
  - 19 cases × 4 chunk sizes × readID Number × readMapNumber, with a
    one-permit pool on alternate runs.
  - Repeated 30 more times: 0 failures.
  - ThreadSanitizer build of the harness: no reports.
- The harness now also compares the "end of input stream" Log.out lines.
  Cases added: `mate2_no_final_newline`, `line_limits_fasta` (650-byte FASTA
  lines), `line_limits_fastq` (650-base reads, header lines near 50,000 bytes)
  and `fasta_line_over_limit`.
- Candidate STAR builds with no new warnings. `libstar_suite.a` contains
  `FastxMateReaders.o`, through `STAR_HOST_OBJECTS` (`OBJECTS` minus
  `STARmain.o`), not `OBJECTS_NO_MAIN`.
- `parametersDefault.xxd` equals `xxd -i parametersDefault`.
  `tests/test_parameters_default_generation.py` passes.
- `--readFilesMateThreads bogus` is rejected with the intended message.
- `tests/run_fastx_contract_star_smoke.sh` passes with the candidate (default
  `auto`; readers active in 5 runs; output in
  `$W/smoke/fastx_contract`). This is a smoke only; it is not an identity
  gate.

## Code changes made in this session (in `75608c3`)

- **End-of-input log line.** The new filler wrote
  "end of input stream, nextChar=-1" where v1.10.0 wrote nothing: for a file
  with no final newline, and for every FASTA input, where the last peek
  reaches EOF. v1.10.0's loop tests mates 1 and 2 for `good()` before each
  pair and ends silently in those cases. Each end batch now records whether
  the stream had already failed (`endSilent`), and the filler follows the
  same rule. The harness catches the old behaviour (checked with the rule
  disabled). Log.out only; chunk text was never affected.
- A reader that has been asked to stop no longer waits for a new permit.
- The Log.out summary says "N reads" rather than "N read pairs".
- `.gitignore`: `core/legacy/source/fastx_mate_reader_harness`.

## Builds (M0/M1/M3 inputs)

`W=/mnt/pikachu/star_rna_mate_threads_20260929`. All three are local clones
checked out detached, built on the host with `nice -n 10`, `-j16`,
`make STAR STAR_SUITE_COMMIT_SHA=<full sha>`. Each binary reports "1.10.0"
and embeds its commit (C was rebuilt from clean after the blank-line commit;
an incremental rebuild keeps the old revision string).

| Label | Source | Commit | Binary |
|---|---|---|---|
| B0 | `$W/src_b0` | `0f9701a4d5145f211dc51797b0b21c3d551912c0` | `$W/src_b0/core/legacy/source/STAR` |
| B1 | `$W/src_b1` | `130af40050b954cb3acc905946918c832603b7dc` | `$W/src_b1/core/legacy/source/STAR` |
| C | `$W/src_c` | `4bc3b8f2db8cf889c226117ce19b8b0531d93dd7` | `$W/src_c/core/legacy/source/STAR` (also `libstar_suite.a`) |

## Gate tooling (in `$W/tools`, untracked)

- `stage_fixtures.sh`: plain, BGZF, single-file and manifest variants of F1;
  plain F2 GEX (M0 diagnostic); plain F6.
- `make_argv.py`: writes `$W/argv/<id>.args` (STAR-only argv, one argument
  per line): f1a (zcat, lanes), f1a1 (1 thread), f1b (internal gzip), f1c
  (BGZF), f1d (plain lanes), f1e (plain single file), f1f (manifest), f2,
  f2plain, f3, f4, f5a (zcat), f5b (no command), f6a (sorted BAM +
  GeneCounts), f6a1 (1 thread), f6b (unsorted), f6c (two-pass), f6d
  (BySJout), f6e (SE), f6f (plain), f7vanilla, f7modern, f7modernbam. F2 and
  F3 come from `run_dogmaplex_lane.py --dry-run` (in `$W/runs/dryrun*`).
- `run_star.sh BIN ARGS OUT [trace]`: holds the benchmark lock per run; if
  another STAR/multiomics process runs, releases it and retries; exit 75
  after 30 minutes blocked that way.
- `run_queue.sh QUEUE`: queues in `$W/queues/` (m0, m1, m3); outputs in
  `$W/runs/<id>/<label>_r<N>`.
- `compare_run.py` (runbook 4.2), `compare_traces.py` (G-R1),
  `compare_pair.sh`; reports in `$W/compare/`.

## Milestones

| M | State |
|---|---|
| M0 | **Done.** Fixtures staged; argv in `$W/argv`; B0 pairs run (queue `m0`, all exit 0). Variance classes and diagnostic below. |
| M1 | Committed (`130af40`). **B1 = B0 on every variant (F1-F7, F3), including F1 and F2 as required.** B1 traces collected for all. The first B1 f6c run used a misspelled option (`--twoPassMode`; STAR's is `--twopassMode`), exited 102 and is kept as `runs/f6c/failed_b1_twoPassMode_typo`; rerun with the fixed argv. |
| M2 | Committed (`75608c3`), blank-line rule `4329a6a`; G-R0 passes (272 identity runs, behaviour, blank-line part). |
| M3 | **Done: G-R1, G-R2, G-R3 pass.** 21 variants (below): C trace = B1 trace row for row (both CRC32s), C outputs = B0 outputs under runbook 4.2 and the variance classes, B1 = B0; C's Log.out says the readers were active (twice for two-pass). G-R3: B0 and C both reproduce all 121 non-log outputs of the preserved v1.9.5a Flex reference; C logs "Fastx mate readers: not active (Flex runs read their own lanes)". |
| M4 | **Done: G-R1 and G-R2 pass on F3** (391 chunks, 20M reads; three C runs vs B0). Informal measurements below. D5 (M4b) is the author's call. |
| M5 | **Done: G-S2 and G-M1 pass.** G-S2 on C (`$W/m5/gs2_c`): PASS, readers active in all 31 STAR/host runs, including the host with a dummy External domain. G-M1 (Multiomics development build against C, `$W/g_m1`): 5 of 5 fixtures pass under the committed lists; details below. |
| M6 | Not started; G-S1 not run (stop before G-S1 and the release). Remaining M6 work: README option list if it lists input options, final release notes, G-R5 (G-S1 plus Tier A) with the new default. |

### M0 results (B0 repeatability pairs, `$W/compare/`)

- F1a (8 threads), F6a (8 threads), F7 vanilla/modern/modern_bam (8 threads):
  every output identical between the two B0 runs (runbook 4.2 rules,
  including `Log.final.out` count lines).
- F2 (16 threads): differences only in `cr_assign/` feature-arm files (25
  order-only, 13 content, including three `assignBarcodes.api_run.txt`
  timing lines) and `outs/feature_analysis/celltag_dp01/*/matrix.mtx`; all
  within the committed Multiomics lists for `dogmaplex_lane1_2m_sidecar`
  (comparator verdict PASS). `Solo.out`, `SJ.out.tab` and `Log.final.out`
  counts identical.
- F4 (8 threads): identical except `assignBarcodes.api_run.txt` processing
  time (normalized by the Multiomics comparator).
- F5a (8 threads): guide-arm `cr_assign/` order-only files and
  `feature_sequences.txt` match positions, within the CAT-ATAC lists (paths
  without the recipe's `star_run/` prefix).
- Variance classes for C: F1, F6, F7 must be byte identical; F2, F3, F4, F5
  may differ only in the Multiomics-listed feature-arm files, as judged by
  `tools/verdict.sh` (strict comparison plus the Multiomics comparator; a
  `Log.final.out` difference is never accepted).

### M0 diagnostic (informal; shared host; no claims)

F2 with B0, 16 threads, 52 chunks: "Avg chunk read time" 62.08 and 63.72 ms
with `.gz` input (two runs) against 41.27 ms with the same GEX mates as
plain FASTQ. Plain input is not as slow as `.gz`, and the plain-input chunk
read (no decompression wait) is still about two thirds of the `.gz` value,
so the parse is the larger part. The runbook's M0 stop condition is not met.

### M3 results (G-R1, G-R2; `tools/m3_report.sh`, reports in `$W/compare/`)

All PASS: f1a (zcat, lanes), f1a1 (1 thread), f1b (internal gzip), f1c
(BGZF pipe group, "BGZF raw input: active"), f1d (plain lanes, `cat`), f1e
(plain single file, direct), f1f (manifest), f2, f4, f5a (zcat,
`--readMapNumber` on full files), f5b (no command), f6a (sorted BAM +
GeneCounts), f6a1 (1 thread), f6b (unsorted BAM), f6c (two-pass; 36 trace
rows over both passes), f6d (BySJout), f6e (single-end), f6f (plain),
f7vanilla, f7modern, f7modernbam.

- 1-thread runs (f1a1, f6a1) and F1, F6, F7 at 8 threads: every output byte
  identical to B0, including `Log.final.out` count lines.
- F2, F4, F5: differences only in feature-arm files the Multiomics lists
  cover; `Solo.out`, `SJ.out.tab` and `Log.final.out` counts identical.
- f6b unsorted BAM: B0, B1 and C each wrote the same 1,099,014 records in a
  different order (thread scheduling); compared as a multiset per 4.2.

### M4 results (F3: DOGMA-plex lane 1, 20M, STAR only, 32 threads)

- G-R1: trace equal (391 rows, 20,000,000 reads). G-R2: PASS for C runs 1-3
  against B0 run 1, and C run 2 against B0 run 2.
- Variance class found at 20M: `cr_assign/CellTag/celltag_dp01/ct/
  feature_sequences.txt` differs between the two B0 runs (and B1) in the
  "Match Position" field of mismatched sequences, as the ADT and HTO files do
  in the committed 2M list. Recorded from the B0 pair in
  `$W/compare/lists/content_varying_f3_extra.txt`
  (method `feature-sequences`) and used for F3 only. Not a C difference.
- Informal measurements (shared host, order B0, C, C, B0; no claims):

  | run | wall | avg chunk read ms | mutex wait thread-s | map chunk thread-s | MAP/FEATURE/idle permit occupancy |
  |---|---|---|---|---|---|
  | B0 r1 | 3:12 | 93.57 | 699.1 | 445.2 | 11.8 / 0.51 / 19.7 |
  | C r2 | 3:16 | 52.23 | 124.4 | 539.7 | 24.9 / 0.01 / 7.1 |
  | C r3 | 3:14 | 53.48 | 127.0 | 552.5 | (similar) |
  | B0 r2 | 3:20 | 93.57 | 684.6 | 459.8 | (similar) |

  Reader summary (C r2): mate 1 (R2) parse 7.4 s, input wait 13.1 s, filler
  waited 17.6 s for mate 1; mate 2 blocked 13.7 s on a full queue. The
  ceiling has moved to mate 1's decompression (R2 on one core), as the
  runbook predicted; wall time is unchanged at this scale. MAP occupancy
  rose toward `runThreadN` (the runbook's expected effect). FEATURE permit
  acquisitions fell from about 1.57M to about 0.1k with identical feature-arm
  phase times and outputs; the cause was not established (observation only).
- D5 (M4b fast parser): the mate-1 parse (7.4 s) is not the ceiling here;
  the decompressor is. Recommendation for the author: no M4b for 1.11.0.

### M5 results

- **G-S2** (`tests/host_api/run_host_api_tests.sh` in `$W/src_c`, under the
  lock): PASS. The `.gz` fixture made the readers active in all 31 runs.
- **G-M1** (`$W/tools/run_g_m1_c.py`: the committed perf-lane driver with
  only ROOT, OUT, BIN, HERE and `MULTIOMICS_MANIFEST` changed; outputs
  `$W/g_m1`): exit codes match the reference for all five fixtures.
  - PBMC 3k, HIV DOGMA four-arm, CAT-ATAC: comparator PASS as run.
  - DOGMA-plex lane 1 (2M) and five-arm: as run, 10 files each differed
    only in absolute paths to the reference's Multiomics checkout
    (`/mnt/pikachu/multiomics-suite-single-binary-20260928`; each size gap
    exactly the 45 bytes between that path and `<REPO0>`). With
    `--reference-repo-root` (the comparator option added for this in 0.10.0,
    per its release notes) both PASS, no list changed
    (`comparison_with_reference_repo_root.{log,json}`).
  - `.gz` representation audit: one entry, `atac/peak_mex/matrix.mtx.gz` in
    lane 1 (compressed bytes differ, decompressed identical). It is written
    by the ATAC peak layer, not STAR, and the perf-lane G-M1 of another
    Multiomics build flagged the same file; not attributable to this change.
- The first G-M1 attempt failed to start four fixtures because the copied
  driver pointed `MULTIOMICS_MANIFEST` at the clone's committed manifest
  (STAR v1.10.0) while the build used the scratch manifest; kept as
  `$W/g_m1_failed_manifest_env`, fixed, rerun.

## Next steps

1. Author review of the results and of the items in "For the author" below.
2. M6 (after approval): README option list if it lists input options, final
   1.11.0 release notes, G-R5 (G-S1 plus Tier A) with the new default in
   force, final handoff. No tag, push, merge or release.

## For the author

- **D5 (M4b fast parser):** at F3 (20M, 32 threads) the mate-1 reader parses
  in 7.4 s and waits 13.1 s on its input; the ceiling is mate 1's
  decompression, not the parse. Recommendation: no M4b for 1.11.0.
- **Observation (informal):** with C, MAP permit occupancy rises toward
  `runThreadN` (11.8 -> 24.9 of 32 on F3), as the runbook expected; FEATURE
  permit acquisitions fall from about 1.57M to about 0.1k with identical
  feature-arm phase times and outputs. Cause not established; relevant to
  the hosted lane, where the host's floors govern the split.
- **Log.out "not active" checks for SLAM, TranscriptVB and CBQ** (runbook
  4.2) were not run here; those modes are G-S1 rows (M6).

## Author answers to the M5 report (30 Sep)

1. **FEATURE permit drop: find the cause first** (about 1.57M acquisitions on
   B0 against about 0.1k on C at F3). Use existing logs and permit traces
   plus at most one or two instrumented runs under the lock. Establish what
   changed; confirm that the feature arm still takes permits as designed, so
   every worker goes through the allocator. If the cause needs a change to
   the allocator design, Flex, CBQ or the feature readers, or a permit path
   is being bypassed: stop and report without changing anything. Put the
   explanation, with evidence, in this handoff.
2. **D5 (fast parser) is still open.** Until the coordinator sends the
   answer, it is not in 1.11.0; do not start M4b.
3. **M6 is approved** once the permit question is resolved without a stop:
   first the parts that do not depend on D5 (README option list, final
   release notes with the D5 line left open), then G-R5 (G-S1 plus Tier A)
   unless D5 is in by then. Stop and report at the end of M6. No tag, push,
   merge or release (the STAR 1.11.0 release is a separate go).
4. Lock: batching under the lock is fine; wait without a cutoff; stop only
   if a single wait passes about 3 hours.

## Author decisions of 30 Sep (implemented)

- Blank lines (runbook 3.7): a single blank line followed by a header is
  skipped with a WARNING and a count; a second blank line in a row before the
  end of the file is fatal at the second blank line; blank lines at the end of
  a file (then end of input or the next file) are skipped and counted, any
  number; a blank line is never read as a header.
- `tests/run_flex_tiny_public_smoke.sh` is not run. G-R3 uses the half-probe
  100k smoke with local fixtures.
- The multi-file FASTA marker limitation (a `FILE` marker right after a
  FASTA read's sequence is read as sequence, in v1.10.0 and here alike) stays
  as it is, documented.
- Lock waits: queue behind the lock with no cutoff; the single-run agent's
  gate runs use it too. Only if a single wait passes about 3 hours, update
  this handoff and stop to report.

## D6 conditions (approved)

- Run the Multiomics host gate in M5.
- A local, unpushed tag is allowed, only in a scratch clone (not
  `/mnt/pikachu/STAR-suite` or any shared checkout), and only for the
  development build.
- Use a scratch manifest via `--manifest`.
- Delete the scratch clone afterwards, and never push the tag.

## Key finding (unchanged)

- The lane-1 RNA arm is limited by the serial parse under STAR's input lock,
  not by decompression. Each mate is already decompressed in its own process.
- The lock was held for 2,476 of the 2,498 s of RNA mapping.

## Rules in force

- No push, tag (except the D6 scratch tag), merge or release. Stop before
  G-S1.
- Never read 10x code. Follow the exclusions in the maintainers' private
  notes.
- Implement only the reader design in the runbook.
- Keep private names and policy text out of every file and commit.
- Any run that loads a genome index or may exceed about 16 GB RSS holds
  `flock /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock`. Builds use
  `nice -n 10` and at most `-j16`.
- There are no timings and no benchmark claims.
- Stop on any output difference from v1.10.0, on a needed change to Flex, CBQ,
  the feature readers or the allocator design, on a permission denial (record
  the command, no workaround), on a single lock wait of more than about 3
  hours, or near 90% of the usage limit.
