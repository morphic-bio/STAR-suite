# Handoff: RNA mate reader threads (STAR Suite 1.11.0)

Updated: 2026-09-29, second implementation session. Runbook:
`docs/runbooks/RUNBOOK_RNA_MATE_THREADS_20260929.md`. This handoff supersedes
the runbook where they differ. The runbook edits that were pending are done
(status, decisions, sections 3.3 steps 2, 7 and 8, 3.5, 3.6, 3.7, 4.1).

## Where the work stopped (read first)

**Stopped at 00:59 UTC, 30 Sep, in M3: blocked on the benchmark lock for 30
minutes** (the coordinator's limit). My queue's next run (`b0 f6d`) waited on
`/mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock` from 00:29:54. The
holder was the single-run agent's `validation/single_run/run_gates.py` (held
31 minutes by then), with its `gs4.py` and `gs5.py` queued behind it. My queue
chains were stopped; no process of this work is running or waiting, and the
empty `runs/f6d/b0_r1` was removed. An earlier wait (23:45-00:05, behind its
`gs1.py`) cleared after 20 minutes.

**Resume:** `$W/tools/run_queue.sh $W/queues/m3_resume.txt` (remaining B0
variants, the B1 F6c rerun, every C run with traces), then
`$W/tools/m3_report.sh f1a f1a1 f1b f1c f1d f1e f1f f2 f4 f5a f5b f6a f6a1
f6b f6c f6d f6e f6f f7vanilla f7modern f7modernbam`, then
`tools/run_flex_smoke.sh b0` and `c` (G-R3), then queue `m4` (F3).

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
| C | `$W/src_c` | `4329a6ab3b4147688875a09ea75d73a4a3cd5b8d` | `$W/src_c/core/legacy/source/STAR` (also `libstar_suite.a`) |

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
| M1 | Committed (`130af40`). **B1 = B0 on F1 and F2: PASS** (F1a fully identical; F2 within the Multiomics lists, Solo.out exact). Also PASS on every other variant with both runs: f1a1, f1b-f1f, f4, f5a, f5b, f6a, f6a1, f6b, f7vanilla, f7modern, f7modernbam. B1 traces collected for every variant except f6c: its first B1 run used a misspelled option (`--twoPassMode`; STAR's is `--twopassMode`), failed with exit 102 and is kept as `runs/f6c/failed_b1_twoPassMode_typo`; the argv is fixed and the rerun is in `m3_resume`. |
| M2 | Committed (`75608c3`), blank-line rule `4329a6a`; G-R0 passes (272 identity runs, behaviour, blank-line part). |
| M3 | **Partly run.** B0 done for f1a1, f1b-f1f, f5b, f6a1, f6b, f6c (two-pass ran both passes); f6d-f6f and all C runs not yet (stopped on the lock). No G-R1/G-R2 result yet. |
| M4-M5 | C built at `4329a6a` (`libstar_suite.a` too); queue `m4`, G-R3 wrapper `tools/run_flex_smoke.sh`, report script `tools/m3_report.sh` ready. Not run. |
| M6 | Not started; stop before G-S1 and the release. |

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

## Next steps

1. Queue `m3_resume` (see "Where the work stopped"), then
   `tools/m3_report.sh <variants>`.
2. G-R3: `tools/run_flex_smoke.sh b0` and `c`; check C's `Log.out` for
   "Fastx mate readers: not active (Flex ...)".
3. M4: queue `m4` (F3), G-R1/G-R2 on F3, informal measurements.
4. M5: G-S2 in `$W/src_c` (under the lock), then the Multiomics development
   build against C and G-M1 under the D6 conditions.
5. Stop before G-S1 and the release.

## For the author

- **Blank-line rule (author decision, 29 Sep): implemented** (runbook 3.7,
  draft `docs/RELEASE_NOTES_v1.11.0.md`). Two implementation choices within
  the decision need confirming:
  - a run of consecutive blank lines counts as one gap: every line is
    counted, and the header check applies to the first non-blank line;
  - blank lines at the end of a file (then end of input or the next input
    file) are skipped and counted, not fatal.
- **F8 second script not run.** `tests/run_flex_tiny_public_smoke.sh` clones
  a third-party repository (`minoda-lab/universc`) to get test data derived
  from the vendor's tiny test set. Under the clean-room rule it was not run
  without asking. G-R3 uses `tests/run_flex_half_probe_100k_smoke.sh` (local
  fixtures, byte comparison with the preserved v1.9.5a reference outputs).
- Pre-existing, unchanged: with multi-file FASTA through lane markers, the
  FASTA loop reads the next file's `FILE` marker line as sequence, in
  v1.10.0 and here alike.

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
  the command, no workaround), or near 90% of the usage limit.
