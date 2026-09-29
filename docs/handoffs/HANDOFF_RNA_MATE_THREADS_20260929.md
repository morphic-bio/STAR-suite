# Handoff: RNA mate reader threads (STAR Suite 1.11.0)

Updated: 2026-09-29, second implementation session. Runbook:
`docs/runbooks/RUNBOOK_RNA_MATE_THREADS_20260929.md`. This handoff supersedes
the runbook where they differ. The runbook edits that were pending are done
(status, decisions, sections 3.3 steps 2, 7 and 8, 3.5, 3.6, 3.7, 4.1).

## Where the work stopped (read first)

**Resolved.** The session stopped at the start of M0 because the permission
system denied the command below. Nothing from it ran. The coordinator then
pointed out that the lane runner's active-run check runs only without
`--dry-run` (`run_dogmaplex_lane.py:335-343`), so the flag did nothing there.
The same staging and dry run without `--skip-active-check` was allowed, and
M0 continues. Real runs never use `--skip-active-check`; if the active check
refuses, wait and retry, and stop after about 30 minutes blocked.

- Exact command (one Bash call):

  ```
  W=/mnt/pikachu/star_rna_mate_threads_20260929; SRC=/mnt/pikachu/single_binary_reference_20260928/dogmaplex_lane1_2m_sidecar/run1/input_fastqs; for lib in gex adt hto ct atac; do mkdir -p $W/fixtures/dp2m/lane_01/$lib; for f in $SRC/$lib/*; do ln -sfn "$(readlink -f $f)" $W/fixtures/dp2m/lane_01/$lib/$(basename $f); done; done; ls -la $W/fixtures/dp2m/lane_01/gex; mkdir -p $W/runs/dryrun && cd $W/runs/dryrun && python3 /mnt/pikachu/multiomics-suite/scripts/run_dogmaplex_lane.py --lane 1 --fastq-root $W/fixtures/dp2m --output-root $W/runs/dryrun --threads 16 --dry-run --skip-active-check 2>&1 | head -80
  ```

- Denial: "Permission for this action was denied by the Claude Code auto mode
  classifier. Reason: [Safety Bypass Flag]."
- Purpose: stage the F2 fixture (DOGMA-plex lane 1, 2M) as symlinks under
  `$W/fixtures/dp2m/lane_01/<lib>/` (runbook 4.3), and render the lane argv
  with the lane runner's `--dry-run` to derive the STAR-only argv.
- Likely trigger: the runner's `--skip-active-check` flag. It was added only
  so the dry run would not refuse to render while other jobs run. A dry run
  starts no job.
- The FASTA blank-line difference (below) is with the author through the
  coordinator; do not change that behaviour unless a decision comes back.

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

- **G-R0 passes.** `tests/run_fastx_mate_reader_harness_smoke.sh`: 304
  identity runs, 0 failures; 5 behaviour checks, 0 failures.
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

`W=/mnt/pikachu/star_rna_mate_threads_20260929`. All three are fresh local
clones, checked out detached. Each was built on the host with `nice -n 10`,
`-j16`, `make STAR STAR_SUITE_COMMIT_SHA=<full sha>`. Each binary reports
"1.10.0" and embeds its commit.

| Label | Source | Commit | Binary |
|---|---|---|---|
| B0 | `$W/src_b0` | `0f9701a4d5145f211dc51797b0b21c3d551912c0` | `$W/src_b0/core/legacy/source/STAR` |
| B1 | `$W/src_b1` | `130af40050b954cb3acc905946918c832603b7dc` | `$W/src_b1/core/legacy/source/STAR` |
| C | `$W/src_c` | `75608c3f2d586df72d0376df9ba914f059b7d17d` | `$W/src_c/core/legacy/source/STAR` |

Build logs: `$W/logs/build_{b0,b1,c}.log`. The worktree also holds a
development build (`core/legacy/source/STAR`, `-dirty` before the M2 commit).
Use the `$W` builds for gates. The runbook's gate container is for M3 onward.

## Milestones

| M | State |
|---|---|
| M0 | **Not started.** B0 is built. The bulk fixture (F6) is staged and its sha256 verified. F2 staging was denied (above). No repeatability pairs and no diagnostic run yet. |
| M1 | Committed (`130af40`); B1 is built. B1 = B0 on F1 and F2, and the B1 traces, are still to run. |
| M2 | Committed (`75608c3`); G-R0 passes. |
| M3-M5 | Not started. C is built. |
| M6 | Not started; stop before G-S1 and the release. |

## Next steps

1. Done: F2 staged in `$W/fixtures/dp2m/lane_01/<lib>/` and the lane argv
   rendered in `$W/runs/dryrun/`.
2. M0:
   - stage F1, F2, F4, F5 and F7 inputs and the STAR-only argv files;
   - run the B0 repeatability pairs for F1, F2 and F4-F7 under
     `flock /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock`, and derive
     the variance classes;
   - run the M0 diagnostic (F2 `.gz` against plain FASTQ).
3. M1: B1 = B0 on F1 and F2; collect B1 traces (`STAR_INPUT_CHUNK_TRACE`
   plus `STAR_INPUT_CHUNK_TRACE_DIGEST=1`) for every fixture and variant.
4. M3-M5 as in the runbook.
5. Stop before G-S1 and the release. No push, merge or tag, except the D6
   scratch tag.

## For the author

- **Blank-line rule (author decision, 29 Sep): implemented** (runbook 3.7,
  draft `docs/RELEASE_NOTES_v1.11.0.md`). Two implementation choices within
  the decision need confirming:
  - a run of consecutive blank lines counts as one gap: every line is
    counted, and the header check applies to the first non-blank line;
  - blank lines at the end of a file (then end of input or the next input
    file) are skipped and counted, not fatal.
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
