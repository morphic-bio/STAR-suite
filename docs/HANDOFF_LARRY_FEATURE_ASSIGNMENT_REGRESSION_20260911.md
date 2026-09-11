# Handoff: LARRY feature assignment regression (MSK ES, STAR Suite 1.9.2)

Date: 2026-09-11. Status: **resolved in STAR Suite 1.9.3** (merge 0143c16). On the full MSK ES P02 test, LARRY assignment took 1,440 s (April 1,443 s) with a LARRY matrix identical to April, and the whole workflow took 1,706 s (April 1,811 s). The record below is kept as the original report.

## Summary

In the paper benchmark arm P02 (MSK 30polyKO ES: GEX + PolyIII guides + LARRY), the whole workflow took
**about 30 minutes in April (v1.4.3)**. On STAR Suite 1.9.2 it was **still in LARRY feature assignment
after more than 3 hours** on 32 cores, having read about 20% of the LARRY FASTQs, and was stopped.

The guide library (PolyIII, 30 features) also slowed: **190 s against 146 s in April (+30%)**. GEX
alignment and cell calling were normal.

The same slowdown shows on A375 (arm P01, which completed): CRISPR guide-library assignment
`stats.processing_time_sec` was **58.4 s on 1.9.2 against 11.9 s** in the March v1.4.x run
(`/storage/A375/paper_bench_20260326_134444`), same inputs and feature reference. A375 feature FASTQs are
BGZF (1.9.2 read them with `--crAssignBgzfMode auto`); the MSK feature FASTQs are plain gzip. Records:
`/storage/paper_bench_v192_20260911/out/P01/cr_assign/CRISPR_Guide_Capture/A375_CRISPR_Guide_Capture_1/assignBarcodes.api_run.txt`
and the same path under the March run.

The owner's reading: **the LARRY slowdown points to a bug in the feature hash load and/or compare, so reads
fall through to the fastHamming search.** Nothing below has confirmed or ruled that out; it is the lead to
test first.

## What was observed

| Stage | 1.9.2 (P02, 2026-09-11) | April v1.4.3 (2026-04-30) |
| --- | --- | --- |
| Whole workflow | stopped at 3 h 10 min | 1,811 s (30.2 min) |
| GEX mapping + Solo | 10:05 -> 10:31 (26 min), 32,898 cells | 33,226 cells |
| PolyIII assignment (30 guides) | 189.7 s | 146.4 s |
| LARRY assignment (245,979 features) | 10:09:25 start; ~20% of each LARRY FASTQ read by 13:12 | 1,442.9 s |

- LARRY progress was measured from the read offsets of the 16 open LARRY FASTQs (`/proc/<pid>/fdinfo`):
  18.8-20.7% of each file after about 3 h. STAR's stdout showed only "Thread N has read 1 million records"
  for eight reader threads over that period.
- CPU: all 32 threads at ~100% for the whole LARRY phase (whole-machine sample ~3,160%); resident memory
  steady at ~41 GB. It was working, not blocked.
- The run was stopped at ~13:14 UTC (STAR exit 143); it is not a valid timing.

## Stack sample (one snapshot, 13:13 UTC, `sudo gdb -p <STAR> -batch -ex "thread apply all bt 8"`)

Frames at the top of the 32 busy threads (counts are frames at depth 0-3 across threads):

| Function | Count |
| --- | --- |
| `pf_anchor_matches_at` | 22 at #2, 19 at #1, 5 at #3, 2 at #0 |
| `pf_find_anchor_position` | 17 at #3, 2 at #2 |
| `pf_anchor_base_matches` | 9 + 8 at #0, 5 at #1 |
| `simple_hash_search` / `feature_lookup_seq` / `simple_search` / `process_feature_sequence_internal` | 2 each |
| futex waits (idle helpers) | 10 |

Full dump: `/mnt/pikachu/star_suite_paper/analysis/paper_benchmark_refresh_20260910/v192/P02/gdb_all_threads_13h13Z.txt`.
Most time is in the anchor search, not in the hash lookup. Whether that means the hash misses and the
search runs for most reads is the question to answer.

## Where the regression can be

- 1.9.2 made no process_features changes: `git diff 21b6fbe 68eb92f -- core/features/process_features`
  is empty. The window is **v1.4.3 bench `e4292f1` (2026-07-01) .. v1.9.1 `21b6fbe`**.
- process_features commits in that window (`git log e4292f1..21b6fbe -- core/features/process_features`):
  - `8c49e3e` 2026-09-10 Preserve PF sparse extents, direct input boundaries and failure cleanup
  - `b768003` 2026-09-10 Process leased BGZF batches with PF direct workers and remove duplicate view copies
  - `dd57aac` 2026-09-10 Add native BGZF ingestion to process_features with shared permits
  - `ec3d0ab` 2026-09-04 Support Ubuntu 22.04 Debian builds
  - `cc62337` 2026-07-24 Fix binary regression harnesses
- Also check the STAR-side feature driver (`core/legacy/source/PfMulti*.cpp`) in the same window; it was not
  listed here.
- New fields in the 1.9.2 assignment record that did not exist in April: `input_format=fastq`,
  `cbq_mode_requested=auto`, `readBufferLines=-1`, `bgzfModeRequested=auto`, `bgzfReaderThreads=32`,
  `bgzfCrcCheck=1`, `hash*`, `output_mode=default`. The MSK FASTQs are **plain gzip, not BGZF**.
- The assignment parameters are otherwise the same as April: `maxHammingDistance=1`,
  `featureConstantOffset=0`, `limitSearch=-1`, `featureN=0`, `barcodeN=1`, `consumerThreadsPerSet=31`,
  `searchThreads=1`, `minCounts=0`, `maxBarcodeMismatches=5`. `stats.*` fields are 0 in both runs (not filled).

## Inputs and commands

- 1.9.2 run: `/storage/paper_bench_v192_20260911/out/P02/` (`RUN_COMMAND.sh`, `multi_config.csv`,
  `Log.out`, `cr_assign/`). Driver logs: `/storage/paper_bench_v192_20260911/logs/P02.*`. Evidence copy:
  `/mnt/pikachu/star_suite_paper/analysis/paper_benchmark_refresh_20260910/v192/P02/`.
- Binary: STAR Suite 1.9.2, `/storage/paper_bench_v192_20260911/stage/STAR`, sha256 `def647c0b027e7d6...`,
  `--source-revision` `68eb92ff3cc107179d32ed52d95229b121debae9` (frozen build
  `/mnt/pikachu/benchmark_build_v192_20260911/src`, git archive of `v1.9.2`).
- Wrapper: `analysis/paper_benchmark_refresh_20260910/wrappers/run_msk_30polyko_benchmark.ES.v192.CORRECTED.sh`
  (corrections: ES paths, ES labels, 28-file gate, zcat reader for plain gzip). Run through
  `run_perturb_arm_v192.sh P02`.
- LARRY FASTQs (16 files, 50,681,432 reads, NVMe, sha256 verified at staging):
  `/home/lhhung/paper_bench_stage_20260911/MSK-perturb-comparison/msk30ko_full_3lib_ES_20260430_054814/fastqs/LARRY/`.
- LARRY feature reference: `/storage/paper_bench_v192_20260911/stage/msk/ref_feature_larryBC.csv`, 245,979
  rows, byte-identical to April's `/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv`. Every row:
  `read=R2`, `pattern=(BC)ATGTCTGGATCCGATATCGC`, a 40-nt `sequence`, `feature_type=Custom`. Whitelist
  `3M-february-2018_TRU.txt`, `star_max_hamming=1`.
- April run (accepted): `/mnt/pikachu/storage-relocated/20260528T165718Z/MSK-perturb-comparison/paper_bench_ES_20260430_054814/`
  (LARRY outputs in `cr_assign/Custom/larry_de/`, record `assignBarcodes.api_run.txt`). Built from
  `/mnt/pikachu/STAR-suite-v1.4.3-bench` (HEAD `e4292f1`); that binary still exists:
  `core/legacy/source/STAR`, version 1.4.3, sha256 `eb49c289...`.

## Suggested approach

1. Reproduce small: take the first ~1-2 M LARRY read pairs of one lane and run feature assignment alone
   (in-process through STAR, or the standalone assign tool) with the v1.4.3 bench binary and with 1.9.2.
   Record wall time and, if available, how many reads resolve by hash and how many fall through to the
   fastHamming/anchor search.
2. Test the owner's hypothesis directly: is the feature hash built completely for 245,979 x 40-nt features
   (entries, keys, anchor offset), and does a read that matches exactly get a hash hit?
3. If the hash is fine, bisect `e4292f1..21b6fbe` on process_features and `PfMulti*`; the three
   2026-09-10 commits are the first suspects.
4. Acceptance for a fix: LARRY assignment within 30% of 1,443 s on the full ES inputs with the same
   parameters, and a LARRY matrix identical to April's or with every difference explained. The PolyIII
   guide time should also return to within 30% of 146 s.

## Constraints

- Do not modify the frozen trees (`benchmark_build_v192_20260911`, `STAR-suite-v1.4.3-bench`) or the
  stage/output directories under `/storage/paper_bench_v192_20260911`; work in a new worktree.
- Commit messages carry no AI attribution. Do not push or release without the owner.
- Benchmark rule (owner): a run more than 30% slower than its previous accepted number is a regression;
  stop it and debug. One benchmark job at a time on this host.
- P02 is rerun only after the fix is released; the owner has asked that paper benchmarks use a released
  build.

## Other arms on 1.9.2 so far

- P01 (A375 GEX + CRISPR): 3:12.34 wall, 1,170 cells; comparator checks in progress.
- PBMC 10K Solo: 17:12.18 wall, 11,863 cells; comparator checks in progress.
- The remaining arms (B01, B02, B01-EXT, B02-EXT, L01, L02, F-arms) have not been run on 1.9.2.
