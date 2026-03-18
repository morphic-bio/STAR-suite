## Flex H0/H1 Deny-Class Analysis (2026-03-15)

Artifact root:
- `/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/deny_analysis_20260315_185922/`

Context:
- Full H0/H1 cache on the 100K `SC2300771` fixture.
- Cache-builder ordering bug was fixed before this analysis.
- Decision surface under test:
  - `KEEP` for cache classes `0/1`
  - `DENY` for cache class `2`
  - `PASS` for cache miss / fallback

Observed denied reads on 100K:
- `NEG_PROBE_AMBIG`: `5396`
- `NEG_NO_GENE`: `1873`
- `NEG_NO_CANDIDATES`: `136`
- `NEG_UNMAPPED`: `15`
- Total denied reads: `7420`

Rescue test:
- Each deny class was extracted into its own FASTQ subset and run through the same isolated Flex baseline binary at `16` threads.
- Raw MEX rescue totals:
  - `NEG_PROBE_AMBIG`: `12`
  - `NEG_NO_GENE`: `1672`
  - `NEG_NO_CANDIDATES`: `125`
  - `NEG_UNMAPPED`: `1`
  - Total deny rescue: `1810`

Interpretation:
- The full parity undercount from the `KEEP + PASS + DENY` split was `1809`, which matches the deny-rescue total almost exactly.
- `NEG_NO_GENE` and `NEG_NO_CANDIDATES` are not safe terminal classes. Most of their denied reads are recovered by the normal aligner path.
- `NEG_PROBE_AMBIG` is the only deny class that looks close to a true hard-drop candidate.

Baseline BAM observations:
- `NEG_NO_GENE` rescued reads are mostly ordinary `50M40S`, `NM=0`, probe/genome alignments.
- `NEG_NO_CANDIDATES` rescued reads show the same pattern.
- These do not look like simple low-quality or pathological CIGAR edge cases.
- This suggests the synthetic precompute is misclassifying these classes rather than the runtime aligner rescuing obvious garbage.

Conclusion:
- Keep `NEG_PROBE_AMBIG` as the only hard `DENY` candidate for the next iteration.
- Route `NEG_NO_GENE`, `NEG_NO_CANDIDATES`, and `NEG_UNMAPPED` to fallback alignment.

Follow-up test: `NEG_PROBE_AMBIG` as the only hard deny
- Artifact root:
  - `/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/probe_ambig_only_20260315_190952/`
- Decision split on 100K:
  - `KEEP`: `670450`
  - `DENY`: `5396`
  - `PASS`: `124154`
- Full raw MEX compare vs same-binary baseline:
  - `baseline_nonzero_entries`: `619230`
  - `exact_entries`: `618931`
  - `mismatch_entries`: `299`
  - Total counts: baseline `625225`, combined `625500`, delta `+275`
- Interpretation:
  - The prior deny-driven undercount is effectively gone.
  - Remaining drift is dominated by small positive overcounts, consistent with split-run UMI dedup across separate `PASS` and `KEEP` STAR jobs.
  - This makes `NEG_PROBE_AMBIG` the best current hard-deny policy if the fast path remains external to Flex.

## Internal Flex Hash-Screen E2E (2026-03-15)

Artifact roots:
- `/storage/downsampled_100K/SC2300771/results/flex_hash_screen_internal_20260315_193538/`
- `/storage/downsampled_100K/SC2300771/results/flex_hash_screen_internal_exactonly_20260315_194025/`

Context:
- The screen was moved inside Flex and enabled by default when a cache is available.
- `--no-hash-screen yes` forces the legacy align-everything path.
- Both runs used the same branch-local STAR binary and `16` threads.
- Cache:
  - `/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache.bin`

Timing:
- Internal hash-on, `KEEP` on cache classes `0/1`, `DENY` on `NEG_PROBE_AMBIG`:
  - elapsed `47.59s`
- Internal hash-on, exact-only `KEEP`, `DENY` on `NEG_PROBE_AMBIG`, pass otherwise:
  - elapsed `48.36s`
- Legacy internal control, `--no-hash-screen yes`:
  - elapsed `58.46s`

Correct comparison surface:
- Compare raw and filtered MEX on `(feature_id, barcode)` keys, not raw MatrixMarket column indices.
- The internal hash path changes barcode ordering, so direct `matrix.mtx` line diffs overstate drift.

Results vs internal legacy control:
- H0/H1 internal keep (`cacheClass 0/1` keep):
  - raw nonzero union entries: `620192`
  - mismatches: `962`
  - under: `113`
  - over: `849`
  - count delta: `+736`
- Exact-only internal keep (`cacheClass 0` only, offsets `0,+1,-1`, `DENY` secondary):
  - raw nonzero union entries: `620004`
  - mismatches: `902`
  - under: `110`
  - over: `792`
  - count delta: `+682`

Per-sample filtered MEX drift for exact-only variant:
- `BC004`: `141` mismatches, delta `+121`
- `BC006`: `203` mismatches, delta `+180`
- `BC007`: `294` mismatches, delta `+218`
- `BC008`: `279` mismatches, delta `+209`

Interpretation:
- Moving the screen inside Flex removes the split-run dedup artifact, but the current direct-insert fast path is still not byte-identical to legacy.
- The residual is now dominated by small positive overcounts and a small number of extra barcodes (`9` extra raw CB+TAG barcodes in the exact-only run).
- Tightening `KEEP` from `H0/H1` to exact-only helps, but only modestly (`962 -> 902` raw mismatches).
- The next implementation step should reuse the `process_features` pre-screen / constant-offset sampling primitives and instrument keep decisions by relative offset to isolate where the residual overcounts are entering.

## Internal Flex Hash-Screen E2E: Full 16-Sample Whitelist Fix (2026-03-15)

Artifact root:
- `/storage/downsampled_100K/SC2300771/results/flex_hash_screen_internal_full16_20260315_203745/`

Context:
- The earlier internal sample-aware runs used the filtered 4-entry sample whitelist:
  - `/storage/SC2300771_filtered_2M/sample_whitelist.tsv`
- That was the wrong sample-index surface for this experiment.
- The correct whitelist for the internal hash-screen run is the full 16-entry whitelist:
  - `/mnt/pikachu/flex/tables/sample_whitelist_full_16.tsv`
- The BC004 cache was restamped on that same 16-entry surface:
  - `/storage/downsampled_100K/SC2300771/results/flex_h01_full_cache_20260315_153914/reclassified/sequence_cache_bc004_full16.bin`

Results:
- Same-binary 100K internal E2E with `16` threads:
  - hash-on elapsed: `50.48s`
  - legacy elapsed: `59.67s`
- Keyed raw-MEX comparison:
  - `raw`: `619448` baseline nonzero entries, `619448` hash nonzero entries, `619448` exact entries, `0` mismatches
  - `ATGTTGAC / BC004`: `99286` exact, `0` mismatches
  - `ATCCCAAC / BC006`: `125830` exact, `0` mismatches
  - `AAGTAGAG / BC007`: `214330` exact, `0` mismatches
  - `AGCTGTGA / BC008`: `179784` exact, `0` mismatches
- Additional low-count tags present under the full whitelist also matched exactly:
  - `AACGGGAA`, `ACAGACCT`, `ACAGTCTG`, `AGAGGCAA`, `AGTAGGCT`, `AGTGAGTG`, `ATACGTCA`, `ATCATGTG`, `ATTCGGTT`, `AACGCCGA`

Conclusion:
- The remaining BC004 drift was not a hash-classification problem.
- It was a sample-whitelist surface bug: using the filtered 4-entry whitelist changed sample indices and broke the internal sample-aware screen.
- With the full 16-entry whitelist and the cache stamped on the same surface, the current internal hash-screen path is byte-identical to legacy on the 100K fixture.
