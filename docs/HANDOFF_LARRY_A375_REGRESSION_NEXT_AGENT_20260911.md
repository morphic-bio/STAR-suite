# LARRY/A375 process_features regression: next-agent handoff

Date: 2026-09-11. This is the current entry point. Detailed measurements and
artifact paths are in [the validation report](HANDOFF_PF_SEARCH_REVERT_A375_VALIDATION_20260911.md).

## Current outcome

A one-line search-policy revert is implemented and clean-built, but **the
regression is not closed**. It restores A375's historical assigned-read total
and nearly all historical feature counts; it does not recover integrated A375
feature-phase time. Full LARRY performance and real-read large-cache behavior
have not been validated with the patch. Nothing has been committed, merged,
pushed, or released during this work.

The owner wants one defect fixed at a time, with direct instrumentation first
and timings as a final test. Focus on process_features cache resolution. Do not
add unrelated optimizations or change permit scheduling on the strength of wait
counters: the owner specifically notes that waits may be a side effect.

## Checkout and exact change

- Active worktree: `/mnt/pikachu/STAR-suite-larry-regression-20260911`.
- Branch: `fix-larry-feature-regression-20260911`.
- Base: `68eb92ff3cc107179d32ed52d95229b121debae9`, v1.9.2.
- Production change: `core/legacy/source/PfMultiProcess.cpp:2429`:

```diff
- assignOpts.featureModeBootstrapReads =
-     (assignOpts.featureConstantOffset >= 0) ? 0 : 100000;
+ assignOpts.featureModeBootstrapReads = 100000;
```

This restores feature-position learning for explicit offsets, including the
offset 0 used by A375 and MSK. It is unrelated to cell-calling bootstrap/ties.
Environment overrides and the split-read `free` mode override remain intact.
The disabling change was introduced by `864b04f` on May 31, in the CBQ direct
range integration. Restoring learning also selects the existing BGZF ring
handoff, because learning is order-sensitive.

The only other working-tree changes are documentation: `docs/todos`,
`tests/ARTIFACTS.md`, and the two handoff files. The primary checkout
`/mnt/pikachu/STAR-suite` remains on master; unrelated edits there were preserved.

The production binary is `core/legacy/source/STAR`, SHA256
`c0e399de851638f6f3eeab6f176fe95b393dee19a0503636aa16d651a722d941`, reporting
`68eb92ff3cc107179d32ed52d95229b121debae9-dirty`. Both process_features and core
were cleaned before the Chromap-enabled build. A larger anchor-grouping
optimization was drafted earlier, then fully discarded at the owner's request;
it was never built or included in these results.

## Verified results

Full A375 GEX + CRISPR, 32 threads, BGZF auto, shared permits, no BAM:

| Measurement | Saved v1.9.2 | Policy revert | Historical March 26 |
| --- | ---: | ---: | ---: |
| Assigned feature reads | 9,063,551 | **9,274,337** | 9,274,337 |
| Deduplicated feature UMIs | 3,191,995 | **3,221,811** | 3,221,810 |
| Feature API time (s) | 58.4219 | 59.3113 | 11.8753 |
| Whole wrapper time (s) | 192.34 | 195.24 | — |
| Peak RSS (GiB) | 44.704 | 44.706 | — |
| GEX cells | 1,170 | 1,170 | — |

The wrapper completed with exit 0. Raw and filtered GEX MEX files are
byte-identical to saved v1.9.2. Semantic feature-MEX comparison finds one
historical difference: an additional UMI for `IL1B_sg2_HEK` at
`GGGTCAGCACCTCACT`, which is outside the current filtered GEX cell set. The
reason for that one UMI remains unknown. Versus v1.9.2, the patch changes 3,032
feature/barcode entries and adds a net 29,816 UMIs.

The owner also requested a quick **uninstrumented v1.4.3 feature-only** run.
Using the actual release tag `e4292f1`, all 9,748,584 feature read pairs, the
original gzip reader and 31 consumers, it completed in **9.66366 s API time**,
**11.30 s total**, with 6.74 s in read processing. Learning was 0, as in the
tag's explicit-offset policy. Its raw feature matrix exactly matches saved
v1.9.2 semantically. The earlier current-version feature-only control was also
fast: 7.4115 s API time using BGZF auto. These reader modes differ, so do not
claim a controlled version speedup ratio. The approximately 59-second integrated
A375 delay is not reproduced by either feature-only test; its cause is unproven.

## What the instrumentation establishes

Six instrumented feature-only runs processed the same first 200,000 A375 lane-1
read pairs with one consumer and gzip:

- Pre-May-29 source: `74b6e05^` = `f2ac1e2`.
- July v1.4.3 tag: `e4292f1`.
- Current v1.9.2 source: `68eb92f`.
- Each was run once with learning 0 and once with learning 100,000.

All six completed successfully. For each setting, **all recorded counters are
identical across the three revisions**. Learning 0 yielded zero offset-cache
hits and 207,881 successful resolver returns rejected by the enclosing
per-feature anchor-index check. Learning 100,000 yielded 1,537,396 offset-cache
hits; the trace includes bootstrap replay. See the detailed report for totals.

Critical limitation: **A375 makes zero full feature-hash lookups** in these
traces. Its 11 guides take the small-reference direct-comparison path, and two
non-common guide lengths disable Hamming prehash and the uniform-length tiered
path. It tests position-cache/acceptance behavior, but cannot clear LARRY's
245,979-entry exact hash or its Hamming-1 table.

Earlier LARRY diagnostics found:

- All 245,979 reference sequences self-hit the exact hash.
- The normal API path on the first 20,000 pairs, learning 0 and BGZF OFF, still
  produced **244,880,980 resolver calls**. Plain gzip was detected correctly.
- The required-anchor branch searches the shared suffix anchor repeatedly and
  only accepts a global resolver result when it equals the current loop's
  feature index. This can repeat even successful lookups.

These support a bypass/repeated-work mechanism, but are not a real-read H0/H1
hit-rate comparison against a confirmed successful historical build. Do not
claim the large-cache investigation is complete.

The owner reported that CATATAC and multiomics did not show the problem, then
noted that low guide counts might mask it. CATATAC's split-read `free` mode
explicitly disables required-anchor matching and follows another branch.
Multiomics' specific control run/path has not yet been identified.

## Source history and provenance corrections

The original handoff incorrectly treated the present v1.4.3 benchmark checkout
as proof of the binary used in April. Use execution logs and checksums, not a
checkout's current HEAD or directory date.

The exact and simple-Hamming resolver bodies are unchanged since March 17.
Relevant later changes are May 29 (`74b6e05`: serialized/synchronized position
learning and a match-position fix) and May 31 (`864b04f`: disable learning for
explicit offsets). The single-consumer traces do not test May 29 concurrency.
The detailed report records the earlier function-body diff audit.

I initially overgeneralized from a March archive in the public provenance
ledger. The owner correctly reminded me that later tests exist. Verified run
records include:

| Run | Recorded build | Relevance |
| --- | --- | --- |
| March 18 A375/MSK | Clean `ebc1362` | Older aggregate paper record; not the complete run history |
| March 26 A375 | `a70a039-dirty` | Historical feature-count/API-time comparator used here |
| April 30 MSK ES full | `cc84382-dirty` | 1,811 s / 33,226 cells; current paper's 30.2-minute Perturb row |
| June 30 MSK ES GEX-only | Clean v1.4.2, `996eea0` | 1,894 s / 33,346 cells; current paper's 31.6-minute Solo row |
| July 25 A375 compatibility test | v1.5.0, `26e6d05-dirty` | Later testing, but 16 threads and BAM; not an equivalent paper timing |

The public `STAR-suite-provenance/records/paper-2026/perturbseq-full-depth.json`
labels March-linked aggregates `1.7.0-candidate`. This metadata needs
reconciliation, but **does not prove later runs were absent**. That repository
is a curated public ledger, explicitly not a complete run registry. The
inspected paper wrappers retain run-local provenance without automatically
registering every execution centrally. No public record was edited.

The owner requested a separate follow-up in [docs/todos](todos): every run must
have provenance, with separate `dev/` and `production/` directories permitted.
It covers actual binary/revision/dirty patch, settings, input/reference/cache
identities, validation and terminal status, automatic registration with durable
retry/receipt, historical backfill, and exact paper-row/run links. This TODO is
recorded, not implemented. Continue the regression first.

## Artifacts and warnings about earlier attempts

All new artifacts are under `/home/lhhung/pf_larry_regression_20260911/`:

| Relative path | Contents / status |
| --- | --- |
| `a375_revert_star/` | Successful full patched STAR run; source patch, binary hash, wrapper, logs, time, completion, semantic feature and bytewise GEX comparisons |
| `a375_v143_feature_only/` | Clean uninstrumented v1.4.3 sources/build, full feature-only timing, command/provenance, matrix comparison |
| `cache_trace_a375/` | Three clean instrumented source snapshots, six successful runs, counters, input manifest and source comparisons |
| `setup_cache_trace.py` | Diagnostic-only instrumentation/setup code; not production source |
| `compare_a375_revert.py` | Semantic feature-MEX and bytewise GEX comparator |
| `api_current_off_20k/` | Clean LARRY API control, learning 0, BGZF OFF; repeated resolver-call evidence |
| `a375_feature_only_auto/` | Earlier clean current-version full feature-only control; learning 0, BGZF auto, no shared permits |
| `api_bootstrap_20k/` | LARRY mechanism diagnostic using **1,000** learning reads; not validation of the production 100,000 setting |
| `baseline_20k/` | Early CLI run omitted required-anchor flags; do not use as a matched STAR control |
| `historical_bootstrap_20k/` | Early CLI attempt timed out; not a valid historical control |

`baseline_anchor_20k/` is another CLI diagnostic; prefer the clean API controls
for interpreting STAR settings. No full LARRY run of the production patch was
launched. The subsequently proposed tiny instrumented LARRY test was not
launched before the owner redirected to the quick v1.4.3 A375 test and provenance.
No test started here is intentionally left running.

## Next work, in order

1. Verify the latest comparable successful full A375/LARRY run from actual
   provenance. The exact v1.4.3 release test is now available, but do not relabel
   older or later runs as that release without execution evidence.
2. Instrument the large-reference path on a tiny real LARRY subset: queries,
   H0 hits, H1 hits/ambiguities, misses, accepted hits, hits rejected by the
   enclosing anchor check, offset-cache use, and actual fallback entry. Include
   initialization and context-switch cache state if necessary. A reference
   self-hit check or A375's zero-hash path is insufficient to rule out the
   owner's hypothesis. Preserve the full LARRY reference for this diagnostic.
3. Localize the first divergence against a confirmed baseline with otherwise
   matched inputs/settings. Keep any correction narrow and preserve unrelated
   scheduling/reader behavior. The existing one-line revert remains a candidate
   fix with count validation, not a fully accepted performance fix.
4. Validate counts, then time the uninstrumented corrected path. Full LARRY
   acceptance in the original handoff is within 30% of 1,442.9 s, PolyIII within
   30% of 146.4 s, with historical count differences explained. The official
   full P02 paper rerun is held for a released fix.

## Constraints to retain

- Clean room: **never inspect Cell Ranger source**. Outputs, logs and documented
  interfaces may be used.
- Read `AGENTS.md` and the primary checkout's `AGENTS.local.md` for host rules.
  Clean builds are mandatory for regression diagnosis. Serialize runs and use
  fresh output directories; wrapper completion plus exit status is authoritative.
- Preserve frozen trees `/mnt/pikachu/benchmark_build_v192_20260911`,
  `/mnt/pikachu/STAR-suite-v1.4.3-bench`, and saved stage/results under
  `/storage/paper_bench_v192_20260911`.
- The local repeat rule forbids an absolutely identical execution without fresh
  authorization. Testing an authorized code fix after rebuilding is allowed.
- The owner requested a handoff, not a commit, merge, push, release, or provenance
  implementation. Retain the working-tree changes for review/integration.
