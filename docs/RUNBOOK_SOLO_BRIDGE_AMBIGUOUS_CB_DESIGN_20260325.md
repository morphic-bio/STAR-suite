# Runbook: Solo Bridge Ambiguous-CB Design

Date: 2026-03-25  
Audience: ongoing Solo optimization work on the non-Flex direct bridge and spool->HASH_MEX follow-on paths

## Purpose

This runbook captures the current design decisions and working assumptions for
the bridge-side Solo mapping data structures, with primary emphasis on the
ambiguous-CB path but also covering the simplified exact-read structure that
should accompany it.

This is a design-tracking document, not a benchmark summary. Benchmark roots
and measured timings should continue to live in:

- [HANDOFF_SOLO_OPTIMIZATION_20260324.md](/mnt/pikachu/STAR-suite/docs/HANDOFF_SOLO_OPTIMIZATION_20260324.md)
- [HANDOFF_SOLO_OPTIMIZATION_NEXT_AGENT.md](/mnt/pikachu/STAR-suite/docs/HANDOFF_SOLO_OPTIMIZATION_NEXT_AGENT.md)
- [ARTIFACTS.md](/mnt/pikachu/STAR-suite/tests/ARTIFACTS.md)

## Current Problem

Profiling of the non-Flex direct bridge mapping path on the full corrected UCSF
GEX-only run showed that the largest measured bridge-local mapping cost is the
deferred multi-CB bookkeeping path:

- [captureBridgeDeferredAccounting()](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature_record_base.cpp#L87)
- storage in:
  - [bridgeDeferredAccounting_](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature.h#L71)
  - [bridgeDeferredCandidates_](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature.h#L72)

This path currently fires whenever the direct bridge is enabled and
`cbMatch > 1`:

- [bridgeUsesDeferredAccounting()](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature_record_base.cpp#L30)

The current shape is expensive because it stores per-read deferred metadata
beside the already-existing ambiguous-key aggregate in:

- [pendingAmbiguous_](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature.h#L61)

## Working Design Direction

### 0. Simplify the exact-read structure at the same time

The current direct bridge exact-read path uses:

- a main tuple hash:
  - packed `(compactCB, umi, compactGene) -> slot id`
- a packed payload vector:
  - `slot id -> packed(umi, gene, count, flags)`
- two per-thread compaction hashes:
  - whitelist CB index -> compact CB id
  - full gene id -> compact gene id

Current working decision:

- remove redundant payload duplication
- remove per-thread compaction hashes
- use global ids directly in the main exact-read hash key

Target exact-read shape:

- main hash key:
  - packed `(wlCb, umi, geneFull)`
- hash value:
  - count only
  - optionally packed with overflow/flags if needed

This means the bridge should stop storing `umi` and `gene` again in a payload
vector when they are already encoded in the tuple key.

Rationale:

- the payload vector/grouping work was not the top mapping bottleneck
- CB compaction and gene compaction were both measurable mapping costs
- whitelist CB index and gene id are already globally stable identifiers
- if the packed key budget works with global ids, there is no need to pay for
  dynamic per-thread remapping

Corollary:

- [bridgePackedSlots_](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature.h#L83)
  should not be part of the long-term exact-read design
- [bridgeCbCompactByWl_](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature.h#L77)
  and
  [bridgeGeneCompactByFull_](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature.h#L79)
  should not be part of the long-term exact-read design either, unless later
  measurement proves that direct global-id keys are too large

### 1. Aggregate by ambiguous key, not by ambiguous read

The main design decision is to push as much state as possible into the
ambiguous-key aggregate and remove the per-read deferred side logs.

Target shape:

- one hash keyed by ambiguous key
- aggregate state attached to that key
- no `bridgeDeferredAccounting_`
- no `bridgeDeferredCandidates_`

### 2. Reuse the ambiguous-key aggregate as the only in-mapping structure

The bridge already has the right conceptual owner:

- [ReadAlign::AmbiguousEntry](/mnt/pikachu/STAR-suite/core/legacy/source/ReadAlign.h#L136)
- [SoloReadFeature::ExtendedAmbiguousEntry](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature.h#L50)

The design goal is to extend that aggregate rather than create another
sidecar log.

### 3. Prefer a representative best-quality CB string over additive QC evidence

Earlier discussion considered storing aggregated per-position log evidence.
That is mathematically clean, but it may overweight PCR duplicates by treating
them as independent observations.

Current working decision:

- keep one representative `cbQual` per ambiguous key
- define "best" by a deterministic total quality score
- when another read for the same ambiguous key arrives:
  - if the new `cbQual` scores better, replace the stored one
  - otherwise keep the existing one

This is intended as a pragmatic biological approximation:

- repeated PCR duplicates should not keep increasing confidence
- the best observed CB quality string is a better estimator of the underlying
  barcode accuracy than the sum of many likely-duplicate low-information reads

### 4. UMI histogram is not currently required for CB resolution

The current resolver accepts a UMI histogram:

- [CbBayesianResolver::resolve()](/mnt/pikachu/STAR-suite/flex/source/solo/CbBayesianResolver.cpp#L105)

However, the current UMI likelihood term is candidate-invariant:

- [CbBayesianResolver.cpp](/mnt/pikachu/STAR-suite/flex/source/solo/CbBayesianResolver.cpp#L160)

Therefore, under the current resolver logic, the UMI histogram does not help
choose between candidate CBs. It can be removed from the in-mapping ambiguous
aggregate unless later work intentionally introduces a candidate-dependent UMI
model.

Working decision:

- do not carry a dynamically growing `umiCounts` histogram solely for current
  CB disambiguation
- any UMI-related information needed for later tuple emission should come from
  the aggregated `(ambigKeyId, umi, gene)` observation store

### 5. Compact 64-bit observation key is preferred for ambiguous observations

Current preferred packed key for ambiguous observations:

- `ambigKeyId` : 16 bits
- `umi` : 32 bits
- `gene` : 16 bits

Packed into one 64-bit key:

- `[ambigKeyId:16][umi:32][gene:16]`

Rationale:

- avoids oversized custom keys
- avoids storing 2-bit packed CB sequence in the hot observation key
- ambiguous-key cardinality staying below 65,535 is considered safe for this
  workload
- 16 bits for gene is acceptable for current intended scope

This implies two tables:

1. `AmbigKey -> ambigKeyId`
2. `ambigKeyId -> metadata`

And one main ambiguous-observation aggregate hash:

3. packed `(ambigKeyId, umi, gene) -> count / best-obs payload`

### 6. Per-key metadata should remain small

Per ambiguous-key metadata should contain only:

- candidate list
- observed CB sequence
- representative best `cbQual`
- optional deterministic best-quality score
- aggregated read-accounting totals needed after resolution

Avoid:

- per-read deferred logs
- per-key dynamic UMI histogram unless the resolver is changed to actually use
  it

### 7. Prefer one authoritative hash per logical purpose

Current working principle:

- exact resolved reads:
  - one tuple hash keyed by `(wlCb, umi, geneFull)` with count as value
- ambiguous reads:
  - one ambiguous-key aggregate
  - one packed ambiguous-observation hash keyed by `(ambigKeyId, umi, gene)`

Avoid layering redundant structures unless profiling shows a specific need.

In particular, avoid:

- tuple hash + self-describing payload vector for the same exact tuple
- tuple hash + dynamic compaction hashes where global ids already exist
- ambiguous-key aggregate + separate per-read deferred logs

## Proposed Replacement for Current Deferred Bookkeeping

Replace:

- [bridgeDeferredAccounting_](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature.h#L71)
- [bridgeDeferredCandidates_](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature.h#L72)

With:

- ambiguous-key aggregate metadata
- packed ambiguous observation hash keyed by `(ambigKeyId, umi, gene)`

At resolve time:

1. resolve ambiguous key to winning CB
2. apply aggregated read-accounting totals to the winner
3. emit aggregated `(CB, umi, gene)` tuples from the packed observation hash

This keeps mapping-time cost proportional to:

- number of ambiguous keys
- number of unique `(umi, gene)` observations per ambiguous key

instead of:

- number of ambiguous reads
- number of candidate entries appended per ambiguous read

## Proposed Exact-Read Structure

Replace the current exact-read bridge shape:

- packed tuple hash -> slot id
- slot id -> packed payload vector
- CB compaction hash
- gene compaction hash

With:

- packed exact tuple hash keyed by `(wlCb, umi, geneFull)`
- count as hash value

If overflow/flags are needed, store them in the hash value next to count rather
than duplicating tuple fields into a separate slot payload.

This is intended to reduce:

- dynamic per-thread remap cost
- extra vector growth
- duplicate tuple metadata storage
- downstream dependency on self-describing slots

## Open Questions

### Gene width

Current working assumption is that `16` bits for gene is acceptable for the
intended GEX scope. If that assumption changes, the packing plan must be
revisited.

### Exact tuple key packing budget

The preferred exact-read design uses direct global ids, not compacted ones.
This requires the packed key budget for `(wlCb, umi, geneFull)` to remain
practical. If that budget becomes too tight, revisit static/global remapping
before reintroducing per-thread dynamic compaction hashes.

### Deterministic tie-break for best `cbQual`

Need a fixed policy if two quality strings score equally. Acceptable examples:

- keep existing
- lexicographically smaller `cbQual`
- smaller read id if already available cheaply

This should be explicit in implementation and docs.

### Whether read-flag stats can be fully aggregated per ambiguous key

The design expects that the necessary read-accounting stats can be accumulated
at the ambiguous-key level instead of replaying per-read deferred records.
This should be verified carefully during implementation.

## Next Implementation Target

The next coding step should focus on:

1. removing per-read deferred bridge bookkeeping
2. extending the ambiguous-key aggregate with:
   - best-quality CB string
   - aggregated read-accounting totals
3. introducing a compact packed observation hash keyed by
   `(ambigKeyId, umi, gene)`
4. simplifying the exact-read path toward:
   - direct global-id tuple key
   - count as hash value
   - no redundant payload vector
   - no per-thread compaction hashes
5. preserving current resolver behavior except for the planned best-quality
   representative policy

## Non-Goals for This Design Slice

- redesigning the post-map HASH_MEX path
- changing the downstream `1MM_CR` / `MultiGeneUMI_CR` logic
- switching the current resolver to a different probabilistic model
- adding a candidate-dependent UMI prior model

## Full-sample benchmark — EBs2_2 tuple redesign (2026-03-25)

**Intent:** measure the current non-Flex bridge tuple + per-key ambiguous
aggregate design (with **Bayes-resolved ambiguous accounting**) on the corrected
full UCSF **EBs2_2** GEX input family, using the **same STAR CLI harness** as
the archived direct-bridge full run `ucsf_ebs2_2_standard_solohash_optimized_v4`
(GEX-only `readFilesIn`; `pfMultiConfig` + CR-assign flags unchanged).

- **Worktree binary:** `/tmp/star-suite-v10-redesign-20260325` after
  `make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR`
- **Environment:** `STAR_SOLO_NONFLEX_HASH_BRIDGE=1`
- **Artifact root:**
  `/storage/100K/ucsf_solo_bridge_redesign_20260325/ucsf_ebs2_2_gexonly_bridge_redesign_v1/`
- **Driver:** `RUN_COMMAND.sh` in that directory (also `time.txt`,
  `BENCHMARK_SUMMARY.txt`, `pf_multi_config.csv`)

**Captured metrics (this run):**

- Wall **20:07.92**; max RSS **70985716 kB**
- Mapping **10:19:27 → 10:28:18** (**531 s**)
- `collapseUMIall_fromBridgeHash` **277.749 s**; `countCBgeneUMI` **279.954 s**;
  `processRecords` **505.697 s**
- `Summary.csv`: GeneFull unique **0.85468**; cells **13715**; UMIs **254181669**;
  Total GeneFull **33750**

**vs prior full-sample anchors (handoff):** bridge **v4** wall **22:38.81** /
RSS **70575848 kB** / mapping log **678 s**; baseline mastermerge v2 wall
**19:48.45** / RSS **68149756 kB**. On this execution the redesign sat **between**
baseline and v4 on wall time and showed a **much shorter** mapping log interval
than v4; **collapse** time stayed ~**278 s** (unchanged bottleneck).

Serialize full-sample benchmarks (do not overlap with other bench jobs on the
same host). See [HANDOFF_SOLO_OPTIMIZATION_20260324.md](HANDOFF_SOLO_OPTIMIZATION_20260324.md)
for tables and interpretation.
