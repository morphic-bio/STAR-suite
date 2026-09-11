# Allocation audit: node-based maps and one-container-per-item patterns

Date: 2026-09-11. Audit only — no code was changed. Base: master `68eb92f` (v1.9.2).

## Why

Astra's `CbCorrector` change (uncommitted, `flex/source/solo/CbCorrector.*`) replaced
`std::unordered_map` exact/variant tables and an `unordered_map<uint32_t, vector<uint32_t>>`
of ambiguous variants — built for the 3,686,400-entry A375 whitelist — with khash and one
candidate array. That removed the ~27 s that every A375 run spent after `ALL DONE` freeing
those nodes, about 25% of A375 total wall. This audit looks for the same pattern elsewhere:

- **Node-based containers at data scale** (`std::map`, `std::set`, `unordered_map`,
  `unordered_set`) where khash or a flat array would do.
- **One small container per item** (a vector, hash table or string per barcode, UMI, cell or
  variant) where the total size is known or can be counted first, so one pre-sized buffer
  with offsets would do.

Two facts make both patterns expensive here:

- A 16-base barcode does not fit libstdc++'s 15-character inline string buffer, so every
  barcode held as `std::string` is a separate heap allocation.
- Long-lived and static containers are freed at exit, one node at a time. That cost appears
  after `ALL DONE`, outside every phase timer.

## Scope and method

Code compiled into the STAR binary, taken from the object files of a v1.9.2 build:
`core/legacy/source`, the VPATH sources from `flex/source`, `slam/source`, `vbem` and
`bamsort` (a same-named file in `core/legacy/source` takes precedence), process_features,
libscrna and libchromap_contract. Upstream STAR 2.7.11b lines were identified by exact
match and excluded unless noted.

The scan found 878 STAR Suite container declarations of these kinds; 482 are named after
data-scale items. Each of those candidates was read to establish key count, lifetime and
whether it runs by default. Counts below are from the v1.9.2 P01 run
(`analysis/paper_benchmark_refresh_20260910/v192/P01`) unless stated.

## Tier 1 — default path, data-scale

### 1. process_features: one khash table per unique barcode–UMI

- `core/features/process_features/src/assignBarcodes.c:2339-2372` (`update_umi_counts`).
- Every new barcode–UMI pair gets a pooled struct (good) plus its own `u32u32` khash
  (`s->counts = kh_init(u32u32)`). The table holds the feature and a key-0 "visited" marker
  used by connected components.
- Each table costs about four heap allocations (the struct from `kh_init`, then flags, keys
  and values on the first insert). They are freed one at a time by
  `destroy_feature_umi_counts` (`assignBarcodes.c:1133-1140`).
- **Scale:** A375 has 3,296,220 unique barcode–UMIs (P01 `stats.txt`). The LARRY 20K-read
  subset has 15,316 unique barcode–UMIs from 16,366 assigned reads — about one table per
  read. Full LARRY will be in the tens of millions.
- Same pattern per barcode: `s->counts` at `assignBarcodes.c:2293` and `:2357` (108,922
  barcodes in A375). The merge (`:6346`, `:6390`) copies each table entry by entry into a
  new table.
- **Direction:** keep a fixed slot (feature, count) and a visited bit inside the pooled
  struct, and move only the rare multi-feature UMIs to a shared table.

Related: tables created with `kh_init` and never pre-sized, although the size is known or
can be estimated: `sequence_umi_hash` and `filtered_hash` (`assignBarcodes.c:6759-6761`), and
`whitelist_hash` (`pf_api.c:925`; the whitelist line count is known). STAR already
estimates reads per library (`estimateReadsAnchored`, `PfMultiProcess.cpp:882`).

### 2. Feature gather loads the whole GEX whitelist as `unordered_set<string>`

- `PfMultiProcess.cpp:1174` (`loadWhitelistBarcodeSet`), called at `:3656-3662` whenever GEX
  filtered barcodes come from Solo. That is the default; P01 logs "Using GEX filtered barcodes
  from Solo (normalized to TRU, in-memory)" at `Log.out:1445`.
- It builds 3,686,400 string nodes with no `reserve`, each needing a heap string, to
  normalize about 1,170 called barcodes.
- It runs inside the serial feature gather: P01 merge/finalize took 9 s (09:40:53–09:41:02).
  This step was not timed separately.
- Also at `:2736-2739`, when explicit filtered barcodes are given.
- **Direction:** look the ~1K barcodes up in the existing Solo `cbWLhash` (khash of packed
  barcode to whitelist index, `ParametersSolo.h:119`, already built with 8,388,608 buckets).

### 3. Whitelist re-read to check for a second column

- `PfMultiMexStub.cpp:43` (`parseWhitelistOutputMap`), called from `copyBarcodesTsv` (`:299`)
  and `remapFeaturePerCellCsv` (`:331`).
- `copyBarcodesTsv` is on the default deferred filtered-MEX path (`PfMultiProcess.cpp:2244`;
  P01 `Log.out:1451`).
- For a one-column whitelist it reads all 3.69M lines, making a substring per line, and
  inserts nothing. pf-multi already logs `hasTwoColumnWL=no` for A375.
- **Direction:** pass the known column count, or decide from the first line.

### 4. The same whitelist held several times as heap strings

In an A375 run:

- `pSolo.cbWLstr` — upstream, 3.69M strings.
- `CbCorrector::whitelist_` — a full copy (`CbCorrector.cpp:12`, `whitelist_ = whitelist;`).
  The member is still present in Astra's working-tree header.
- The gather set from item 2.
- `InlineCBCorrection::whitelistHash_` when inline correction is on (item 5).

The process_features copy is a single block (good).

**Direction:** one packed whitelist shared by reference. Solo already holds `cbWL` (uint64)
and `cbWLhash`; keep strings only for output.

## Tier 2 — conditional paths or moderate size

5. **`InlineCBCorrection`** (`flex/source/InlineCBCorrection.cpp:12-16`, `:320-376`). Off by
   default (`soloInlineCBCorrection no`, `parametersDefault:2036`), but presets can switch it
   on (`ParametersSolo.cpp:574-584`). P01 did not run it. Its main tables are pre-sized khash
   (good), but it also builds:
   - `whitelistHash_`, an `unordered_set<uint64_t>` duplicate of `exactHash_` covering the
     whole whitelist ("kept for legacy uses (N-path)"), not reserved;
   - `variantCollisions_`, an `unordered_map<uint64_t, vector<uint32_t>>` — the exact shape
     Astra replaced in `CbCorrector`.

   Both are static, so they are freed at exit.

6. **Ambiguous cell-barcode entries.** Each entry holds a candidate vector, an
   `unordered_map` of UMI counts, two 16-character `std::string`s and two `vector<double>`s:
   - `ReadAlign.h:150-160` (`AmbiguousEntry`);
   - `SoloReadFeature.h:62-72` (plus an observations vector) and `:130-139` (orphan entries);
   - `InlineCBCorrection.h:121-129` (`MergedAmbigEntry`).

   They sit in per-thread `unordered_map`s that are merged later. A375 is small (about 1,000
   per resolution pass, per P01 `Log.out`). The Flex-scale count was not measured.
   **Direction:** fixed arrays for the barcode sequence and qualities, candidates in a shared
   pool, UMI counts in a flat vector sorted at merge.

7. **`bridgeImmediateReadCounts_`** (`SoloReadFeature.h:141`). An `unordered_map` keyed by
   whitelist index, which is dense and bounded; 294,633 observed barcodes in A375. A dense
   array or khash would do. The snapshot code already flattens it to pairs.

8. **One vector per barcode when the counts are already known.**
   `SoloFeature_collapseUMI_fromBridgeHash.cpp:693-694` has `cbRowCounts(nCB)` next to a
   `vector<vector<uint32_t>> cbRows(nCB)`. A375 has nCB=294,633 and 7,587,372 rows, and the
   step takes 1.58 s total, so impact is low on A375. **Direction:** one array plus offsets
   (CSR).

9. **Velocyto only.** `SoloFeature_countVelocytoBridge.cpp:281, 358, 440, 681` builds
   `vector<unordered_map<uintUMI, vector<trTypeStruct>>>(nCB)`, a per-cell map of per-UMI
   vectors, plus a per-cell `std::map` (`:718`). Reserve helpers exist, but the structure is
   still millions of small containers at scale.

10. **Flex tag occupancy (serial cell-calling gather).** Up to 320K called cells on GEM-X:
    - `libscrna/src/OccupancyGuard.cpp:229-230`: `unordered_map<string, unordered_set<string>>`
      and `unordered_map<string, vector<uint32_t>>`, keyed by 16-base barcode strings;
    - `libscrna/include/ObservedTagOccupancy.h:19-27`: `std::map<string, std::set<string>>`,
      used from `FlexFilterTagAware.cpp:156`.

    **Direction:** packed barcode keys, and tags as a bitmask (at most 16 tags).

11. **FlexFilter file input** (`libflex/FlexFilter.cpp:1397`, via `runFromFiles`, `:1484`).
    It builds `vector<vector<pair<>>>(numCells)` one triplet at a time, although the MTX
    header gives the entry count; two passes into an array plus offsets would do. This is the
    file-input path, not the in-process count path.

12. **Merge and import string maps.** `PfMultiMerge.cpp:449, 646, 765` (reserved) and the
    `PfMultiTableImport.cpp` barcode maps key raw barcodes as strings (about 294K in A375).
    They cost gather-step time. **Direction:** packed keys.

## Tier 3 — low impact or no action

- **Debug and trace only:** `g_iReadToQname` (`flex/SoloReadFeature_record_flex.cpp:107`,
  used only with reject-log qname tracing), the trace barcode and read sets in
  `SoloFeature_collapseUMIall.cpp:30-108` and `SoloReadFeature_record_base.cpp:347-615`, and
  the `SlamQuant.h:85-87` distribution maps.
- **Small sizes:** parameter and config maps, per-tag maps, `FlexProbeIndex` gene sets (built
  once), and `SpatialR1Decoder` variant generation (reserved, build time only).
- **Multiome ATAC peak MEX:** has string-keyed barcode maps
  (`multiome_atac_peak_mex.cpp:940, 1005, 1089, 1129`), but integrated STAR uses the binary
  sidecar path with uint64 keys (`star_chromap_orchestration.cpp:710, 736`).
- **`bucketSegments_`** (`SoloReadFeature.h:192`): a few large buckets handed off by swap.
  Fine as is.
- **Upstream STAR structures** (for example `umiGeneMapCount` in `collapseUMIall`,
  `cbFeatureUMImap`): not used in inline-hash runs ("inline-hash mode completed (skipping
  legacy output)").

## Checks for every replacement

- **Iteration order.** `std::map` iterates in sorted order; unordered containers and khash
  iterate in hash order. Any output or tie-break that iterates the container must sort
  explicitly first. Require byte-identical matrices, cell calls and guide calls.
- **Timing.** Measure the gap from `Log.final.out` "Finished on" to process exit, then each
  phase including its gather step. The ranking is total wall, then tail, then average CPU.
- **Scale.** Log the key count next to each replaced structure, so the next audit has
  measurements instead of estimates.

## Not covered

- **Per-read temporaries** (`std::string` or `vector` built per read in hot loops) are a
  separate pattern that needs its own pass.
- **Rest of process_features:** only the khash sites were read, not the other allocation
  sites (`call_features.c`, `nbem.c`, `gmm.c`).
- **Other components** (`vbem`, SLAM hot loops, libscrna internals beyond container
  declarations) were not read.
- **Flex-scale counts** for items 5, 6 and 10 were not measured.

## Addendum after the MSK P02 test of 86795d1 (2026-09-11)

A new finding, not covered above, and pre-existing in v1.4.3, v1.9.2 and 86795d1.

- **Where:** `core/features/process_features/src/assignBarcodes.c:2525` (`find_deduped_counts`),
  line numbers as of 86795d1.
- **What:** `find_deduped_counts` clears `clique_counts[number_of_features+1]` once per unvisited
  barcode–UMI group. It then scans every feature again for counts above `min_counts`
  (`:2674`).
- **Scale:** LARRY has 245,979 features, so each group costs about 1 MB of memset plus a full
  scan. That is roughly 15M groups on the full MSK ES LARRY library.
- **Measured:** in the 86795d1 MSK test, `umi_dedup` took 1,079 s on one thread — 75% of LARRY's
  1,440 s — while read/assign took 339 s. A gdb sample during the phase was in this
  `memset` (983,920 bytes).
- **Direction:** clear and scan only the features the group touched (a short touched-list or
  sparse accumulator), or keep the clique counts in a small hash for large references.

Byte-identical LARRY output is required: the current LARRY matrix matches April exactly.
