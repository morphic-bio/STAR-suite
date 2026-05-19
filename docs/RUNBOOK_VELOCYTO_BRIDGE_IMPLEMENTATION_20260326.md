# RUNBOOK: Velocyto Bridge Implementation and UCSF 2M Validation

Date: 2026-03-26
Branch intent: `feature/velocyto-optimizations-20260326`

Related docs:
- `docs/RUNBOOK_VELOCYTO_COUNT_RESOLUTION_20260326.md`
- `docs/RUNBOOK_VELOCYTO_OPTIMIZATION_PORT_20260326.md`
- `docs/RUNBOOK_SCRNA_OCM_MULTI_MEX_MATERIALIZER_IMPLEMENTATION_20260519.md`
- `docs/HANDOFF_SCRNA_DOWNSTREAM_MEX_VELOCYTO_FINDINGS_20260313.md`
- `tests/run_perturb_velocyto_mex_smoke.sh`
- `tests/run_ucsf_velocyto_gexonly_exact_100k.sh`
- `tests/external_fixtures_env.sh`

## Branch status (do not overstate progress)

**Stage 1** (sorted global replay + stream) and **Stage 2** (CB-bucket deterministic merge, env `STAR_VELOCYTO_INTEGRATED_HASH=1` with `STAR_VELOCYTO_DETERMINISTIC_REPLAY=1`) are implemented in `SoloFeature_countVelocytoBridge.cpp` / `SoloFeature_countVelocyto.cpp`. UCSF validation uses `tests/run_ucsf_velocyto_exact_*.sh` (Stage 1 gates + Phase 6 `Gene`/`GeneFull` parity via `compare_velocyto_mex.py --mode genes`) and `tests/run_ucsf_velocyto_hash_*.sh` (Stage 2 vs Stage 1 + thread parity + the same Gene/GeneFull checks). **Phase 4** (external `velocyto.py`) remains blocked without full BAMs (`tests/run_ucsf_velocyto_external_compare.sh`).

**Status:** Implementation / code review for this bridge is effectively **closed**. **Formal acceptance** stays **open** until a host has valid UCSF fixture paths (or overrides below) and a chosen `UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB` for the **2M** exact and hash harnesses (unless using the dev-only `UCSF_VELOCYTO_ALLOW_UNCAPPED_2M=1` opt-out).

**OCM materialization note:** New OCM production runs should not call the legacy
`prepare_velocyto_mex.py` helper. STAR now packages run-level raw/filtered
Velocyto MEX through `VelocytoMexWriter`, and the OCM materializer writes
per-sample `raw_velocyto_feature_bc_matrix` and
`filtered_velocyto_feature_bc_matrix` mirrors. The per-sample subsetting maps
GeneFull columns to Velocyto columns by normalized barcode key, so it does not
depend on GeneFull and Velocyto barcode order being identical.

## Goal

Implement exact legacy-style `Velocyto` support on the optimized Solo path in
two stages:

1. first as a separate post-mapping / post-recording path that does **not**
   interfere with the normal `Gene` / `GeneFull` optimized count flow
2. then as an integrated hash-backed path that is shaped like the optimized
   `GeneFull` count rescue / collapse path

Validation must include the UCSF perturb fixtures, culminating in UCSF 2M
testing against:

- legacy STAR `Velocyto`
- packaged raw/filtered `Velocyto` MEX output
- real external `velocyto.py` output when full retained BAMs are available

## Preferred First Debug Loop

Before going back to the full perturb / full-sample harnesses, use the smallest
useful correctness surface:

- **GEX-only**
- **corrected `EBs2_2/GEX` source**
- **100K downsample**
- **1-thread vs N-thread only on that small fixture**

Canonical small-loop harness:

- `tests/run_ucsf_velocyto_gexonly_exact_100k.sh`

That harness:

1. creates a small GEX-only fixture from corrected `EBs2_2/GEX`
2. runs `stream_t1`
3. runs deterministic replay `det_t1`
4. runs deterministic replay `det_tN`
5. checks exact `Velocyto` packaged/raw parity plus `Gene` / `GeneFull` parity

This is the preferred first step for correctness debugging because it avoids:

- feature-calling / `pfMultiConfig` confounders
- full-sample wall times
- paying the cost of 1-thread full-sample runs just to debug exactness

## Non-Goal

This document is not redefining `Velocyto` semantics.

Exact legacy behavior is already pinned down in:

- `docs/RUNBOOK_VELOCYTO_COUNT_RESOLUTION_20260326.md`

This runbook is about implementation order, validation surfaces, and what to
compare at each stage.

## Exact Legacy Target

Legacy STAR uses:

- STAR transcriptome ids and transcript-to-gene mapping from `Transcriptome`
- per-read `trVelocytoType` entries
- a `(CB, UMI) -> vector<trTypeStruct>` aggregate during counting

The exact working state is:

- key: `(CB, UMI)`
- value: surviving transcript set
- payload per transcript:
  - `transcript_id`
  - `4-bit mask`

The exact update rule is:

1. for the first read of `(CB, UMI)`, store the full transcript vector
2. for later reads with the same `(CB, UMI)`:
   - intersect transcript ids
   - keep only survivors
   - merge masks on survivors with bitwise `OR`

The exact finalize rule is:

1. if the surviving transcript set is empty, drop the UMI
2. if surviving transcripts map to multiple genes, drop the UMI
3. otherwise collapse the merged masks for that single gene to exactly one
   class:
   - `spliced`
   - `unspliced`
   - `ambiguous`
4. add `+1` UMI to that `(gene, class)` bucket

Important constraint:

- this logic is transcript-level until finalization
- it is not recoverable from final `Gene` or `GeneFull` matrices alone

## Implementation Strategy

### Stage 1. Separate exact sidecar path

Build a new exact bridge-compatible `Velocyto` path that is logically separate
from the optimized normal count path.

The purpose of Stage 1 is:

- prove exactness
- keep the implementation easy to reason about
- avoid destabilizing `Gene` / `GeneFull`
- establish UCSF 2M validation before optimizing the data structure

Recommended shape:

1. preserve or derive per-read `trVelocytoType` exactly as legacy does now
2. maintain a dedicated `(CB, UMI) -> vector<trTypeStruct>` sidecar aggregate
3. finalize into raw `spliced`, `unspliced`, `ambiguous` matrices
4. package raw/filtered MEX the same way the current helper path does

This stage should deliberately avoid:

- sharing collapse logic with `GeneFull`
- reusing `GeneFull` counts as a shortcut
- any approximation from `Gene` / `GeneFull` output

Why this stage comes first:

- it isolates semantics from optimization
- it gives a direct parity target against legacy STAR `Velocyto`
- any mismatch is easier to debug before the data structure is fused with the
  optimized hash path

### Stage 2. Hash-backed integrated path

Only after Stage 1 matches legacy behavior should the state be integrated into a
hash-backed optimized path shaped like the normal count flow.

Target shape:

- build the `Velocyto` aggregate while the optimized path is already handling
  per-read / per-UMI information
- keep the `Velocyto` state separate in semantics, but colocated in the hot path
- convert global regrouping into CB-first contiguous work where practical

The working design should still preserve exact information content:

- `(CB, UMI)` key
- surviving transcript vector
- per-transcript 4-bit mask

Counts are optional debug / instrumentation state only. Exact classification
depends on transcript survival and merged masks, not transcript abundance.

The integrated path should mimic the same structural wins used successfully in
the optimized `GeneFull` / Flex work:

- CB-first contiguous spans
- direct emit where possible
- no repeated materialization of equivalent intermediate structures

But optimization comes only after exact Stage 1 parity is proven.

## Recommended Build Order

### Phase 0. Baseline and instrumentation

Before touching code:

1. confirm the current legacy `Velocyto` path still runs on the UCSF perturb
   fixtures
2. record wall time, peak RSS, and phase timing
3. save one canonical set of legacy outputs for later exact diff comparison

Canonical local anchors:

- **Primary UCSF acceptance / benchmark surface (repo default today):** corrected
  full-sample **EBs2_2** perturb layout — same family as
  `scripts/paper/run_ucsf_ebs2_2_benchmark.sh` and `README.md` CR9 parity rows:
  - `DATASET_ROOT=/mnt/pikachu/ucsf-perturb-seq-corrected/EBs2_2`
  - `GEX` / `guides` under that root, plus a `pf_multi_config.csv` that lists
    those directories (see the paper script’s emitted multi-config shape).
- UCSF fixture env defaults (may be stale on some hosts; **override explicitly**
  for Velocyto gates):
  - `tests/external_fixtures_env.sh`
- UCSF **100K-only** perturb Velocyto MEX smoke (packaging invariants):
  - `tests/run_perturb_velocyto_mex_smoke.sh` (optional fast loop when a staged
    100K tree exists; **not** the same as the EBs2_2 acceptance surface)
- Canonical STAR parameter surface + exact/hash harnesses (`100k` and `2m`
  profiles are harness names; `2m` is used for **full-sample** EBs2_2 runs):
  - `scripts/run_star_velocyto_canonical.sh`
  - `tests/run_ucsf_velocyto_exact_100k.sh`, `tests/run_ucsf_velocyto_exact_2m.sh`

**Legacy / historical (do not treat as canonical for new acceptance):** the old
`/storage/ucsf-2M/.../iPSC2_1_AALG2` and bundled 100K pfconfig under
`/storage/ucsf-2M/fixtures/ucsf2m_iPSC2_AALG2_100k_pfconfig` — kept only for
`external_fixtures_env.sh` defaults and older notes. Validators targeting current
paper-style UCSF work should use **corrected EBs2_2** paths and explicit
`UCSF_2M_*` exports (see “Minimum fixture env” below), not the stale tree.

### Phase 1. Exact separate sidecar

Implement the new exact sidecar path with no attempt at optimization.

Minimum requirements:

1. derive the same per-read transcript / mask records as legacy
2. aggregate them by `(CB, UMI)`
3. apply transcript-set intersection plus mask `OR`
4. finalize to `spliced`, `unspliced`, `ambiguous`
5. emit raw matrices on the full raw barcode axis
6. gate filtered output with the same `Gene` cell filter behavior as legacy

Validation at this phase:

- exact matrix parity vs legacy STAR `Velocyto` on a small reproducible fixture
- same feature axis
- same raw barcode axis
- same filtered barcode subset
- same `spliced`, `unspliced`, `ambiguous`, and total matrices

### Phase 2. UCSF 100K exact parity (optional)

When a staged 100K perturb fixture is available, use it for fast iteration. Skip
this phase if the host only has the full-sample tree (go to Phase 3).

Required checks (automated in `tests/run_ucsf_velocyto_exact_100k.sh` when fixtures are available):

0. **Independent reference (default):** `UCSF_VELOCYTO_BASELINE_OUTDIR` must point to outputs from an unmodified STAR build / pre-refactor commit (`scripts/save_velocyto_baseline.sh`). Dev-only without baseline: `UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1` (prints WARNING; not publication-grade).
1. run STAR `Velocyto` stream path (default, no `STAR_VELOCYTO_DETERMINISTIC_REPLAY`)
2. run the deterministic sorted-replay path (`STAR_VELOCYTO_DETERMINISTIC_REPLAY=1`) at the same thread count
3. verify native STAR Velocyto `outs/` MEX on both (canonical runner: `--prepare-mex`)
4. verify with `scripts/compare_velocyto_mex.py --mode all`:
   - identical `Solo.out/Gene/raw` feature and barcode axes (solo mode always checks these; Velocyto/raw often has only `*.mtx`)
   - matrix shapes `(n_features, n_barcodes)` on each layer vs those axes
   - identical packaged raw/filtered axes and layers (`features.tsv.gz`, `barcodes.tsv.gz`, per-layer `.mtx.gz`, totals)
   - matrix equality for `spliced`, `unspliced`, `ambiguous`, and total (`matrix.mtx.gz` / reconstructed total in Solo.out)
5. then verify deterministic 1-thread vs N-thread with the same `--mode all` comparison

**Sorted-replay memory (100K):** optional cap via `UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB`; `report_velocyto_sorted_replay_rss.py` reads `Log.out`.

Optional **before** full sample: if a staged 100K tree is available, run
`run_ucsf_velocyto_exact_100k.sh` for fast iteration. If it is **not** available,
skip straight to full-sample gates (below); that is still valid acceptance when
the full corrected surface passes.

### Phase 3. Full-sample exact parity (`run_ucsf_velocyto_exact_2m.sh`)

Run the exact sidecar on the **full-sample** perturb surface (canonical:
corrected **EBs2_2** via `UCSF_2M_*` overrides and `UCSF_2M_PFCONFIG`). The
harness name is `exact_2m`; it is **not** tied to the legacy iPSC2 path.

That script **requires** `UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB` unless
`UCSF_VELOCYTO_ALLOW_UNCAPPED_2M=1` (dev-only); production acceptance must not
use the opt-out.

Required checks:

1. legacy STAR `Velocyto` (stream) on the chosen full-sample inputs
2. deterministic sorted-replay path on the same inputs
3. native raw/filtered Velocyto `outs/` MEX packaging on both
4. `compare_velocyto_mex.py --mode all` plus deterministic 1 vs N threads
5. full exact diff of:
   - features
   - raw barcodes
   - filtered barcodes
   - `spliced.mtx`
   - `unspliced.mtx`
   - `ambiguous.mtx`
   - `matrix.mtx`

Success criterion:

- exact equality, not just correlation

### Phase 4. Real external `velocyto.py` comparison

This is a distinct validation target from legacy STAR parity.

Purpose:

- show how STAR legacy-style `Velocyto` compares with the widely used external
  `velocyto.py` workflow on the same UCSF perturb input surface

This step requires full retained aligned BAMs. Header-only cleanup stubs are not
usable.

If full BAMs are available, run:

1. the retained STAR-aligned BAMs through external `velocyto.py`
2. normalize the output to comparable feature and barcode axes
3. compare against:
   - legacy STAR `Velocyto`
   - new exact sidecar output

Metrics to record:

- raw feature count
- raw barcode count
- filtered barcode count
- total nnz
- layer-wise nnz
- barcode Jaccard on filtered output
- per-layer matrix equality if axes and semantics align exactly
- otherwise per-gene/per-barcode correlations plus disagreement summaries

If full BAMs are not available on a given local UCSF run, document that the
external `velocyto.py` comparison is blocked rather than silently skipping it.

In-repo placeholder: `tests/run_ucsf_velocyto_external_compare.sh` prints **BLOCKED** to stderr and exits 0
(so it is not mistaken for a passing scientific check).

### Phase 5. Hash-backed integration

**Implemented (Stage 2):** `STAR_VELOCYTO_INTEGRATED_HASH=1` with
`STAR_VELOCYTO_DETERMINISTIC_REPLAY=1` runs `countVelocytoSortedReplayCBuckets`.
**Default:** records are **sharded to temporary binary spill files** under
`Solo.out/Velocyto/` (`mkdtemp`), `iCB % STAR_VELOCYTO_INTEGRATED_HASH_SPILL_BUCKETS`
(default **128**, max 4096). After the scan, each shard is loaded, sorted by
`(iCB, UMI, readId)`, merged with `applyVelocytoMerge`, then deleted—peak RAM is
roughly **one shard’s records** plus merge state, not the full dataset at once.
`STAR_VELOCYTO_INTEGRATED_HASH_INMEMORY=1` restores the all-in-RAM per-CB vector
path for debugging only.

`scripts/report_velocyto_sorted_replay_rss.py` logs **both** Stage 1 materialization
and the **post-spill** line `RAM after Velocyto integrated-hash spill …`, and prints
`PER_LOG_MAX_VM_RSS_KB[<run_subdir>]` per `Log.out` so Stage 2 memory can be
compared to Stage 1 even when a single global `MAX_VM_RSS_KB` is dominated by
another run in the same invocation.

This does **not** yet fuse Velocyto recording with the GeneFull bridge hash;
streams are still read from the existing Velocyto temp files.

### Phase 6. Hash-backed UCSF validation

The integrated hash-backed version must be compared against both:

- legacy STAR `Velocyto`
- the exact Stage 1 sidecar implementation

Required checks:

1. exact matrix equality on UCSF 100K
2. exact matrix equality on UCSF 2M
3. wall / RSS / phase timing comparison vs Stage 1
4. confirm that `Gene` / `GeneFull` outputs remain unchanged

**Automated in-repo (partial):**

- `tests/run_ucsf_velocyto_hash_100k.sh` / `run_ucsf_velocyto_hash_2m.sh` run `scripts/compare_velocyto_mex.py --mode genes` between `stream_t1` vs `det_sort_t1`, `stream_t1` vs `det_hash_t1`, and `det_hash_t1` vs `det_hash_tN` (exact `barcodes.tsv`, `features.tsv`, `matrix.mtx` under `Solo.out/Gene` and `Solo.out/GeneFull`, raw + filtered).
- `tests/run_ucsf_velocyto_exact_100k.sh` / `run_ucsf_velocyto_exact_2m.sh` run the same `--mode genes` checks for Stage 1 (`stream_t1` vs `det_t1`, `det_t1` vs `det_tN`).
- Harness `SUMMARY.txt` includes the first four lines of each run’s `Log.final.out` (started / finished timestamps) for wall-clock comparison; Velocyto VmRSS remains via `report_velocyto_sorted_replay_rss.py`.

Full wall-time regression vs Stage 1 still relies on those artifacts plus host notes; Phase 4 (external `velocyto.py`) stays **blocked** until retained BAMs exist.

## Testing Matrix

### A. Semantic correctness

Use exact diff where possible:

- feature table equality
- raw barcode equality
- filtered barcode equality
- `spliced` matrix equality
- `unspliced` matrix equality
- `ambiguous` matrix equality
- `total` matrix equality

**Stream vs deterministic (required, same binary):** same refactored binary and fixture, `--runThreadN 1`, once without `STAR_VELOCYTO_DETERMINISTIC_REPLAY` (stream) and once with it (sorted replay); `compare_velocyto_mex.py --mode all` must pass. This does **not** prove the refactor matches pre-change STAR if both arms share a bug; see independent baseline below.

**Independent baseline (harness default):** `UCSF_VELOCYTO_BASELINE_OUTDIR` must be set to a frozen STAR `outFileNamePrefix` (+ packaged `outs/`). Harness compares it to `stream_t1` before same-binary checks. `UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1` skips baseline (dev only). Canonical STAR runs require a **fresh** `--out-prefix` unless `UCSF_VELOCYTO_REUSE_STAR_OUTDIR=1`.

**Thread determinism (required):** `STAR_VELOCYTO_DETERMINISTIC_REPLAY=1` at 1 vs N threads; `--mode all` must pass. Legacy stream 1 vs N remains optional (`UCSF_VELOCYTO_CHECK_LEGACY_THREAD_PARITY=1`).

### B. Packaging correctness

Reuse the current packaging invariants already exercised in:

- `tests/run_perturb_velocyto_mex_smoke.sh`

Exact-parity harnesses also diff packaged trees via `scripts/compare_velocyto_mex.py --mode all`: per-run `velocyto_feature_bc_matrix_manifest.json` is checked against that run’s gz MEX axes and nnz, then cross-run diffs. `scripts/save_velocyto_baseline.sh` now expects STAR's native Velocyto `outs/` writer so baselines include packaged outs without Python backfill.

Packaging invariants:

- `matrix == spliced + unspliced + ambiguous`
- raw/filtered feature tables match
- filtered barcodes are a subset of raw barcodes
- manifest dimensions match emitted matrices

### C. Performance correctness

For each implementation stage, record:

- wall time
- peak RSS
- key Solo phase timings from `Log.out`
- any dedicated phase timing added for the new `Velocyto` path

Do not benchmark in parallel with other jobs on this host.

## Suggested Script Layout

Keep publication / audit scripts in-repo.

In-repo (initial wiring):

- `scripts/run_star_velocyto_canonical.sh` — pinned STAR args for `100k` / `2m` profiles
- `scripts/save_velocyto_baseline.sh` — canonical STAR + mandatory native Velocyto `outs/` MEX for baseline dirs used with `--mode all`
- `tests/run_ucsf_velocyto_exact_100k.sh` — default: frozen baseline vs stream_t1, then stream vs det, det@1t vs Nt (`UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1` skips baseline)
- `tests/run_ucsf_velocyto_exact_2m.sh` — same; **requires** `UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB` or `UCSF_VELOCYTO_ALLOW_UNCAPPED_2M=1`
- `tests/run_ucsf_velocyto_hash_100k.sh` / `tests/run_ucsf_velocyto_hash_2m.sh` — Stage 2: baseline/stream/stage-1-det vs CB-bucket det + 1 vs N threads (`STAR_VELOCYTO_INTEGRATED_HASH=1`). When `UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB` is set, the **2M (and optional 100K) cap checks `PER_LOG_MAX_VM_RSS_KB[det_hash_t1]` and `[det_hash_tN]` only**, not the global `MAX_VM_RSS_KB` across `det_sort_t1`.
- `scripts/compare_velocyto_mex.py` — exact diff: `Solo.out/Gene/raw` axes + `Solo.out/Velocyto/raw` layers, and/or packaged `outs/`
- `scripts/report_velocyto_sorted_replay_rss.py` — VmRSS after sorted-replay materialization from `Log.out`

Planned extensions:

- `tests/run_ucsf_velocyto_hash_*.sh` — Stage 1 sorted replay vs Stage 2 CB-bucket + thread parity (see script headers)
- `tests/run_ucsf_velocyto_external_compare.sh` — prints **BLOCKED** until full retained BAMs + `velocyto.py` wiring exist (not a passing check)

## Validation run order (implementation is done; acceptance is not)

No further Velocyto bridge **implementation** is required for acceptance;
**execution** of these gates is. Serialize benchmark-style runs (one job at a
time). Use a **fresh** `--baseline-dir` / harness outdir per attempt.

**Prerequisites (all runs):**

- Pre-refactor **baseline** `STAR` binary (separate checkout or worktree at the
  frozen commit); `export STAR_BIN=/path/to/baseline/STAR` when capturing
  baselines.
- Implementation checkout `STAR` for exact/hash harnesses:
  `export STAR_BIN=/path/to/STAR-suite/core/legacy/source/STAR`.
- `export UCSF_VELOCYTO_PARITY_THREADS=8` (or your target thread count).
- After each `save_velocyto_baseline.sh` completes:
  `export UCSF_VELOCYTO_BASELINE_OUTDIR=/path/to/baseline/star_run`.
- **2M RSS budget:** set `export UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB=<kB>`
  before `run_ucsf_velocyto_exact_2m.sh` and `run_ucsf_velocyto_hash_2m.sh`
  (Stage 1: global `MAX_VM_RSS_KB` from deterministic logs; Stage 2: cap
  applies to `PER_LOG_MAX_VM_RSS_KB[det_hash_t1]` and `[det_hash_tN]` only).
  Until this value is chosen (and fixtures exist), acceptance remains open.

### Minimum fixture env (override when defaults in `tests/external_fixtures_env.sh` are wrong)

Set host-specific GEX/guide roots and pfMultiConfig paths; keep genome,
feature ref, and whitelist pinned if these defaults match your site:

```bash
export UCSF_100K_GEX_DIR="..."
export UCSF_100K_GUIDE_DIR="..."
export UCSF_100K_PFCONFIG="..."
export UCSF_100K_FEATURE_REF="/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv"
export UCSF_100K_CB_WHITELIST="/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt"
export UCSF_100K_GENOME_DIR="/storage/autoindex_110_44/bulk_index"

export UCSF_2M_GEX_DIR="..."
export UCSF_2M_GUIDE_DIR="..."
export UCSF_2M_PFCONFIG="..."
export UCSF_2M_FEATURE_REF="/mnt/pikachu/ucsf-perturb-seq/cellranger_feature_ref_hCRISPRa_v2_like_AALG2_pattern.csv"
export UCSF_2M_CB_WHITELIST="/home/lhhung/cellranger-9.0.1/lib/python/cellranger/barcodes/translation/3M-february-2018_NXT.txt"
export UCSF_2M_GENOME_DIR="/storage/autoindex_110_44/bulk_index"
```

`scripts/run_star_velocyto_canonical.sh` resolves `UCSF_2M_FEATURE_REF` and
`UCSF_2M_CB_WHITELIST` the same way as STAR: if unset, it falls back to
`UCSF_100K_FEATURE_REF` / `UCSF_100K_CB_WHITELIST`, then to built-in paths under
`/mnt/pikachu/...`. **Export them explicitly** when your site differs, so you
do not accidentally inherit the wrong ref or whitelist.

### Rerun sequence (flexible; full-sample acceptance is the bar)

**Acceptance** for current UCSF work is **full corrected EBs2_2** (or equivalent
full-sample surface) passing baseline + Stage 1 exact + Stage 2 hash with frozen
baseline STAR. A staged **100K** run is **optional** for speed when that fixture
exists; hosts that only have the full sample should run the **full-sample**
steps directly with explicit `UCSF_2M_*` (and `UCSF_2M_PFCONFIG`) — do **not**
require 100K first.

**A. Optional fast path (only if 100K fixture + baseline binary available)**

1. **100K frozen baseline:**  
   `bash scripts/save_velocyto_baseline.sh --profile 100k --threads 1 --baseline-dir /path/to/baseline_100k`  
   then `export UCSF_VELOCYTO_BASELINE_OUTDIR=/path/to/baseline_100k/star_run`.
2. **Stage 1 exact — 100K:** `bash tests/run_ucsf_velocyto_exact_100k.sh`
3. **Stage 2 hash — 100K:** same `UCSF_VELOCYTO_BASELINE_OUTDIR` as step 1:  
   `bash tests/run_ucsf_velocyto_hash_100k.sh`

**B. Full-sample path (required for acceptance; use EBs2_2 overrides above)**

1. **Full-sample frozen baseline:**  
   `bash scripts/save_velocyto_baseline.sh --profile 2m --threads 1 --baseline-dir /path/to/baseline_full`  
   (set all relevant `UCSF_2M_*` first, including **feature ref** and **CB whitelist**),  
   then `export UCSF_VELOCYTO_BASELINE_OUTDIR=/path/to/baseline_full/star_run`.
2. **Stage 1 exact — full sample:**  
   `bash tests/run_ucsf_velocyto_exact_2m.sh`  
   (**requires** `UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB` or dev opt-out).
3. **Stage 2 hash — full sample:**  
   `bash tests/run_ucsf_velocyto_hash_2m.sh`  
   (same baseline `star_run` and RSS env as step B2).

You may run **A then B** (recommended when both surfaces exist) or **B alone**.

**Later:** external `velocyto.py` — still blocked without retained BAMs
(`tests/run_ucsf_velocyto_external_compare.sh`).

**Immediate next step:** set the **minimum fixture env** (including
`UCSF_2M_FEATURE_REF` / `UCSF_2M_CB_WHITELIST`), choose fresh outdirs, and run
**B**; add **A** only if a 100K tree is on the host.

## Acceptance Criteria

### Stage 1 accepted when

- **Full-sample** (e.g. corrected EBs2_2): frozen baseline vs `stream_t1`
  (`UCSF_VELOCYTO_BASELINE_OUTDIR`, default-on in harness; opt out only with
  `UCSF_VELOCYTO_ALLOW_SAME_BINARY_ONLY=1`). Optional **100K** run is an extra
  confidence check when that fixture exists, not a substitute for full-sample
  gates.
- Refactored stream matches deterministic sorted replay at 1 thread (`compare_velocyto_mex.py --mode all`)
- **Thread determinism:** `STAR_VELOCYTO_DETERMINISTIC_REPLAY=1` runs match at 1 vs N threads with `--mode all`
- **2M memory (Stage 1 exact harness):** `UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB` enforced vs global `MAX_VM_RSS_KB` from deterministic logs (or `UCSF_VELOCYTO_ALLOW_UNCAPPED_2M=1` dev-only)
- packaged raw/filtered MEX invariants pass (`tests/run_perturb_velocyto_mex_smoke.sh` on 100K)

### Stage 2 accepted when

- Stage 2 hash-backed output is exactly equal to Stage 1 and legacy STAR
- **Full-sample** hash harness passes (`run_ucsf_velocyto_hash_2m.sh` with the
  same surface as Stage 1). **100K** hash (`run_ucsf_velocyto_hash_100k.sh`) is
  optional when the staged 100K fixture is unavailable.
- **Thread determinism** for the hash-backed path (1 vs N threads) passes
- normal `Gene` / `GeneFull` outputs remain unchanged
- **2M memory (Stage 2 hash harness):** `UCSF_VELOCYTO_MAX_SORTED_REPLAY_RSS_KB` enforced per `PER_LOG_MAX_VM_RSS_KB[det_hash_t1]` and `[det_hash_tN]` (or `UCSF_VELOCYTO_ALLOW_UNCAPPED_2M=1`); compare Stage 2 per-log RSS to Stage 1 for memory claims

## Failure Modes To Watch

- collapsing to gene too early and losing transcript-intersection semantics
- combining masks with `AND` instead of `OR`
- accidentally treating one UMI as contributing to multiple splice classes
- silently changing filtered-barcode gating away from `Gene` cell filtering
- using final `GeneFull` matrices as a shortcut for splice classification
- mixing optimization work into Stage 1 before exactness is proven

## Bottom Line

The safe order is:

1. exact separate sidecar first
2. **Full-sample** frozen baseline + exact parity on corrected **EBs2_2** (or
   equivalent), with explicit `UCSF_2M_*` including **feature ref** and **CB
   whitelist**, and a chosen full-sample RSS cap
3. **Optional:** if a staged 100K tree exists, repeat baseline + exact (+ hash)
   on 100K for faster iteration — not required for acceptance if step 2 passes
4. Stage 2 CB-bucket deterministic path (`STAR_VELOCYTO_INTEGRATED_HASH=1`),
   validated vs Stage 1 via `run_ucsf_velocyto_hash_2m.sh` on the **same**
   full-sample surface (and optionally `run_ucsf_velocyto_hash_100k.sh`)
5. external `velocyto.py` comparison when full BAMs exist

Further fusion of Velocyto with the GeneFull bridge hot path (single pass, no
Velocyto temp stream) remains future work if profiling warrants it.
