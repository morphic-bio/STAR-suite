# Runbook: reuse the TranscriptVB auto-detect window reads (close the integrated-quant gap)

**For:** a coding agent working in the STAR-suite repo. **Goal:** make
`--quantVBLibType A` (auto-detect) count the library-format detection-window reads,
so STAR Suite's *integrated* quant matches Salmon on small data (and is exact at
scale). This is the fix for the parity gap diagnosed in
[`HANDOFF_SALMON_PARITY.md`](HANDOFF_SALMON_PARITY.md).

This runbook originally specified a narrow detection-window replay fix. That patch
was tested and is **not sufficient**. Treat the four-edit recipe below as historical
context only; the replacement plan is the Salmon-semantics parity checklist in the
new section immediately after the execution note.

> **Execution note (2026-06-25): do not commit the four-edits-only patch as-is.**
> The patch builds, but the chr22 validation did **not** meet the pass criteria.
> Run `runs/chr22_20260625_174250` was generated from a clean rebuild with the four
> edits applied:
>
> - auto `A` still reported integrated tx Pearson `0.9920`, gene Pearson `0.9998`,
>   TranscriptVB total `12,573` vs same-BAM Salmon `12,925` (gap `-352`).
> - The transcriptome BAM was not duplicated; same-BAM Salmon stayed effectively
>   exact (`0.9999976` tx Pearson vs the reference-chain Salmon).
> - STAR Suite auto-detected `OSF` for normal `R1 R2` input
>   (`OSF=138, ISF=7`), while Salmon on the STAR Suite transcriptome BAM detected
>   `ISR`. The pre-fix baseline also detected `OSF`, so the dominant chr22 gap is
>   not explained by the first 1000 reads being skipped.
> - Fixed `IU` preserved total counts (`12,936`, gap `+11` vs same-BAM Salmon) but
>   only reached tx Pearson `0.9929`; fixed `ISF`/`ISR` discarded most counts for
>   this fixture. EM mode improved only slightly (`IU` tx Pearson `0.9938`).
> - Swapping input order to `R2 R1` made auto-detect choose `ISF`
>   (`ISF=138, OSF=7`) but still left a same-BAM gap of about `-148` reads and tx
>   Pearson `0.9924`.
>
> Conclusion: the pure-vote + rewind idea may still be useful cleanup, but it is
> not sufficient to close this harness gap. The next investigation should focus on
> TranscriptVB mate-order/orientation interpretation and same-BAM EC/effective-length
> parity before reviving this patch.

## 0. Updated diagnosis from Salmon 1.10.3 source comparison (2026-06-26)

Salmon is permissively licensed, so the parity target was checked against the
Salmon 1.10.3 source (`COMBINE-lab/salmon` tag `v1.10.3`, commit
`a2f6912b3f9f9af91e3a4b0d74adcb3bdc4c9a32`). The gaps below are what the four-edit
patch missed.

### A. Auto-detection is inline and permissive, not a separate skipped window

In `src/SalmonQuantifyAlignments.cpp`, when the detector is active Salmon sets
`incompatPrior = LOG_1`, samples `aln->libFormat()` for alignments as they stream,
and still quantifies those same alignments. When `detector->canGuess()` is true,
it calls `mostLikelyType()` and restores the normal incompatibility prior. There is
no skipped detection window and no replay/rewind.

In `include/LibraryTypeDetector.hpp`, Salmon samples **50,000 observed alignment
formats** by default, not 1,000 reads. It then chooses orientation and strandedness
from aggregated counts:

- orientation is the maximum of total inward / outward / same;
- strandedness is forward ratio `<0.3` = reverse, `<0.7` = unstranded, otherwise
  forward.

STAR Suite currently runs a separate first-`autoDetectWindow` read pass, votes on
one selected transcript per read, and uses a winner/ambiguity policy with a special
inward-collapse rule. This explains why the chr22 run saw `OSF=138, ISF=7` in STAR
Suite while Salmon detected `ISR` from the same transcriptome BAM.

### B. Paired-end orientation needs Salmon's dovetail-aware hit type

Salmon has two `hitType()` overloads in `src/SalmonUtils.cpp`. The length-aware
paired-end overload accepts read lengths and `canDovetail`; for opposite-strand
mates it allows a mate to "stretch" by the other read length before declaring the
pair outward. That can convert overlap/dovetail cases that the simple position-only
test calls `OSF/OSR` into inward `ISF/ISR`.

STAR Suite currently uses only the simple position/orientation version in
`LibFormatDetection.cpp` and `ec_builder.cpp`. This is a likely cause of the
spurious outward votes and mate-order sensitivity.

### C. Rich EC weights are accumulated but never finish-normalized

Salmon's `EquivalenceClassBuilder<TGValue>::addGroup()` increments the EC count and
adds per-alignment weights; `finish()` then calls `TGValue::normalizeAux()` so each
EC's accumulated auxiliary weights sum to 1 before EM/VB.

STAR Suite's inline path in
`core/features/vbem/source/TranscriptQuantEC.cpp::addReadAlignments()` does the same
count increment and weight accumulation, but `TranscriptQuantEC::finalize()` is
currently a no-op. EM/VB then sees both `ec.count` and unnormalized accumulated
weights. The first low-risk correction is to normalize each EC's weights in
`finalize()` after merging all thread-local ECs.

### D. Range-factorized EC keying is not actually implemented in the inline path

The same-BAM Salmon output reports:

```json
"eq_class_properties": ["range_factorized", "gzipped"]
```

Salmon's default alignment-mode run uses 4 range-factorization bins, and `--useErrorModel`
forces at least 4 bins. The important implementation detail is that Salmon appends
range-bin IDs to the **EC key** (`TranscriptGroup::txps`) but keeps the real
transcript IDs and weight vector length separate (`getNumTranscriptsForClass()`
returns `weights.size()`). The optimizer consumes only the first `N` real transcript
IDs and their `N` weights.

STAR Suite's `getParams()` sets `use_range_factorization = true`, but the live inline
path only applies `applyRangeFactorization()` in trace output. It cannot simply append
the pseudo-bin IDs to `EC::transcript_ids`, because the current EM/VB representation
would treat them as real transcript indices. The fix needs a separate EC signature
key that may include range bins while the stored `EC` keeps only real transcript IDs
and weights.

### E. Fragment length / effective length semantics still differ

Salmon's paired-end `ReadPair::fragLengthPedantic()` uses BAM read positions and
sequence lengths, clamps to transcript length, and returns `abs(p1 - p2)` for
opposite-strand pairs. STAR Suite currently computes fragment length from inferred
5-prime/3-prime endpoints with a `+1` span. In the chr22 same-BAM run Salmon's
fragment-length mean was `303.25`, while STAR Suite fixed-`IU` reported `202.47`;
this effective-length difference can distort TPM and some count allocation even when
total assigned fragments are close.

## 0a. Replacement implementation order

1. Fix same-BAM EC semantics first:
   normalize accumulated EC weights in `TranscriptQuantEC::finalize()`, then rerun
   fixed-`IU` and same-BAM Salmon comparisons.
2. Add a separate range-factorized EC signature:
   build a key from real transcript IDs plus Salmon-style range bins, but store only
   real transcript IDs/weights in `EC`. Keep the optimizer interface unchanged.
3. Replace the TranscriptVB auto-detect model with Salmon's inline permissive detector:
   sample observed alignment formats during the main pass with incompatibility
   disabled until the detector finalizes, then continue without rewinding.
4. Port the dovetail-aware paired-end `hitType()` and use it consistently for
   detection and compatibility filtering.
5. Reconcile `fragLengthPedantic()` and effective-length calculation against Salmon's
   alignment-mode path.
6. Rerun the chr22 harness. Passing now means:
   same-BAM Salmon parity remains ~1.0, integrated auto `A` no longer drops the
   detection-window reads or mis-detects as `OSF`, and total TranscriptVB reads match
   Salmon within the small residual expected from alignment differences.

### 0b. Implementation note from 2026-06-26

The Salmon-semantics path above has now been implemented in the working tree,
including EC weight normalization, range-factorized EC signatures, fixed
paired-end orientation/fragment-length handling, strict fixed-libtype controls,
and the auto-detect vote/rewind path. One additional parity bug was found during
validation: STAR Suite was treating pedantic fragment lengths over `1000` as
invalid. Salmon keeps the true fragment span for start-position probability and
only clamps inside the FLD histogram lookup/update, so long fragments must remain
positive in `RawAlignment::fragment_len`.

Validation on the cached chr22 fixture:

| Comparison | NumReads Pearson | TPM Pearson | total NumReads | half-L1 moved |
|---|---:|---:|---:|---:|
| STAR fixed `ISR`, no error, no FLD/no effective correction vs Salmon same flags | 0.99999999999 | 0.999999999998 | 12,771.996 vs 12,771.998 | 0.010 |
| STAR fixed `ISR`, no error vs Salmon fixed `ISR`, no error, `-p 1` | 0.999995 | 0.999997 | 12,772.004 vs 12,772.005 | 9.351 |
| STAR auto `A`, no error vs Salmon fixed `ISR`, no error, `-p 1` | 0.999991 | 0.999991 | 12,771.994 vs 12,771.999 | 13.034 |

Salmon's alignment-mode online FLD is thread-order sensitive on this small
fixture: Salmon `-p 1` vs Salmon `-p 32` gives only about `0.9957` NumReads
Pearson and reproduces the same large outlier pattern. This is because Salmon
updates its alignment-mode FLD online as alignment records arrive; changing worker
thread count changes the early observation order, which then changes the FLD and
effective-length weights for a small number of multimappers. STAR Suite's
integrated quantification also builds an online FLD from threaded alignment chunks,
so changing STAR `--runThreadN` can move a small number of fractional
read-equivalents on this tiny fixture. For the chr22 smoke, the stable semantic
target is therefore STAR at 32 threads, matching the retained artifact, versus
Salmon `-p 1`; `run_compare.sh` now separates STAR threads from Salmon quant
threads and defaults Salmon to one thread only in `--mode chr22`.

The retained local artifact is
`benchmarks/nfcore_rnaseq_compare/runs/chr22_20260626_003522/`. Keep it for
diagnostic provenance, but use `tests/run_transcriptvb_chr22_parity_smoke.sh` as
the integration/smoke test surface. The smoke gates on tolerances rather than
exact last-digit equality because STAR's threaded online FLD can move a few
fractional read-equivalents between 32-thread reruns on this tiny fixture.

---

## 1. Historical background (read once)

With `--quantVBLibType A`, STAR runs a **separate detection pass** over the first
`quantVBAutoDetectWindow` reads (default **1000**, `core/legacy/source/Parameters.h`
`autoDetectWindow = 1000`) to vote on the library format
(`core/legacy/source/STAR.cpp`, the `if (P.quant.transcriptVB.yes && libType=="A")`
block, ~line 1817). Two problems in the current code:

1. In `core/features/vbem/source/TranscriptQuantEC.cpp::addReadAlignments`, the EC
   **commit is NOT gated by `inDetectionMode`** — so detection-window reads commit to
   the EC table with the *un-gated / not-yet-determined* format.
2. After detection, the read stream is **not rewound** (the in-code comment confirms
   "…will run again … with no remaining reads"), so the main pass continues *past*
   the detection window.

Net: the first ~1000 reads are mis/under-counted. On the small chr22 fixture this is
visible (integrated transcript NumReads Pearson **0.9922**, total assigned reads
**12,580 vs Salmon 12,924**, a −344 gap). At production scale it is negligible
(<0.1% at 10M+ reads — see `core/features/vbem/docs/STAR_Salmon_Quantification_Parity_Report.md`),
and `--quantVBLibType ISF` already gives a zero gap. We are fixing the auto path so
it counts those reads too.

**Approach — mirror the existing SLAM detection pass** (this is why SLAM was the
reference): make the TranscriptVB detection pass a **pure-vote, no-output** pass, then
**rewind the read files** so the main pass re-reads ALL reads (including the window)
through the *existing, validated* build→commit→error-model/FLD path with the detected
format. Do **not** refactor `addReadAlignments`'s commit/update logic — the
error-model/FLD/burn-in is stateful and order-dependent, so replaying through the real
path is the low-risk choice. Re-aligning ~1000 reads is negligible.

**Canonical source (AGENTS.md):** `core/legacy/` is "the single source of truth for
the STAR core." Edit `core/legacy/source/{STAR.cpp,BAMoutput.cpp,ReadAlign_outputTranscriptSAM.cpp}`
and `core/features/vbem/source/TranscriptQuantEC.cpp` (compiled via the Makefile
`VPATH`). **Do not** edit the duplicate `STAR.cpp`/etc. copies under
`core/features/vbem/source/` — `legacy` wins in VPATH and is what builds
`core/legacy/source/STAR` → `/usr/local/bin/STAR`.

---

## 2. Historical branch

```bash
cd /mnt/pikachu/STAR-suite
git checkout -b fix/transcriptvb-detection-window-reuse   # (or reuse an existing fix/ branch)
```

## 3. Historical four-edit patch (do not apply as-is)

### Edit A — `core/features/vbem/source/TranscriptQuantEC.cpp`
Make the detection pass vote-only: early-return right after the vote block, before
`getParams()`. Anchor on:
```cpp
    // Get params (gating OFF during detection, ON after)
    ECBuilderParams params = getParams();
```
Insert immediately **above** that:
```cpp
    // Detection pass is pure observation: vote only (above), then stop. The first
    // autoDetectWindow reads are re-processed by the main mapping pass (after the
    // read-file rewind in STAR.cpp) with the detected library format, so each is
    // committed exactly once with the correct format -- instead of being committed
    // here under an undetermined format.
    if (P_.quant.transcriptVB.inDetectionMode) {
        return;
    }

```

### Edit B — `core/legacy/source/BAMoutput.cpp` (TWO identical sites)
Suppress BAM output during the detection pass. Both occurrences read:
```cpp
    if (P.quant.slam.autoTrimDetectionPass) return;
```
Change **both** to:
```cpp
    if (P.quant.slam.autoTrimDetectionPass || P.quant.transcriptVB.inDetectionMode) return;
```

### Edit C — `core/legacy/source/ReadAlign_outputTranscriptSAM.cpp`
Suppress transcriptome SAM/BAM output during detection. Change:
```cpp
    if (P.quant.slam.autoTrimDetectionPass) return 0;
```
to:
```cpp
    if (P.quant.slam.autoTrimDetectionPass || P.quant.transcriptVB.inDetectionMode) return 0;
```

### Edit D — `core/legacy/source/STAR.cpp`
Rewind the read files after the TranscriptVB detection pass (mirrors the SLAM rewind
later in the same file). Anchor on the end of the detection block:
```cpp
        extern std::atomic<uint64_t> global_processed_fragments;
        global_processed_fragments.store(0, std::memory_order_relaxed);
        Parameters::global_fld_obs_count.store(0, std::memory_order_relaxed);
    }
```
Insert the rewind **between** the `global_fld_obs_count.store(...)` line and the
closing `}`:
```cpp
        Parameters::global_fld_obs_count.store(0, std::memory_order_relaxed);

        // REWIND the read files so the main mapping pass re-reads ALL reads
        // (including the detection window) and quantifies them with the now-detected
        // library format. The detection pass votes only (output suppressed, no EC
        // commit), so each read is processed/output/counted exactly once by the main
        // pass. Mirrors the SLAM detection-pass rewind below.
        P.closeReadsFiles();
        P.iReadAll = 0;                       // Reset read counter
        g_threadChunks.chunkInN = 0;          // Reset chunk counter
        g_threadChunks.chunkOutN = 0;         // Reset output chunk counter
        P.readFilesIndex = 0;                 // Reset file index
        g_bamRecordIndex.store(0);            // Reset BAM record index
        P.openReadsFiles();
        // Reset global mapping stats so Log.final.out reflects only the main pass
        g_statsAll.resetN();
        time(&g_statsAll.timeStartMap);
    }
```
(`g_threadChunks`, `g_bamRecordIndex`, `g_statsAll` are already used by the SLAM block
later in this file, so they are in scope — copy that block's symbols if unsure.)

## 4. Build

```bash
make -C core/legacy/source -j"$(nproc)" STAR WITH_CHROMAP=1
```
(The `core` target only exists at the repo root; from the subdir the *default* goal is
a phony `FORCE` — you must name the `STAR` target. Binary: `core/legacy/source/STAR`.)
Confirm: the binary mtime updates and `./core/legacy/source/STAR --version` prints `1.4.1`.

## 5. Test 1 — chr22 parity (the fix target)

The harness defaults `STAR_SUITE_BIN` to `core/legacy/source/STAR`, so it picks up the
new binary automatically. The chr22 read fixture is cached under `fixtures/chr22/`, so
this is fast.
```bash
cd benchmarks/nfcore_rnaseq_compare
STAR_UPSTREAM_BIN=/mnt/pikachu/bin/bin/Linux_x86_64_static/STAR \
  ./run_compare.sh --mode chr22 \
  --salmon-threads 1 \
  --reads-dir runs/chr22_20260625_150321/reads     # any existing reads dir; fixture is cached
cat runs/chr22_*/compare/report.md | sed -n '1,20p'
```
**Pass criteria (vs the pre-fix baseline run `runs/chr22_20260625_154221`):**
| Metric | pre-fix | expected post-fix |
|---|---|---|
| integrated tx Pearson (NumReads) | 0.9922 | **≥ ~0.999** |
| integrated gene Pearson | 0.9999 | ≥ 0.9999 (no regression) |
| total assigned reads: TranscriptVB vs Salmon | 12,580 vs 12,924 | **gap → ~0** against Salmon `-p 1` |
| **exact mode (BAM→Salmon) tx Pearson** | 1.0000 | **must stay 1.0000** (regression guard) |

Also sanity-check **no duplicate output**: `samtools view -c
runs/chr22_*/B_starsuite/Aligned.toTranscriptome.out.bam` should be a single-pass
count (not ~2× of the pre-fix value). If the transcriptome BAM is doubled, output
suppression is incomplete — see §6.

## 6. Test 2 — SLAM regression + fixed-libtype (the shared-code guard)

The output-gate edits add `|| P.quant.transcriptVB.inDetectionMode`, which is **false**
during a SLAM run, and the `addReadAlignments` early-return only fires when TranscriptVB
is in detection mode — so SLAM behavior must be byte-identical. Confirm:
```bash
cd /mnt/pikachu/STAR-suite
bash tests/run_slam_parity_smoke.sh        # and/or tests/run_slam_end_to_end.sh
```
Expected: passes unchanged. Also confirm the **fixed-libtype** path is untouched (the
detection block is gated on `libType=="A"`, so with `--quantVBLibType ISF` no detection
pass runs, no rewind): run any existing TranscriptVB fixed-libtype test, or
`tests/run_public_bulk_pe_smoke.sh` (uses an explicit format) and confirm unchanged.

## 7. Done criteria

- chr22: integrated tx Pearson ≥ ~0.999, total-reads gap → ~0, gene ≥ 0.9999, **exact
  mode still 1.0**, transcriptome BAM not duplicated.
- SLAM tests green; fixed-libtype path unchanged.
- Then run `./run_compare.sh --mode full` once to confirm at production scale (and that
  the extra ~1000-read re-alignment doesn't dent wall time meaningfully).
- Commit on the branch; update `HANDOFF_SALMON_PARITY.md` to say the gap is fixed (and
  by what), not a "known tradeoff."

## 8. Risks / caveats (watch for these)

- **RAchunk[0] reuse.** The detection pass reuses `RAchunk[0]`; the rewind resets the
  global read counters + reopens the files. If `RAchunk[0]`'s internal state isn't
  fully reset by that (symptoms: crash, wrong read count, or duplicated reads), switch
  to the **more faithful SLAM mirror**: run detection on a *dedicated throwaway*
  `ReadAlignChunk` (like SLAM's `RAdetect` at `STAR.cpp` ~line 1226) so the main
  `RAchunk[]`s stay pristine, then rewind. This is the fallback if the simple rewind
  misbehaves.
- **Output completeness.** The three gates cover genome BAM + transcriptome BAM (the
  bulk outputs). `g_statsAll.resetN()` keeps `Log.final.out` clean. `SJ.out.tab` may
  double-count splice junctions from the detection window (the genome SJ collector
  isn't gated) — acceptable for quant parity, but note it; if exact SJ output matters,
  gate SJ collection on `inDetectionMode` too, or use the dedicated-chunk variant
  (which produces no SJ during detection). For Solo/scRNA runs the same `inDetectionMode`
  gate would be needed on Solo output paths — **out of scope** for this bulk fix.
- **Flag visibility.** `P.quant.transcriptVB.inDetectionMode` is accessible in
  `BAMoutput.cpp` / `ReadAlign_outputTranscriptSAM.cpp` (both take `Parameters& P`).

## 9. References

- Diagnosis + decision history: [`HANDOFF_SALMON_PARITY.md`](HANDOFF_SALMON_PARITY.md).
- Prior parity analysis (fixed-libtype = zero gap; negligible at scale):
  `core/features/vbem/docs/STAR_Salmon_Quantification_Parity_Report.md`.
- Validation harness: `run_compare.sh` (this dir); pre-fix baseline run
  `runs/chr22_20260625_154221`; the same-BAM diagnostic
  `runs/chr22_20260625_154221/compare/transcriptvb_vs_salmon_sameBAM_topdiff.tsv`.
- SLAM detection-pass pattern to mirror: `core/legacy/source/STAR.cpp` SLAM block
  (~lines 1162–1465).
