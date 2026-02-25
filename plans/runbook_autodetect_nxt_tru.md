# Runbook: Auto-Detect NXT/TRU Chemistry Per Library

**Branch:** `autodetect-NXT-TRU`  
**Prereq:** Merge-boundary normalization (commit `02f43e5`) is in place.  
**Goal:** During early read processing in `assignBarcodes`, trial-match cell
barcodes both raw and NXT-translated against the whitelist, report the match
mode, and let the STAR layer compose `effectiveChem` from match mode + known
whitelist namespace.

**Key design constraints:**

1. **Detection-only in process_features.** `translate_NXT` is never set by
   auto-detection. All actual NXT→TRU translation happens at the merge boundary
   in `PfMultiProcess`. This avoids double-translation.

2. **Match mode, not chemistry.** process_features reports whether barcodes
   matched the whitelist as-is (`RAW_MATCH`) or after NXT translation
   (`TRANSLATED_MATCH`). It does not interpret this as a chemistry label —
   that requires knowing the whitelist namespace, which only the STAR layer has.

3. **Decoupled from feature-offset bootstrap.** Chemistry detection has its
   own config knobs (`autodetect_chemistry` + `autodetect_chemistry_reads`)
   and runs in the consumer loop unconditionally, independent of feature
   offsets, anchors, or `feature_mode_bootstrap_reads`.

---

## Architecture Overview

```
PfMultiProcess (per library)
  │
  ├─ assignOpts.autodetectChemistry = true
  │   (skipped when --crChemistry is explicitly NXT or TRU)
  │
  └─ PfMultiAssign::runAssignBarcodes(whitelist, ref, fastq, out, opts)
       │                                              returns AssignResult
       │
       └─ pf_api: pf_config_set_autodetect_chemistry(cfg, 1)
            │      pf_config_set_autodetect_chemistry_reads(cfg, N)
            │      pf_config_set_autodetect_chemistry_min_hits(cfg, M)
            │
            └─ assignBarcodes consumer loop (line ~4165):
                 │
                 ├─ for first N reads (deterministic atomic ticket):
                 │    extract CB from barcode_lines[0] + barcode_constant_offset
                 │    raw_hit   += whitelist_hash lookup(string2code(CB))
                 │    nxt_hit   += whitelist_hash lookup(string2code(translate(CB)))
                 │
                 ├─ when ticket == N (one thread):  ← inline decision
                 │    chem_detect_decide(ctx)
                 │
                 └─ after all consumer threads finish:  ← drain finalization
                      chem_detect_decide(ctx)  // CAS-guarded, no-op if decided
       │
       └─ pf_get_detected_match_mode(ctx) → "RAW_MATCH" | "TRANSLATED_MATCH"
       │                                     | "AMBIGUOUS" | "UNKNOWN"
       │
       PfMultiProcess:
         wlNamespace = detectChemistryFromWhitelistPath()  // already known
         matchMode   = result.detectedMatchMode

         if matchMode == RAW_MATCH:
           effectiveChem = wlNamespace          // reads match whitelist as-is
         elif matchMode == TRANSLATED_MATCH:
           effectiveChem = opposite(wlNamespace) // reads are in opposite namespace
         else:
           effectiveChem = wlNamespace           // fallback to inferred
```

**Why match mode instead of chemistry:**

| Whitelist namespace | RAW_MATCH means | TRANSLATED_MATCH means |
|---------------------|-----------------|----------------------|
| TRU (1-col TRU)     | reads are TRU   | reads are NXT        |
| NXT (1-col from 2-col col1) | reads are NXT | reads are TRU |

Labeling RAW_MATCH as "TRU" is only correct for TRU whitelists. For NXT
whitelists (extracted col 1 of a 2-column file), RAW_MATCH means reads are
NXT. Returning the match mode and letting the STAR layer compose chemistry
eliminates the ambiguity.

---

## Whitelist Contract

process_features always receives a **1-column whitelist** (pre-normalized by
`PfMultiAssign::normalizeWhitelistIfNeeded`). The namespace of that single
column depends on the original whitelist format:

- **1-column TRU file:** passed through as-is. Whitelist is TRU.
- **2-column NXT file:** column 1 extracted. Whitelist is NXT (read-space).

The STAR layer knows which case applies via `detectChemistryFromWhitelistPath`
(already computed before `assignBarcodes` is called). This is a **hard
requirement** for correct chemistry composition. If the whitelist namespace
is unknown, auto-detection should not be enabled.

When `detectChemistryFromWhitelistPath` returns a result, it is logged. If
auto-detection is enabled and the whitelist namespace is not "NXT" or "TRU"
(e.g., the function fell through to default), the STAR layer should log a
warning and treat match mode as informational only.

---

## Stage 0 — Synthetic Test Fixtures

**What:** Generate minimal synthetic data for fast iteration. Tests exercise
the standalone `assignBarcodes` path (1-column pre-normalized whitelist).

**Files to create:**
- `tests/autodetect_nxt_tru/generate_synthetic.py` — generator
- `tests/autodetect_nxt_tru/test_autodetect.sh` — runner

**Synthetic data design:**
- 50 barcodes (16-mer) in TRU namespace → 1-column TRU whitelist
- Same 50 barcodes translated → 1-column NXT whitelist (for NXT-whitelist tests)
- Feature ref: 2 CRISPR guides, 20-mer, known offset
- 4 FASTQ sets (1000 reads each, R1=barcode+UMI, R2=feature):
  1. `reads_tru/` — barcodes in TRU
  2. `reads_nxt/` — barcodes in NXT (positions 7-8 complemented)
  3. `reads_mixed/` — 50/50 mix
  4. `reads_tiny/` — 5 reads only (tests drain finalization)

**Test matrix (standalone):**

| Case | Whitelist | Reads | Expected match mode |
|------|-----------|-------|-------------------|
| 1 | 1-col TRU | TRU | RAW_MATCH |
| 2 | 1-col TRU | NXT | TRANSLATED_MATCH |
| 3 | 1-col NXT | NXT | RAW_MATCH |
| 4 | 1-col NXT | TRU | TRANSLATED_MATCH |
| 5 | 1-col TRU | mixed | AMBIGUOUS |
| 6 | 1-col TRU | tiny | decision made on drain (5 reads) |

**Test validation:**
- Parse log output for `NOTICE: chemistry auto-detect` line
- Verify match mode string matches expected value

**Note:** These tests call `assignBarcodes` standalone. The whitelist is
1-column (pre-normalized). STAR-wrapper tests (Stage 5) cover 2-column input.

**Estimated effort:** ~1 hour

---

## Stage 1 — Detection State in `pf_context` and `pf_api`

**What:** Add detection counters and result to `pf_context` (not globals).
Add config knobs and a context-based getter. No behavior change yet.

**Files:**
- `core/features/process_features/src/pf_api.c`
- `core/features/process_features/include/pf_api.h`

**Changes to `pf_config`:**

```c
struct pf_config {
    // ... existing fields ...
    int autodetect_chemistry;           // 0=off, 1=on
    int autodetect_chemistry_reads;     // N reads to sample (default 10000)
    int autodetect_chemistry_min_hits;  // minimum total hits for decision (default 50)
};
```

Defaults: `autodetect_chemistry = 0`, `autodetect_chemistry_reads = 10000`,
`autodetect_chemistry_min_hits = 50`.

**Changes to `pf_context`:**

```c
struct pf_context {
    // ... existing fields ...
    unsigned long long chem_detect_raw_hits;
    unsigned long long chem_detect_nxt_hits;
    unsigned long long chem_detect_ticket;      // atomic ticket counter
    int chem_detect_done;                       // 0=sampling, 1=decided
    int detected_match_mode;                    // 0=unknown, 1=RAW_MATCH,
                                                // 2=TRANSLATED_MATCH, 3=AMBIGUOUS
    int autodetect_chemistry;                   // copied from config
    int autodetect_chemistry_reads;             // copied from config
    int autodetect_chemistry_min_hits;          // copied from config
};
```

**New `pf_api.h` functions:**

```c
void pf_config_set_autodetect_chemistry(pf_config *config, int enabled);
void pf_config_set_autodetect_chemistry_reads(pf_config *config, int n_reads);
void pf_config_set_autodetect_chemistry_min_hits(pf_config *config, int min_hits);
const char* pf_get_detected_match_mode(pf_context *ctx);
```

**In `pf_init()`:** Copy config values to context, zero counters.

**Getter implementation:**

```c
const char* pf_get_detected_match_mode(pf_context *ctx) {
    if (!ctx) return "UNKNOWN";
    switch (ctx->detected_match_mode) {
        case 1: return "RAW_MATCH";
        case 2: return "TRANSLATED_MATCH";
        case 3: return "AMBIGUOUS";
        default: return "UNKNOWN";
    }
}
```

**Acceptance:** Builds cleanly. No behavior change. A375 smoke passes.

**Estimated effort:** 30 minutes

---

## Stage 2 — Trial-Matching in the Consumer Loop

**What:** In the consumer loop where CB is available, sample the first N reads
using a deterministic atomic ticket counter. For each sampled read, do a
whitelist hash lookup both raw and NXT-translated. Increment context counters.
After all threads finish, drain-finalize if N was not reached.

**File:** `core/features/process_features/src/assignBarcodes.c`

**Insertion point:** In the consumer loop (function around line 4060), right
after `barcode_lines` is copied from the ring buffer (line 4131) and before
the feature-matching call (line 4172). This is where the raw cell barcode is
available in `barcode_lines[0] + barcode_constant_offset`.

**Inline sampling logic:**

```c
// Fast-path: skip entirely once decided (avoids unnecessary ticket increments).
if (ctx->autodetect_chemistry && !ctx->chem_detect_done) {
    unsigned long long ticket = __sync_add_and_fetch(&ctx->chem_detect_ticket, 1);
    if (ticket <= (unsigned long long)ctx->autodetect_chemistry_reads) {
        char *barcode_seq = barcode_lines[0] + barcode_constant_offset;
        if (strlen(barcode_seq) >= (size_t)barcode_length) {
            char bc_buf[barcode_length + 1];
            memcpy(bc_buf, barcode_seq, barcode_length);
            bc_buf[barcode_length] = '\0';

            // Raw lookup
            unsigned char raw_code[barcode_code_length];
            string2code(bc_buf, barcode_length, raw_code);
            uint32_t raw_key = *(uint32_t*)raw_code;
            if (kh_get(u32ptr, whitelist_hash, raw_key) != kh_end(whitelist_hash)) {
                __sync_add_and_fetch(&ctx->chem_detect_raw_hits, 1);
            }

            // NXT-translated lookup
            translate_nxt_inplace(bc_buf, barcode_length);
            unsigned char nxt_code[barcode_code_length];
            string2code(bc_buf, barcode_length, nxt_code);
            uint32_t nxt_key = *(uint32_t*)nxt_code;
            if (kh_get(u32ptr, whitelist_hash, nxt_key) != kh_end(whitelist_hash)) {
                __sync_add_and_fetch(&ctx->chem_detect_nxt_hits, 1);
            }
        }

        // Inline decision at ticket N (CAS-guarded inside chem_detect_decide)
        if (ticket == (unsigned long long)ctx->autodetect_chemistry_reads) {
            chem_detect_decide(ctx);
        }
    }
}
```

**Drain finalization:** After the consumer loop exits (all sets done), before
the thread returns:

```c
// After while(!done) loop exits:
if (ctx->autodetect_chemistry && !ctx->chem_detect_done) {
    chem_detect_decide(ctx); // CAS-guarded inside, safe from any thread
}
```

This handles libraries with fewer than N reads. `chem_detect_decide` uses an
internal CAS to ensure exactly one thread runs the decision; redundant calls
are no-ops.

**`chem_detect_decide()` (static helper):**

```c
static void chem_detect_decide(pf_context *ctx) {
    // CAS guard: only the first caller runs the decision; all others are no-ops.
    if (!__sync_bool_compare_and_swap(&ctx->chem_detect_done, 0, 1)) {
        return;
    }

    unsigned long long raw = ctx->chem_detect_raw_hits;
    unsigned long long nxt = ctx->chem_detect_nxt_hits;
    unsigned long long total = raw + nxt;
    int min_hits = ctx->autodetect_chemistry_min_hits;
    int n_reads = (int)ctx->chem_detect_ticket;

    fprintf(stderr, "NOTICE: chemistry auto-detect: raw_hits=%llu "
            "nxt_hits=%llu total=%llu reads_sampled=%d\n",
            raw, nxt, total, n_reads);

    if (total < (unsigned long long)min_hits) {
        ctx->detected_match_mode = 3; // AMBIGUOUS
        fprintf(stderr, "WARNING: chemistry auto-detect: total hits %llu "
                "below minimum %d, result AMBIGUOUS\n", total, min_hits);
        __sync_synchronize();
        return;
    }

    double raw_frac = (double)raw / (double)total;
    fprintf(stderr, "NOTICE: chemistry auto-detect: raw_frac=%.3f\n", raw_frac);

    if (raw_frac >= 0.8) {
        ctx->detected_match_mode = 1; // RAW_MATCH
    } else if (raw_frac <= 0.2) {
        ctx->detected_match_mode = 2; // TRANSLATED_MATCH
    } else {
        ctx->detected_match_mode = 3; // AMBIGUOUS
        fprintf(stderr, "WARNING: chemistry auto-detect: ambiguous "
                "(raw_frac=%.3f), result AMBIGUOUS\n", raw_frac);
    }
    __sync_synchronize();
}
```

**Minimum evidence gate:** `total < min_hits` → AMBIGUOUS. Default `min_hits`
is 50 (configurable via `pf_config_set_autodetect_chemistry_min_hits`). This
catches degenerate cases (empty FASTQs, all-N barcodes, wrong barcode length).

**Passing `ctx` to the consumer loop:** `sample_args` needs a `pf_context*`
field. If not already present, add it and set it during the processing setup
in `process_files_in_sample`.

**Multi-sample scoping:** In the pf-multi caller path, `PfMultiAssign` creates
a fresh `pf_context` per library call, so counters are naturally per-library.
In standalone `pf_process_fastq_dir` with multiple sample subdirectories,
counters aggregate across all samples within the same run. This is correct
because chemistry should be uniform across samples in the same library — if
it isn't, the ratio will land in the ambiguous range and fall back to inferred.

**Acceptance:** With synthetic TRU whitelist + TRU reads → RAW_MATCH. With
TRU whitelist + NXT reads → TRANSLATED_MATCH. With tiny reads → drain path
triggers, result depends on content.

**Estimated effort:** 1.5 hours (including `sample_args` plumbing)

---

## Stage 3 — PfMultiAssign Wiring

**What:** Add `autodetectChemistry` to `AssignOptions`, change return type to
`AssignResult` struct, wire through `pf_api`, read back detected match mode.

**Files:**
- `core/legacy/source/PfMultiAssign.h`
- `core/legacy/source/PfMultiAssign.cpp`

**Changes to `PfMultiAssign.h`:**

```cpp
struct AssignResult {
    int exitCode = 0;
    string detectedMatchMode; // "RAW_MATCH", "TRANSLATED_MATCH", "AMBIGUOUS", "UNKNOWN"
};

struct AssignOptions {
    // ... existing fields ...
    bool autodetectChemistry = false;
    int autodetectChemistryReads = 10000;
    int autodetectChemistryMinHits = 50;
};

AssignResult runAssignBarcodes(const string& whitelist,
    const string& featureRef, const string& fastqDir,
    const string& assignOut, const AssignOptions& options = AssignOptions());
```

**Changes to `PfMultiAssign.cpp`:**

```cpp
// In applyAssignOptions:
if (options.autodetectChemistry) {
    pf_config_set_autodetect_chemistry(cfg, 1);
    pf_config_set_autodetect_chemistry_reads(cfg, options.autodetectChemistryReads);
    pf_config_set_autodetect_chemistry_min_hits(cfg, options.autodetectChemistryMinHits);
}

// After processing completes:
AssignResult result;
result.exitCode = /* existing return code */;
result.detectedMatchMode = pf_get_detected_match_mode(ctx);
return result;
```

**All existing callers** that check the `int` return value now check
`result.exitCode`. Minimal change.

**Acceptance:** Builds. Existing callers updated. No behavior change when
`autodetectChemistry=false` (detectedMatchMode returns "UNKNOWN").

**Estimated effort:** 45 minutes

---

## Stage 4 — PfMultiProcess Integration

**What:** For each feature library, enable auto-detection, read back match mode,
compose `effectiveChem` from match mode + known whitelist namespace, and set it
on the `FeatureRun`. The existing merge-boundary normalization uses
`effectiveChem` — no additional changes needed downstream.

**File:** `core/legacy/source/PfMultiProcess.cpp`

**Where:** In the per-library loop, when building `runAssignOpts` and after
`runAssignBarcodes` returns.

**Chemistry composition logic:**

```cpp
static string oppositeNamespace(const string& ns) {
    return (upperCopy(ns) == "NXT") ? "TRU" : "NXT";
}

// Before calling runAssignBarcodes:
bool useAutodetect = (context.requestedChem == "auto");
string wlNamespace = prepared.inferredChem; // from detectChemistryFromWhitelistPath
bool wlNamespaceConfident = prepared.inferredChemConfident;

// Always enable detection when auto is requested — results are always logged.
if (useAutodetect) {
    runAssignOpts.autodetectChemistry = true;
}

// After:
auto result = PfMultiAssign::runAssignBarcodes(..., runAssignOpts);
if (result.exitCode != 0) { /* existing error handling */ }

// Compose effectiveChem from match mode + whitelist namespace.
// Only override effectiveChem when whitelist namespace is confidently known.
// When confidence is false, detection result is logged but effectiveChem is
// NOT changed — all libraries stay consistent with GEX to prevent merge
// mismatches (since we have no absolute anchor for the namespace).
if (useAutodetect && wlNamespaceConfident) {
    if (result.detectedMatchMode == "RAW_MATCH") {
        run.effectiveChem = wlNamespace;
        P.inOut->logMain << "NOTICE: auto-detect RAW_MATCH for "
                         << run.featureType << " → effectiveChem="
                         << run.effectiveChem << "\n";
    } else if (result.detectedMatchMode == "TRANSLATED_MATCH") {
        run.effectiveChem = oppositeNamespace(wlNamespace);
        P.inOut->logMain << "NOTICE: auto-detect TRANSLATED_MATCH for "
                         << run.featureType << " → effectiveChem="
                         << run.effectiveChem << "\n";
    } else {
        P.inOut->logMain << "WARNING: auto-detect " << result.detectedMatchMode
                         << " for " << run.featureType
                         << ", keeping inferred: " << run.effectiveChem << "\n";
    }
} else if (useAutodetect && !wlNamespaceConfident) {
    P.inOut->logMain << "WARNING: auto-detect " << result.detectedMatchMode
                     << " for " << run.featureType
                     << " but whitelist namespace is uncertain ('"
                     << wlNamespace << "' by default). effectiveChem kept as "
                     << run.effectiveChem << " to stay consistent with GEX. "
                     << "Use --crChemistry NXT/TRU or rename whitelist "
                     << "to enable absolute namespace resolution.\n";
}
```

**`prepared.inferredChem` / `prepared.inferredChemConfident`:** These come from
`detectChemistryFromWhitelistPath` (computed in `buildPfMultiPreparedContext`).
`confident` is `true` for content-based or filename-based detection, `false`
for the default-to-TRU fallback.  When confidence is false, auto-detect still
runs (for logging/diagnostics) but `effectiveChem` is **not overridden** —
all libraries stay in the same namespace as GEX to prevent merge mismatches.

**GEX effectiveChem:** Unchanged — detected from whitelist path/content.
Auto-detection runs on feature libraries only.

**Acceptance:** A375 smoke test passes (TRU whitelist + TRU reads → RAW_MATCH
→ effectiveChem=TRU, no behavior change). Log shows match mode + composition.

**Estimated effort:** 30 minutes

---

## Stage 5 — Synthetic End-to-End Test

**What:** Run the full pf-multi pipeline (through STAR) with synthetic data.
Tests the complete path including 2-column whitelist handling in the STAR wrapper
and the match-mode-to-chemistry composition.

**Test cases:**

| Case | Whitelist | Feature reads | Match mode | wlNamespace | effectiveChem |
|------|-----------|---------------|------------|-------------|---------------|
| A | 1-col TRU | TRU reads | RAW_MATCH | TRU | TRU |
| B | 1-col TRU | NXT reads | TRANSLATED_MATCH | TRU | NXT |
| C | 2-col NXT | NXT reads | RAW_MATCH | NXT | NXT |
| D | 2-col NXT | TRU reads | TRANSLATED_MATCH | NXT | TRU |

In all cases, merge-boundary normalization converts to TRU before merge.
Output barcodes are in TRU (or outputChem if `--crOutputChemistry` is set).

**Additional tests:**
- `--crChemistry NXT` explicit → auto-detect skipped, effectiveChem=NXT
- `--crChemistry TRU` explicit → auto-detect skipped, effectiveChem=TRU
- Low-support (total_hits < 50) → AMBIGUOUS, falls back to inferred
- Tiny library (< N reads) → drain finalization, correct match mode
- `--crOutputChemistry NXT` → output barcodes in NXT namespace

**Validation per case:**
- Check log for `NOTICE: Auto-detected match mode` with correct values
- Check output `barcodes.tsv.gz` — all in expected namespace
- Check feature counts — non-zero (merge joined correctly)

**Estimated effort:** 1–2 hours

---

## Stage 6 — Real-Data Validation (User-Provided)

**What:** Run on real datasets with known chemistry to confirm detection
agrees with ground truth.

**Candidate datasets:**
- UCSF iPSC2 (NXT 2-col whitelist, NXT reads) → RAW_MATCH → effectiveChem=NXT
- A375 (TRU 1-col whitelist, TRU reads) → RAW_MATCH → effectiveChem=TRU
- User-provided mixed-chemistry dataset (3' HT)

**Validation:**
- Detection matches known chemistry
- Output barcodes match prior validated runs
- No regression in cell counts, feature calls, or CRISPR calling

---

## Stage Summary

| Stage | What | Risk | Blocked by |
|-------|------|------|------------|
| 0 | Synthetic fixtures | Low | — |
| 1 | Detection state in pf_context + pf_api | Low | — |
| 2 | Trial-matching + drain in consumer loop | Medium | 1 |
| 3 | PfMultiAssign wiring | Low | 1 |
| 4 | PfMultiProcess integration (match mode → chem) | Medium | 2, 3 |
| 5 | Synthetic E2E test | Low | 0, 4 |
| 6 | Real-data validation | Low | 5 |

Stages 0 and 1 can be done in parallel.  
Stages 2 and 3 can be done in parallel after Stage 1.  
Total estimated effort: ~5–6 hours.

---

## Rollback / Feature Flag

- Auto-detection is gated by `autodetectChemistry` in `AssignOptions` (default
  false).
- `PfMultiProcess` only enables it when `requestedChem == "auto"` (the default
  for `--crChemistry`).
- When `--crChemistry` is explicitly `NXT` or `TRU`, auto-detection is skipped
  and the explicit value is used (existing behavior, no change).
- When detection returns `AMBIGUOUS` or `UNKNOWN`, the existing
  whitelist-inferred chemistry is used (no behavior change from today).
- `translate_NXT` inside process_features is never set by auto-detection.
  All actual translation happens at the merge boundary in `PfMultiProcess`.

---

## Both-Hit Barcodes

Positions 7-8 might already be self-complementary (e.g., AT→TA). Such barcodes
hit both raw and translated lookups. With typical 10x whitelist sizes (~700K
barcodes), self-complementary pairs are a small fraction and inflate both
counters equally, pushing toward 0.5 ratio. The 0.8/0.2 thresholds tolerate
this. The minimum-evidence gate (total_hits >= 50) further protects against
degenerate cases.
