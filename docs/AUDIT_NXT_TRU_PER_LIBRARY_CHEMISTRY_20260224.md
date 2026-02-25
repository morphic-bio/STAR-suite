# Audit: NXT/TRU Per-Library Chemistry Support

**Date:** 2026-02-24  
**Scope:** STAR-suite codebase audit for single global vs. per-library NXT/TRU chemistry assumptions  
**Critical for:** 3' HT experiments where GEX libraries use TRU barcodes and feature libraries use NXT barcodes

---

## Executive Summary

The codebase **assumes a single global NXT/TRU chemistry setting** throughout the pf-multi pipeline. For 3' HT experiments with mixed chemistries (GEX=TRU, Feature=NXT), this is **blocking**. Key findings:

| Component | Global/Per-Lib | Blocking? | Notes |
|-----------|----------------|-----------|-------|
| PfMultiProcess | **Global** | **Yes** | Single effectiveChem/outputChem, single whitelist |
| PfMultiMerge | **Global** | **Yes** | Single input/output chemistry applied to ALL barcodes |
| PfMultiAssign | **Global** | **Yes** | No translate_nxt passed to pf_api; single whitelist |
| PfMultiConfig | **N/A** | **Yes** | No per-library whitelist or chemistry |
| PfMultiMexStub | **Global** | Partial | Uses single whitelist for NXT→TRU mapping |
| call_features.c | **N/A** | No | Barcode-agnostic |
| gmm.c | **N/A** | No | Barcode-agnostic |
| star_feature_call | **Global** | No | Standalone; no NXT/TRU options |

---

## 1. PfMultiProcess.cpp

### 1.1 Chemistry Storage (Lines 50–66, 496–497)

```cpp
struct PfMultiPreparedContext {
    ...
    string effectiveChem;
    string outputChem;
    ...
};

struct PfMultiPreloadInput {
    ...
    string crChemistry;
    string crOutputChemistry;
    ...
};
```

**Finding:** Single global `effectiveChem` and `outputChem` per run. No per-library chemistry fields.

### 1.2 Chemistry Inference and Assignment (Lines 804–846)

```cpp
static PfMultiPreparedContext buildPfMultiPreparedContext(const PfMultiPreloadInput& input) {
    ...
    context.whitelist = input.crWhitelist;  // Single whitelist
    ...
    context.requestedChem = parseChemistryToken(input.crChemistry, "--crChemistry");
    context.requestedOutputChem = parseChemistryToken(input.crOutputChemistry, "--crOutputChemistry");
    context.inferredChem = detectChemistryFromWhitelistPath(context.whitelist, context.inferredReason);
    context.effectiveChem = context.inferredChem;
    if (context.requestedChem == "nxt") {
        context.effectiveChem = "NXT";
    } else if (context.requestedChem == "tru") {
        context.effectiveChem = "TRU";
    }
    context.outputChem = "TRU";  // Default
    if (context.requestedOutputChem == "nxt") {
        context.outputChem = "NXT";
    } else if (context.requestedOutputChem == "tru") {
        context.outputChem = "TRU";
    }
    ...
}
```

**Finding:** Chemistry is inferred from a **single whitelist path** (`detectChemistryFromWhitelistPath`). All libraries share this.

### 1.3 Usage in assignBarcodes Loop (Lines 1116–1117, 1467, 1771–1786)

```cpp
const string& effectiveChem = prepared.effectiveChem;
const string& outputChem = prepared.outputChem;
...
ret = PfMultiAssign::runAssignBarcodes(whitelist, refPath, resolvedFastq, assignOut, runAssignOpts);
...
ret = PfMultiMerge::writeCombinedMex(rawOutDir, ..., effectiveChem, outputChem);
ret = PfMultiMerge::writeCombinedMex(filteredOutDir, ..., effectiveChem, outputChem);
```

**Finding:** Same `effectiveChem`/`outputChem` used for all feature libraries and for both raw and filtered MEX writes.

### 1.4 Required Changes for Per-Library Chemistry

- Add per-library `whitelist` and `chemistry` (or `effectiveChem`/`outputChem`) to `PfMultiConfig::LibraryEntry` or equivalent.
- Extend `PfPreparedFeatureLibrary` with `effectiveChem`, `outputChem`, and `whitelist`.
- Pass per-library chemistry into `runAssignBarcodes` and into merge/write logic.
- Extend multi config format to support per-library `whitelist` and `chemistry`.

---

## 2. PfMultiMerge.cpp

### 2.1 writeCombinedMex Signature (Lines 484–490)

```cpp
int writeCombinedMex(const string& outputDir,
                     const MexData& data,
                     const string& gemWell,
                     ofstream& logStream,
                     const vector<string>& gexBarcodes,
                     const string& inputChemistry,
                     const string& outputChemistry)
```

**Finding:** Single `inputChemistry` and `outputChemistry` for the entire combined MEX.

### 2.2 Translation Logic (Lines 584–596, 671–680)

```cpp
// Step 2: Apply optional barcode namespace translation.
vector<string> namespaceAdjustedBarcodes = filteredBarcodes;
if (needsNamespaceTranslation(inputChemistry, outputChemistry)) {
    for (size_t i = 0; i < namespaceAdjustedBarcodes.size(); ++i) {
        string translated = namespaceAdjustedBarcodes[i];
        translateNxtMiddleTwoBasesInplace(translated);
        ...
    }
}
...
// Derive native/read-space barcode representation
const bool emitNamespaceArtifacts = needsNamespaceTranslation(inputChemistry, outputChemistry);
if (emitNamespaceArtifacts) {
    for (auto& bc : nativeSuffixedBarcodes) {
        translateNxtMiddleTwoBasesInplace(bc);
    }
}
```

**Finding:** Translation is applied uniformly to **all** barcodes. In 3' HT, GEX barcodes are already TRU; feature barcodes may be NXT. Applying NXT→TRU to everything would incorrectly translate GEX barcodes.

### 2.3 mergeMex Barcode Matching (Lines 417–441)

```cpp
// Create barcode map (barcode -> column index in GEX)
map<string, uint32_t> barcodeMap;
for (size_t i = 0; i < merged.barcodes.size(); ++i) {
    barcodeMap[merged.barcodes[i]] = i;
}
...
for (const auto& t : featData.triplets) {
    string bc = featData.barcodes[t.cell_idx];
    auto it = barcodeMap.find(bc);
    if (it != barcodeMap.end()) {
        // Merge triplet
    } else {
        missingCount++;  // Feature barcode not in GEX - DROPPED
    }
}
```

**Finding:** Feature barcodes are matched to GEX by exact string. If GEX is TRU and feature is NXT, they will not match and feature counts are dropped.

### 2.4 Raw vs Filtered MEX

**Finding:** Both raw and filtered MEX use the same `effectiveChem` and `outputChem` (PfMultiProcess.cpp 1771–1786). No separate chemistry for raw vs filtered.

### 2.5 Required Changes for Per-Library Chemistry

- Track chemistry per feature type or per MEX source.
- In `writeCombinedMex`, translate only barcodes from NXT feature libraries to TRU before merging.
- Alternatively, translate feature barcodes during merge (before `barcodeMap.find`) so GEX (TRU) and feature (NXT→TRU) align.
- Ensure GEX barcodes are never translated when they are already in output namespace.

---

## 3. PfMultiMerge.h

### 3.1 writeCombinedMex Declaration (Lines 70–80)

```cpp
 * @param inputChemistry Effective correction chemistry (NXT/TRU)
 * @param outputChemistry Output barcode namespace (NXT/TRU)
 ...
int writeCombinedMex(..., const string& inputChemistry = "TRU",
                     const string& outputChemistry = "TRU");
```

**Finding:** API assumes a single input and output chemistry for the whole MEX.

---

## 4. PfMultiMexStub.cpp / PfMultiMexStub.h

### 4.1 copyBarcodesTsv (Lines 254–298, 342)

```cpp
bool copyBarcodesTsv(const string& barcodesTxt, const string& barcodesTsv, bool force,
                     const string& whitelistPath = "");
...
if (copyBarcodesTsv(barcodesTxt, barcodesTsv, force, whitelistPath)) {
```

**Finding:** Uses a single `whitelistPath` for NXT→TRU mapping (2-column whitelist). Called from `processAssignOutput` with the global `whitelist` (PfMultiProcess.cpp:1573).

### 4.2 processAssignOutput (Lines 301–357)

```cpp
int processAssignOutput(const string& assignOutDir, const string& featureCsvPath,
                       const string& defaultFeatureType, bool force,
                       const string& whitelistPath = "");
```

**Finding:** One whitelist for all feature libraries. For per-library chemistry, each library would need its own whitelist (and possibly mapping).

---

## 5. PfMultiAssign.cpp / PfMultiAssign.h

### 5.1 AssignOptions (PfMultiAssign.h:15–30)

```cpp
struct AssignOptions {
    int maxHammingDistance = -1;
    int featureConstantOffset = -1;
    ...
    string filteredBarcodesPath;
};
```

**Finding:** No `translate_nxt` or chemistry-related options.

### 5.2 applyAssignOptions (PfMultiAssign.cpp:138–187)

**Finding:** Does not call `pf_config_set_translate_nxt`. Chemistry/translate_nxt is never set for pf_api.

### 5.3 runAssignBarcodes (PfMultiAssign.cpp:293–396)

```cpp
int runAssignBarcodes(const string& whitelist,
                     const string& featureRef, const string& fastqDir,
                     const string& assignOut,
                     const AssignOptions& options) {
    ...
    applyAssignOptions(cfg, options);
    ...
    pf_load_whitelist(ctx, whitelistForAssign.c_str());
    ...
    pf_process_fastq_dir(ctx, fastqDir.c_str(), assignOut.c_str(), &stats);
```

**Finding:** Single whitelist per call; `pf_config_set_translate_nxt` is never invoked. assignBarcodes always uses `translate_nxt=0`.

### 5.4 Required Changes for Per-Library Chemistry

- Add `translateNxt` (or equivalent) to `AssignOptions`.
- Call `pf_config_set_translate_nxt(cfg, options.translateNxt ? 1 : 0)` in `applyAssignOptions`.
- In PfMultiProcess, set `assignOpts.translateNxt` per library based on that library’s chemistry.
- Support per-library whitelist (e.g., NXT for feature, TRU for GEX).

---

## 6. Feature Calling: call_features.c

### 6.1 Search Results

**Finding:** No references to `translate_NXT`, `NXT`, `TRU`, or `chemistry`. Barcodes are used as opaque strings. **Compatible** with per-library chemistry as long as merged MEX barcodes are in a consistent namespace.

---

## 7. CRISPR Calling / GMM: gmm.c

### 7.1 Search Results

**Finding:** No references to NXT, TRU, or chemistry. **Compatible** with per-library chemistry.

---

## 8. star_feature_call.cpp

### 8.1 Configuration (Lines 32–69)

```cpp
struct Config {
    ...
    std::string whitelist;
    ...
    // No translate_nxt, chemistry, or NXT/TRU options
};
```

**Finding:** No NXT/TRU or chemistry options. Uses single whitelist. Standalone tool; not part of pf-multi. **Compatible** for single-chemistry runs; would need options for mixed-chemistry use cases.

---

## 9. process_features (assignBarcodes, pf_api, mex_writer)

### 9.1 pf_api.c (Lines 51, 110, 233–234, 318)

```cpp
struct pf_config {
    ...
    int translate_nxt;
    ...
};
config->translate_nxt = 0;  // Default
void pf_config_set_translate_nxt(pf_config *config, int enable) { ... }
translate_NXT = ctx->config->translate_nxt;  // Applied at pf_init
```

**Finding:** `translate_nxt` exists in pf_api but is never set from PfMultiAssign. Default remains 0.

### 9.2 assignBarcodes.c / mex_writer.c

**Finding:** Use global `translate_NXT` for NXT→TRU at output/filter stages. Behavior is correct when `translate_NXT` is set; the gap is that PfMultiAssign does not set it.

---

## 10. PfMultiConfig

### 10.1 LibraryEntry (PfMultiConfig.h:24–32)

```cpp
struct LibraryEntry {
    string fastqs;
    string feature_types;
    string sample;
    string library_type;
    string gem_well;
    // No whitelist, chemistry, or translate_nxt
};
```

**Finding:** No per-library whitelist or chemistry. Config format would need extensions for per-library chemistry support.

---

## Summary Table: Blocking vs Compatible

| Location | File:Line | Global/Per-Lib | Blocking? | Change Required |
|----------|-----------|----------------|-----------|----------------|
| Chemistry context | PfMultiProcess.cpp:59-60 | Global | **Yes** | Per-library effectiveChem/outputChem |
| Chemistry inference | PfMultiProcess.cpp:831-846 | Global | **Yes** | Per-library inference or override |
| Whitelist | PfMultiProcess.cpp:820-828 | Global | **Yes** | Per-library whitelist in config |
| assignBarcodes call | PfMultiProcess.cpp:1467 | Global | **Yes** | Pass per-library whitelist + translate_nxt |
| writeCombinedMex (raw) | PfMultiProcess.cpp:1766-1772 | Global | **Yes** | Per-source chemistry or pre-merge translation |
| writeCombinedMex (filtered) | PfMultiProcess.cpp:1780-1786 | Global | **Yes** | Same as raw |
| writeCombinedMex impl | PfMultiMerge.cpp:584-596 | Global | **Yes** | Translate only NXT-sourced barcodes |
| mergeMex | PfMultiMerge.cpp:437-441 | Implicit | **Yes** | Translate feature barcodes before match |
| processAssignOutput | PfMultiProcess.cpp:1573 | Global | Partial | Per-library whitelist for copyBarcodesTsv |
| AssignOptions | PfMultiAssign.h:15-30 | N/A | **Yes** | Add translateNxt |
| applyAssignOptions | PfMultiAssign.cpp:138-187 | N/A | **Yes** | Call pf_config_set_translate_nxt |
| PfMultiConfig | PfMultiConfig.h:24-32 | N/A | **Yes** | Add whitelist, chemistry per library |
| call_features.c | - | N/A | No | Barcode-agnostic |
| gmm.c | - | N/A | No | Barcode-agnostic |
| star_feature_call | star_feature_call.cpp | Global | No | Optional: add --translate-nxt, --chemistry |

---

## Recommended Implementation Order

1. **PfMultiConfig**: Add `whitelist` and `chemistry` (or `effective_chem`/`output_chem`) to `LibraryEntry`; extend config parser.
2. **PfMultiProcess**: Build per-library chemistry and whitelist; pass them into assign and merge.
3. **PfMultiAssign**: Add `translateNxt` to `AssignOptions`; call `pf_config_set_translate_nxt`; accept per-library whitelist.
4. **PfMultiMerge**: Before `mergeMex`, translate feature barcodes from NXT→TRU when feature chemistry is NXT and output is TRU. In `writeCombinedMex`, do not translate barcodes that are already in output namespace (or pass per-source chemistry and translate selectively).
5. **PfMultiMexStub**: Use per-library whitelist in `processAssignOutput` when available.

---

## Resolution (2026-02-24)

### Approach: Normalize to TRU at the integration boundary

Rather than translating per-library in assignBarcodes and tracking namespaces through
the pipeline, a simpler design was implemented: all barcode namespaces are normalized
to **TRU** at the merge boundary (before `mergeMex`). This means:

1. Each library's `assignBarcodes` works in read-space (no `translate_NXT`)
2. `normalizeMexBarcodesToTru(data, effectiveChem)` translates NXT→TRU for any
   NXT-sourced MexData before cross-library merge
3. `writeCombinedMex` receives `"TRU"/"TRU"` → no further translation needed
4. Output is always in TRU namespace

### Changes Made

| File | Change |
|------|--------|
| `PfMultiAssign.h` | Added `bool translateNxt` to `AssignOptions` (wired for future per-library use) |
| `PfMultiAssign.cpp` | `applyAssignOptions` calls `pf_config_set_translate_nxt` when set; logged in run summary |
| `PfMultiProcess.cpp` | Added `effectiveChem` to `PfPreparedFeatureLibrary` and `FeatureRun`; added `normalizeMexBarcodesToTru()` helper; normalize all MexData to TRU before merge; pass `"TRU"/"TRU"` to `writeCombinedMex` |

### Validation

- **A375 (TRU path)**: Full end-to-end CRISPR smoke test passed — no behavioral change (normalization is a no-op for TRU inputs).
- **UCSF iPSC2 (NXT path)**: 7,286/7,286 barcode NXT→TRU translations verified against prior run output. 100% match.

### Remaining Work

- Auto-detection of NXT/TRU per library during bootstrap phase
- Synthetic test harness for mixed-chemistry experiments (GEX=TRU + Feature=NXT)
- Per-library whitelist support in `PfMultiConfig`

---

## References

- `docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md`
- `core/features/process_features/README.md` (--translate_NXT)
- `parametersDefault` (crChemistry, crOutputChemistry, crWhitelist)
