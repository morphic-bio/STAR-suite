# Runbook: Multi-Feature Libraries (Lineage + Perturb) with Per-Library Chemistry and Dynamic Permits

Date: 2026-02-25  
Branch: `multi-feature`  
Status: Draft v1 (post-review)

## Objective
Support multiple feature libraries in one run (for example, lineage + perturb), where each library can independently use TRU or NXT chemistry and can use different whitelists and feature references, while remaining compatible with dynamic permit scheduling.

## What Is Already Implemented
1. Per-library chemistry override via `star_chemistry` in `pfMultiConfig` (`TRU`, `NXT`, `auto`, empty).
2. Per-library auto-detection of match mode (`RAW_MATCH`, `TRANSLATED_MATCH`, `AMBIGUOUS`) in `assignBarcodes`.
3. Safe merge-boundary normalization to a single output chemistry (`--crOutputChemistry`, default TRU).
4. Dynamic map/PF permit control for overlap and retuning.

## Current Gaps To Close
1. `pfMultiConfig` still carries one global feature reference (`[feature] ref`) and one global whitelist (`--crWhitelist` / `--soloCBwhitelist`) for all feature libraries.
2. No explicit per-library whitelist/ref contract in config parsing and prep validation.
3. Multi-library PF scheduling is still effectively one-library-at-a-time assignment runs.
4. Output provenance is not explicit enough to audit which whitelist/ref/chemistry each feature library used.

## Scope for This Runbook
1. Add per-library whitelist + feature reference support.
2. Keep existing single-library behavior unchanged.
3. Keep chemistry composition rules deterministic and backwards compatible.
4. Add a minimal multi-library PF scheduler policy compatible with existing permit control.

## Serial-First Execution Gate (Required)
1. Phases `0-3` are correctness-first and must run with serial/non-dynamic PF scheduling.
2. Do not enable new intra-PF permit allocation policy until Phase 3 acceptance criteria pass.
3. Required run mode during Phases `0-3`:
   - PF controller disabled (`off`)
   - existing serial assign/merge behavior
   - no new PF scheduler logic beyond current baseline
4. Phase 4 is the first phase allowed to change PF scheduling behavior.

## Non-Goals (This Slice)
1. No full rewrite to a global producer/consumer task queue.
2. No change to CRISPR calling model beyond ensuring correct per-library inputs.
3. No new feature-calling algorithms.

## Target Contract

### Config Contract (`pfMultiConfig`)
1. Add optional columns in `[libraries]`:
   - `star_feature_ref` (per-library feature reference CSV)
   - `star_library_id` (stable output/provenance key; auto-generated if absent)
2. `star_whitelist` is **deferred** to a future slice. All libraries in one
   GEM well share cell barcodes; NXT/TRU variants are already handled by
   `star_chemistry` + merge normalization. Per-library CB whitelists would
   require barcode-set union semantics at merge, which is out of scope.
3. Continue accepting existing files without these columns.
4. Precedence rules:
   - feature ref: `star_feature_ref` (per-library) > `--crFeatureRef` (CLI flag) > `[feature] ref` (config global)
   - chemistry request: `star_chemistry` (library) > `--crChemistry` > `auto`
5. When `star_feature_ref` is provided, skip `filterFeatureRefCsv` for that
   library (assume the per-library ref is already specific to its feature type).
6. Parse robustness requirements:
   - trim whitespace
   - case-insensitive headers
   - trailing-field padding (no row drops)
   - hard error on invalid chemistry token

### Chemistry Contract
1. Auto-detect remains match-mode only in C layer.
2. STAR layer composes absolute chemistry per library from:
   - library whitelist namespace (if confident)
   - detected match mode
   - explicit override path when configured
3. If whitelist namespace is uncertain:
   - do not force opposite namespace
   - keep library chemistry anchored to explicit/global chemistry
   - emit warning with actionable fix (`star_chemistry` or explicit whitelist naming)
4. Auto-detect validation is required in serial mode before permit work:
   - explicit cases: `TRU`, `NXT`, `auto`
   - mixed-library case: one explicit + one auto
   - ambiguous case: deterministic fallback + warning

### Output Contract
1. Keep merged outputs in one normalized namespace (`--crOutputChemistry`, default TRU).
2. Preserve feature type labels and library provenance in logs and metadata.
3. Add per-library provenance file under each assign output:
   - resolved fastq path
   - resolved feature ref
   - resolved whitelist
   - detected match mode
   - effective chemistry
4. Ensure merged feature matrix includes all configured feature libraries that produced valid outputs.

### Dynamic Permit Contract
1. Preserve existing global permit budget (`runThreadN`).
2. During map/PF overlap, retain map-vs-PF retuning behavior.
3. Inside PF domain, allocate work by library with coarse units:
   - at least one producer slot per active library when possible
   - grow consumer capacity per producer up to cap
   - when cap reached and files remain, add producer slot for another library
4. Single-thread hard fallback remains legacy sequential path.

## Implementation Phases

## Phase 0: Baseline Freeze and Fixtures
1. Freeze current single-library behavior and outputs as baseline.
2. Make `defaultPfFeatureSpecs()` data-driven: any non-GEX `feature_types`
   value in the config is routed to assignBarcodes. This is required for
   LARRY/Custom/Lineage library types that aren't in the hardcoded list.
   The hardcoded list becomes a fallback `featureRefType` when the config
   `feature_types` doesn't map to a known feature reference type.
3. Add/confirm synthetic fixture with two feature libraries:
   - library A: CRISPR Guide Capture + NXT
   - library B: Custom (lineage) + TRU
4. Create MSK 30polyKO downsampled fixture (100K reads: mRNA + PolyIII + LARRY).

Exit criteria:
1. Existing single-library tests unchanged (8/8 autodetect, 15/15 config).
2. Data-driven feature specs route LARRY/Custom libraries correctly.
3. Two-library fixture reproducible (data untracked; locations in `tests/ARTIFACTS.md`).

## Concrete Real-Data Fixture (MSK 30polyKO)
1. Input FASTQ root:
   - `/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM`
2. Feature tag references (per library):
   - gRNA (PolyIII): `/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv`
   - lineage (LARRY): `/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv`
3. Naming convention:
   - `LARRY_*_L00N_R{1,2}_001.fastq.gz` for lineage
   - `PolyIII_*_L00N_R{1,2}_001.fastq.gz` for gRNA
   - `mRNA_*_L00N_R{1,2}_001.fastq.gz` for expression
4. Sample token rule:
   - Treat the token between `<LABEL>_` and `_L00N_` as the sample id for parsing/grouping.
5. Initial smoke fixture rule:
   - choose one `R1/R2` pair from lane `L001` for each label (`mRNA`, `PolyIII`, `LARRY`)
   - create deterministic downsampled FASTQs by taking first `N` reads (for example `N=100000`)

Example fixture materialization:
```bash
set -euo pipefail
SRC=/mnt/pikachu/MSK-perturb/scRNAseq_30polyKO_ES_DE_XM
OUT=/tmp/msk_multifeature_fixture_$(date +%Y%m%d_%H%M%S)
NREADS=100000
NLINES=$((NREADS*4))
mkdir -p "$OUT"/{mRNA,PolyIII,LARRY}

for LIB in mRNA PolyIII LARRY; do
  R1=$(ls "$SRC/${LIB}_DE_*_L001_R1_001.fastq.gz" | head -n1)
  R2=${R1/_R1_/_R2_}
  zcat "$R1" | head -n "$NLINES" | gzip -c > "$OUT/$LIB/${LIB}_L001_R1.fastq.gz"
  zcat "$R2" | head -n "$NLINES" | gzip -c > "$OUT/$LIB/${LIB}_L001_R2.fastq.gz"
done

echo "Fixture root: $OUT"
```

6. Config intent for multi-feature tests:
   - `mRNA` row for GEX.
   - `PolyIII` row with `star_feature_ref=/mnt/pikachu/MSK-whitelists/ref_feature_geneBC.csv`.
   - `LARRY` row with `star_feature_ref=/mnt/pikachu/MSK-whitelists/ref_feature_larryBC.csv`.
   - set `star_chemistry` per library as needed (`auto`, `TRU`, or `NXT`) for test cases.

## Phase 1: Config/Data Model
1. Extend `PfMultiConfig::LibraryEntry` with `starFeatureRef` and `starLibraryId`.
2. Parse new columns (`star_feature_ref`, `star_library_id`) with aliases.
3. Implement feature ref precedence: `star_feature_ref` > `--crFeatureRef` > `[feature] ref`.
4. When `star_feature_ref` is provided, skip `filterFeatureRefCsv` (per-library
   ref is assumed to be already specific to the library's feature type).
5. Auto-generate `star_library_id` when absent: `{sampleName}_{featureTypes}_{index}`.
6. Fail-fast on: missing resolved feature ref, duplicate `star_library_id`.

Exit criteria:
1. Config parser unit tests cover old and new formats.
2. Invalid combinations fail with clear errors.
3. Per-library feature ref resolution logged in prep log.

## Phase 2: Per-Library Assign Wiring
1. Build each `PfPreparedFeatureLibrary` with its own resolved whitelist and feature ref.
2. Call `runAssignBarcodes` using per-library inputs.
3. Keep autodetect behavior unchanged except anchored to each library whitelist.
4. Validate auto-detect and chemistry composition per library in serial mode:
   - `detectedMatchMode` is correct per library
   - `effectiveChem` is composed correctly per library
   - no cross-library chemistry leakage

Exit criteria:
1. Mixed chemistry multi-library synthetic tests pass.
2. No regression on current UCSF/A375 single-library fixtures.
3. Serial MSK multi-feature fixture verifies per-library auto-detect behavior.

## Phase 3: Merge and Provenance Hardening
1. Merge all successful library MEX outputs after normalization to output namespace.
2. Add per-library provenance manifests and top-level summary table in `Log.out`.
3. Guard against silent partial merge:
   - if a configured library fails assign, record explicit failure and apply existing run policy (fail-fast or controlled skip).

Exit criteria:
1. Output includes both lineage and perturb features in merged matrices.
2. Provenance shows exact resolved inputs per library.

## Phase 4: PF Intra-Domain Scheduler (Permit-Compatible)
1. Introduce PF library scheduler state with remaining-work estimates per library.
2. Implement coarse allocation policy:
   - producer-first minimum per active library
   - consumer expansion up to `consumerCapPerProducer`
   - then open additional producers while files remain
3. Keep map/PF permit controller unchanged at outer boundary.

Exit criteria:
1. No deadlocks in 1-thread, low-thread, and multi-library high-thread runs.
2. PF uses available permits more fully after map completion.

## Phase 5: Validation and Rollout
1. Add multi-library parity scripts for:
   - raw merged feature matrix parity vs baseline fixture
   - filtered-on-CR-cells parity where applicable
   - chemistry normalization checks (output namespace consistency)
2. Run on:
   - synthetic two-library fixture
   - UCSF 2M downsampled set with two feature libraries
   - one full sample smoke
3. Record wall time and peak RAM vs current implementation.

Exit criteria:
1. Correctness gates pass.
2. Performance is neutral or improved.
3. Default mode unchanged for users not using new columns.

## Edge Cases (Must Be Explicitly Tested)
1. `runThreadN=1`: force legacy sequential path.
2. One feature library has zero readable FASTQs.
3. One library auto-detects `AMBIGUOUS` while another is explicit NXT/TRU.
4. Multiple libraries share same `feature_types` but different refs/whitelists.
5. Duplicate `star_library_id` values (must fail-fast).
6. Mixed confidence states for whitelist namespace inference.
7. One library fails assign while others succeed.
8. Per-library `star_feature_ref` that contains features of a different type than the library's `feature_types`.

## Proposed Test Matrix (Minimum)
1. Single-library regression (existing).
2. Two-library matched chemistries.
3. Two-library mixed TRU/NXT with auto-detect.
4. Two-library with one explicit chemistry and one auto.
5. Per-library feature ref override sanity.
6. Custom/Lineage feature_types routed via data-driven specs.
7. Low-thread stress (`runThreadN=1`, `2`, `4`).
8. Permit retune stress with map completion handoff.

## Rollout Flags and Compatibility
1. Add one gating flag for new per-library feature ref behavior if needed:
   - `--pfMultiPerLibraryInputs 1` (default off for first merge), then flip after validation.
2. Keep parser backward compatible regardless of flag.
3. Keep `star_chemistry` behavior unchanged.
4. `star_whitelist` remains deferred; parser ignores it if present.

## Deliverables
1. New parser/unit tests for config columns and precedence.
2. Multi-library synthetic E2E tests.
3. Updated docs:
   - `docs/feature_barcodes.md`
   - handoff summary for this slice
4. Benchmark summary note appended to dynamic-threads handoff.

## Resolved Questions
1. **`star_library_id`**: Auto-generate from `{sample}_{feature_types}_{index}` when
   absent. Allow explicit override for stable provenance keys.
2. **Single library failure**: Fail the whole run by default (existing behavior).
   Controlled-skip mode is a future extension.
3. **Per-library output chemistry**: No. Keep one global `--crOutputChemistry`
   (default TRU). Per-library output chemistry would produce mixed namespaces
   in the merged MEX, defeating normalization.
4. **`star_whitelist`**: Deferred. All libraries in a GEM well share cell barcodes;
   NXT/TRU is handled by `star_chemistry`. Per-library CB whitelists need
   barcode-set union at merge, which is out of scope for this slice.
