a# Handoff: UCSF Full CR-vs-Compat Gene Pearson Anomaly (2026-02-20)

## Status: ROOT-CAUSED AND FIXED

## Root cause: Poly-G artifact reads mapping to LINC00486

STAR lacked poly-G 3' tail trimming. NovaSeq/NextSeq two-color chemistry
produces poly-G artifact reads (no signal = G call). The reference genome
has a ~90bp G-homopolymer within the gene body of **LINC00486**
(ENSG00000230876, chr2:32,916,243, inside a 101kb lncRNA). Poly-G reads
align there perfectly (90M, MAPQ=255, NH:i:1), passing all quality filters
and inflating this single gene by 147x.

### Quantitative impact

| Metric | With LINC00486 | Without LINC00486 |
|---|---:|---:|
| Gene Pearson (all genes) | 0.304560 | **0.997871** |
| Gene Spearman (filtered) | 0.985121 | (unchanged) |
| LINC00486 UMIs (STAR) | 524,679 | -- |
| LINC00486 UMIs (CR) | 3,560 | -- |
| Delta from this one gene | +521,119 | -- |
| Share of total abs delta | 88% | -- |

Excluding LINC00486, gene Pearson recovers to 0.998. Excluding the top 5
delta genes gives 0.999.

### How we confirmed it

1. Per-gene delta table on common CR-filtered barcodes identified LINC00486
   as 88% of total absolute delta (524K STAR vs 3.5K CR).
2. BAM trace at chr2:32,916,000-32,917,000 showed 3,001 uniquely mapped
   reads piling up at a single position, 99.3% of which are poly-G
   (`GGGGGGGGG...`, 90bp).
3. All 3,316 poly-G reads in the entire 2M BAM map to chr2 at this locus.
4. The reference genome contains a >=90bp G-homopolymer at chr2:32,916,243
   inside LINC00486's gene body, providing a perfect alignment target.
5. Rescue mode is not involved: forward_yes vs forward_no 2M runs give
   identical LINC00486 counts (2,877 vs 2,879).
6. CellRanger trims poly-G; STAR's `--clipAdapterType CellRanger4` only
   trimmed poly-A (`ClipCR4::polyTail3p` checked base==0 only, not base==2).

### Fix applied

New parameter `--clip3pPolyG` (values: `yes`, `no`, `auto`; default `auto`):
- `yes` — always trim poly-G 3' tails.
- `no` — never trim poly-G 3' tails.
- `auto` — trim poly-G when `--clipAdapterType CellRanger4` is active.

`core/legacy/source/ClipCR4.cpp` — `polyTail3p` now trims poly-A
unconditionally (biological mRNA tails) and poly-G conditionally (gated
by the resolved `clipPolyG` flag). The poly-G scan uses the same scoring
algorithm as poly-A (20+ base run at >=70% purity, -2 mismatch penalty).

Parameter flow: `parametersDefault` → `Parameters.cpp` (parArray
registration) → `ParametersClip::initialize()` (resolve yes/no/auto →
bool `clipPolyG`) → `ParametersClip::initializeClipMates()` (propagate
to `ClipCR4::clipPolyG`) → `polyTail3p()` (gate the G scan).

**TODO**: Implement `auto` mode with data-driven detection (count poly-G
fraction in first N reads, similar to fastp's approach). Currently `auto`
is equivalent to `yes` when `--clipAdapterType CellRanger4` is used.

CellRanger comparison (2M BAM trace): CR's BAM at the LINC00486 hotspot
has 138 reads (118 poly-G) vs STAR's 4,335 (4,257 poly-G). Critically,
**zero pure poly-G reads** (100% G content) appear in CR's BAM despite
identical 1M GEX input FASTQs, confirming CR discards them in its Rust
preprocessing upstream of STAR, not via STAR's polyA trimming path (no
`pa:i` tags on surviving reads). On this NovaSeq data, no poly-C or
poly-T artifacts were observed in the STAR BAM at any threshold.

## Original anomaly (preserved for context)

### What was run (pinned)
Run script (committed):
- `scripts/run_ucsf_full_compat_forward_rescue_guides.sh`

Parity script (committed):
- `scripts/run_ucsf_full_compat_forward_rescue_parity.sh`

Run output:
- STAR run: `/storage/ucsf-full/bench_20260218_dynamic_first/runs/star_full_iPSC2_1_AALG2_forward_rescue_guides_20260220_0825`
- Parity output: `/tmp/ucsf_full_parity_after_rerun_20260220_085841/PARITY_GEX_FEATURES_RAW_AND_CR_FILTERED.txt`

Important pinned params:
- `--runThreadN 32`
- `--soloStrand Forward`
- `--soloCrMultimapRescue yes`
- `--soloCBwhitelist .../3M-february-2018_NXT.txt`
- `--crWhitelist .../3M-february-2018_NXT.txt`
- guides included (`--pfMultiConfig`, `--crFeatureRef`, `--crAssign*`)

### GEX metrics vs Cell Ranger (pre-fix)
| Subset | Common barcodes | Cell-level Pearson | Gene Pearson (all genes) | Gene Pearson (filtered genes) | Gene Spearman (filtered genes) |
|---|---:|---:|---:|---:|---:|
| `raw_all` | 196,669 | 0.982826 | 0.320055 | 0.260106 | 0.905978 |
| `raw_restricted_to_cr_filtered` | 7,325 | 0.946994 | 0.304560 | 0.301089 | 0.985121 |
| `filtered_vs_filtered` | 6,891 | 0.945690 | 0.296456 | 0.293051 | 0.982762 |

### Why this was suspicious
- Barcode overlap is high (`6891` common filtered barcodes; Jaccard `0.927456`), so this is not a barcode-intersection failure.
- Cell-level Pearson is high (>0.94), but gene Pearson stays ~0.30.
- Gene Spearman is high (~0.98), classic signature of a few extreme outlier genes distorting Pearson but not rank order.
- STAR GEX totals are consistently higher than CR on common barcodes:
  - `raw_all`: STAR `2,622,570` vs CR `1,957,055` (delta `+649,566`)
  - `raw_restricted_to_cr_filtered`: STAR `1,862,153` vs CR `1,350,854` (delta `+511,292`)
  - `filtered_vs_filtered`: STAR `2,116,069` vs CR `1,553,777` (delta `+587,175`)

## Secondary inconsistency (feature side)
Same run shows:
- Feature raw counts are close to CR (`19,695,899` vs `19,773,321`).
- But `Feature FILTERED vs FILTERED` STAR total collapses to `253,506` with only `1,177` STAR filtered feature barcodes and `88` common filtered barcodes.
- Additional metrics feature-call section also reports low STAR rows (`rows_star=1177`, `total_assigned_umis_star=253038`).

This indicates a likely filtered-feature/call path issue distinct from raw feature counting. Tracked separately.

## Fixes applied and validated (full sample)

### Fix 1: Poly-G trimming (2026-02-20)
`--clip3pPolyG auto` eliminates NovaSeq poly-G artifacts. Gene Pearson
recovered from 0.30 to 0.997 on the full UCSF iPSC2 dataset.

### Fix 2: EmptyDrops BH correction + umiMin (2026-02-20)
Non-Flex libscrna path: `umiMin` 500→100, BH FDR correction enabled.
Rescued 8 additional cells; Flex defaults unchanged.

### Fix 3: Bootstrap OrdMag knee estimation (2026-02-21)
Non-Flex libscrna path now uses CR9-style bootstrapped `recovered_cells`
estimation instead of hardcoded `nExpectedCells=3000`. The bootstrap
(100 samples) searches a log2-spaced range of candidate cell counts and
minimizes `(filtered - recovered)^2 / recovered`, matching Cell Ranger's
`filter_cellular_barcodes_ordmag()`.

MC simulations bumped from 10,000 to 100,000 (matching CR9 `NUM_SIMS`).

Files changed:
- `core/features/libscrna/include/scrna_api.h` — `use_bootstrap` config field
- `core/features/libscrna/src/scrna_api.cpp` — dispatch to bootstrap path
- `core/legacy/source/SoloFeature_emptyDrops_libscrna.cpp` — enable bootstrap
  for non-Flex, set `n_expected_cells=0`, `sim_n=100000`

Flex path (FlexFilter.cpp) is untouched; it keeps `nExpectedCells=3000`.

### Validated results: UCSF iPSC2 full sample (2026-02-21)

Run: `star_full_iPSC2_1_AALG2_forward_rescue_guides_bootstrap_20260221_010635`

| Metric | Before (fixed 3000) | After (bootstrap) | CR9 reference |
|--------|--------------------|--------------------|---------------|
| Filtered cells | 6,609 | **7,286** | 7,325 |
| Cell count gap vs CR | 716 | **39** | — |
| Jaccard (cell overlap) | — | **0.9898** | — |
| Only STAR | — | 18 | — |
| Only CR | — | 57 | — |
| Gene Pearson (all) | 0.997 | **0.9969** | — |
| Gene Pearson (filtered) | — | **0.9969** | — |
| Gene Spearman (filtered) | — | **0.9789** | — |
| Cell Pearson | — | **0.9987** | — |
| Cell Spearman | — | **0.9983** | — |

Bootstrap estimated `recovered_cells = 7,263` (CR: ~7,325), lowering the
OrdMag UMI cutoff from 111 to 88. All 6,609 previously called cells are
retained (proper superset). The 716-cell gap closed to 39.

## Remaining items
1. Investigate remaining top-delta genes (ENSG00000287001, ENSG00000280441, etc.).
2. Investigate feature-filtered barcode collapse (secondary inconsistency above).
3. Implement data-driven poly-G `auto` mode (count poly-G in first N reads).
