# CAT-ATAC Trimodal Guide Arm Runbook

## Goal

Support the CAT-ATAC guide-capture library as the third modality in a single
STAR-suite multiome run:

- GeneFull GEX from STARsolo.
- ATAC peak MEX from libchromap.
- CRISPR guide count/call matrix from `process_features`.

The missing support is not guide calling itself. It is the guide-arm read layout:
cell barcode, UMI, and protospacer are split across three different reads.

```text
R1: capture sequence plus UMI, UMI = R1[0:12]
R2: 24 bp 10x ARC ATAC i5 barcode read, CB = revcomp(R2[8:24])
R3: protospacer read, match spacer anywhere in R3
```

The barcode transform must be the same contract used by the existing chromap
ATAC path:

```text
STAR/chromap read-format: bc:8:23:-
0-based half-open:        R2[8:24)
inclusive bases:          0-based 8-23
```

## Hard Invariants

- Do not use the literal awk interpretation of `substr($seq,9,24)` as a 24 bp
  barcode. The correct barcode is 16 bp. The fixture validates this:
  `revcomp(R2[8:24])` gives high ATAC whitelist hit rate; offset 9 or length 24
  gives zero hits.
- Reuse the chromap barcode layout contract (`bc:8:23:-`) for guide R2. If the
  implementation cannot call chromap internals directly, share a small parser or
  adapter that consumes the same layout string.
- Correct observed guide barcodes in ATAC barcode namespace first, then translate
  to GEX barcode namespace for MEX output using `atac2gex.tsv`.
- Keep CAT-ATAC out of the NXT/TRU middle-two-base translation path. This is an
  ATAC-to-GEX two-column map, not a 10x NXT/TRU namespace conversion.
- Match guides as free k-mers anywhere in R3 with `max_hamming <= 1`, using each
  spacer's actual length. The CRISPRi library has both 19 bp and 20 bp spacers.
- Accept either capture sequence in R1:
  - `CS1=CAAGTTGATAACGGACTAGCC`
  - `CS2=CAAGTTGTAAACGGACTAGCC`
- Compare guide count matrices before comparing per-cell guide assignments.
  Assignment concordance includes caller/model differences; count concordance
  tests extraction, matching, correction, and UMI dedup directly.
- Feature-reference duplicate names are accepted with first-definition-wins
  semantics. Later rows with the same feature name are skipped and emit a
  warning (for example `HIC2_1` / `HIC2_2` in the CAT-ATAC CSV). Do not edit
  the CSV to remove duplicates and do not special-case individual guide names
  in the loader.
- Mixed-length spacers are a correctness requirement, not a performance
  blocker for this library. The CAT-ATAC CRISPRi reference has only 54 unique
  features, so the existing mixed-length fallback is acceptable for the smoke
  and concordance runs.

## Data

- Handoff: `docs/HANDOFF_CATATAC_TRIMODAL_GUIDE_ARM_20260612.md`
- Fixture: `/mnt/pikachu/catatac_gse288996/guide_redump/fixture/`
  - `guide_R1.fastq.gz`
  - `guide_R2.fastq.gz`
  - `guide_R3.fastq.gz`
  - `MANIFEST.txt`
- Full guide FASTQs:
  `/mnt/pikachu/catatac_gse288996/guide_redump/SRR32265756_{1,2,3}.fastq.gz`
- CRISPRi spacers:
  `/mnt/pikachu/catatac_gse288996/guide_ref/ps_ref_crispri.tsv`
- STAR-suite guide feature ref (56 CSV rows, 54 unique names after loader dedup):
  `core/features/process_features/feature_lists/catatac_crispri_guide_capture.csv`
- ATAC whitelist:
  `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt`
- ATAC-to-GEX map:
  `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/atac2gex.tsv`
- GEX whitelist:
  `/mnt/pikachu/GEX_whitelist/737K-arc-v1.txt`
- Concordance target:
  `/mnt/pikachu/catatac_gse288996/catatac_public/guide_caller_calls/DMSO1/DMSO1_calls.rds`
- Full trimodal GEX/ATAC references must stay on the 2020-A ARC reference set;
  do not mix 2024-A GEX with 2020-A ATAC for this benchmark.

## Existing Code To Inspect First

- `core/features/process_features/src/assignBarcodes.c`
  - `checkAndCorrectBarcode()` currently assumes CB+UMI are contiguous in the
    same barcode read.
  - Feature matching paths include exhaustive/free search and per-feature
    lengths, but mixed 19/20 bp spacer behavior needs a regression test.
  - Current cumulative feature prehash is same-length only. If a future guide
    library is large and mixed-length, add length-bucketed prehashes instead of
    creating an assay-specific matcher.
- `core/features/process_features/include/pf_api.h`
  - `pf_read_record_view` has `barcode_sequence`, `feature_sequence`, and
    optional `feature_sequence2`; it does not currently model split CB and UMI
    sources.
- `core/features/process_features/src/pf_api.c`
  - In-memory record processing can be reused if the producer synthesizes an
    in-memory `CB+UMI` barcode sequence.
- `core/legacy/source/PfMultiAssign.cpp`
  - STAR pf-multi currently calls the process-features API by FASTQ directory
    or CBQ views, without per-library guide layout fields.
- `core/legacy/source/PfMultiConfig.h`
  and `core/legacy/source/PfMultiConfig.cpp`
  - Add per-library layout/config columns here.
- `core/legacy/source/PfMultiProcess.cpp`
  - Wire parsed per-library layout into `AssignOptions`.
  - Generalize two-column barcode output mapping beyond NXT/TRU.
  - Avoid forcing anchor-search behavior for unanchored CAT-ATAC spacers.
- `core/legacy/source/PfMultiMexStub.cpp`
  - `copyBarcodesTsv()` already supports `col1 -> col2` output mapping when
    given a two-column whitelist/map. Use this after ATAC-space correction.
- `core/features/libchromap_contract/include/star_chromap_contract.h`
  and `core/features/libchromap_contract/src/star_chromap_contract.cpp`
  - Chromap receives the ATAC barcode transform as `read_format`.

## Implementation Plan

### Future Scaling TODO: Length-Bucketed Feature Prehash

This is not required for the 54-guide CAT-ATAC run, but it is the clean
generalization for large mixed-length guide or lineage libraries.

Current `process_features` behavior:

- Exact matching supports mixed lengths through per-feature lengths and the
  `mismatched_feature_indices` fallback.
- Cumulative Hamming prehash (`<=1`, `<=2`) is enabled only when all features
  share one `common_length`.
- If any feature length differs, mixed-length matching remains correct but uses
  fallback search for fuzzy matching.

Desired implementation:

- Partition features by `feature_lengths[i]`.
- Build one exact / `<=1` / `<=2` prehash family per length bucket.
- At query time, probe only lengths that fit the candidate read window. For
  free R3 search, scan candidate windows by observed length bucket; for fixed
  or learned offsets, apply the same offset set per bucket.
- Preserve existing ambiguity semantics within and across buckets: if two
  features at the same best distance compete, return no-call with the best
  distance, not an arbitrary feature.
- Keep memory budgeting per bucket and report which buckets build or skip.
- Add a regression with 19 bp and 20 bp guides proving matrix parity with the
  fallback path while confirming the prehash path is used.

STAR-Spatial has a useful local analogue for length-bucketed candidate search:
`/mnt/pikachu/STAR-Spatial/native/hd_sw_rescue_harness.cpp` groups oligos by
length and evaluates only candidate length buckets within the edit-distance
window (`group_oligos_by_length()`, `find_half_hits()`). Reuse the idea, not
the spatial-specific code.

### Phase -1: Branch And Worktree

Develop this on a dedicated feature branch:

```text
feature/catatac-trimodal-guide-arm
```

This checkout may contain unrelated local work. Do not clean, reset, or switch
away from a dirty working tree. Prefer a fresh worktree from `master`:

```bash
git fetch origin
git worktree add /mnt/pikachu/STAR-suite-catatac-trimodal-guide-arm \
  -b feature/catatac-trimodal-guide-arm origin/master
cd /mnt/pikachu/STAR-suite-catatac-trimodal-guide-arm
```

If the branch or worktree already exists, reuse it after checking `git status`.
Follow the repo policy: raw development stays on `feature/*`; release-intended
work later merges to `dev-release`, then `master`.

### Phase 0: Build An Oracle

Before changing STAR internals, create a small oracle that streams the fixture
and validates the exact extraction rules:

1. Read synchronized `guide_R1/R2/R3`.
2. Build an in-memory or temporary barcode sequence:
   `revcomp(R2[8:24]) + R1[0:12]`.
3. Use R3 as the feature read.
4. Use ATAC whitelist column 1 for barcode correction.
5. Use `atac2gex.tsv` only for output barcode translation.
6. Produce a guide count matrix with standalone `assignBarcodes`.

This oracle is not the production path. It is a correctness gate for native
split-read support.

Required oracle metrics:

```text
total_reads
atac_whitelist_exact_hit_rate_for_revcomp_R2_8_24
atac_whitelist_exact_hit_rate_for_revcomp_R2_9_25
atac_whitelist_exact_hit_rate_for_revcomp_R2_8_32
capture_cs1_hit_rate
capture_cs2_hit_rate
capture_either_hit_rate
r3_spacer_match_rate
umi_valid_rate
deduped_guide_umi_count
gex_space_output_barcode_count
```

The offset checks should show high hit rate for `R2[8:24]` and near zero for the
negative controls.

### Phase 1: Feature Reference Conversion

Convert `ps_ref_crispri.tsv` to a process-features feature ref:

```text
id,name,read,pattern,sequence,feature_type
ADNP_1,ADNP_1,R3,(BC),CACCCTCTCCGCCGAAGTG,CRISPR Guide Capture
```

Rules:

- Keep guide names (`GENE_1`, `GENE_2`) as separate features.
- Preserve 19 bp and 20 bp spacer lengths.
- Use free R3 search, not a fixed offset.
- Set guide max hamming to 1 for parity with the public `bbduk` caller.
- If `pattern=(BC)` triggers an unwanted fixed-anchor assumption, add an
  explicit config path for unanchored/free search rather than encoding a fake
  anchor.

### Phase 2: Native Split-Read Layout Support

Add a layout adapter that can synthesize the existing process-features barcode
sequence in memory:

```text
synth_barcode_sequence = transformed_CB + UMI
synth_barcode_quality  = transformed_CB_quality + UMI_quality
```

For CAT-ATAC:

```text
CB source:      R2
CB offset:      8
CB length:      16
CB strand:      reverse-complement
UMI source:     R1
UMI offset:     0
UMI length:     12
Feature source: R3
```

Quality handling:

- Reverse the R2 quality substring when reverse-complementing the barcode.
- Append R1 UMI quality after the barcode quality.
- Keep quality length equal to barcode length plus UMI length.

Implementation preference:

- Keep `checkAndCorrectBarcode()` and downstream UMI counting largely unchanged.
- Add a producer/adapter that emits the same `CB+UMI` shape the existing
  consumer already expects.
- Do not materialize derived FASTQ files in the production STAR path.

### Phase 3: Config Plumbing

Add per-library fields to pf-multi config. Suggested generic fields:

```text
star_barcode_read
star_barcode_format
star_umi_read
star_umi_start
star_umi_length
star_feature_read
star_capture_read
star_capture_sequences
star_capture_max_hamming
star_barcode_output_map
star_feature_search_mode
```

CAT-ATAC row values:

```text
star_barcode_read=R2
star_barcode_format=bc:8:23:-
star_umi_read=R1
star_umi_start=0
star_umi_length=12
star_feature_read=R3
star_capture_read=R1
star_capture_sequences=CAAGTTGATAACGGACTAGCC|CAAGTTGTAAACGGACTAGCC
star_capture_max_hamming=0
star_barcode_output_map=/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/atac2gex.tsv
star_feature_search_mode=free
```

If generic fields are too large for the first patch, a named layout such as
`star_layout=catatac_guide` is acceptable only as an initial bridge. The parser
must still expose enough provenance to catch the barcode offset and namespace.

### Phase 4: Namespace And MEX Output

The correct order is:

1. Extract raw ATAC barcode from guide R2 using `bc:8:23:-`.
2. Correct/rescue against the ATAC whitelist namespace, i.e. column 1 of
   `atac2gex.tsv` or `737K-arc-v1_atac.txt`.
3. Count UMIs per `(corrected_atac_barcode, guide)`.
4. When writing MEX barcodes, map corrected ATAC barcode to GEX barcode using
   `atac2gex.tsv` column 1 to column 2.
5. Merge guide MEX with GEX and ATAC surfaces in GEX/TRU barcode namespace.

Do not apply NXT/TRU middle-two-base translation to CAT-ATAC guide barcodes.

Implementation notes:

- `read_whiteList()` in `assignBarcodes.c` uses the first token/column for
  matching.
- `PfMultiMexStub::copyBarcodesTsv()` maps `barcodes.txt` through a two-column
  map when it is given that map path.
- `PfMultiProcess.cpp` currently enables `barcodeOutputMapPath` only for some
  NXT cases. Generalize this with an explicit per-library
  `star_barcode_output_map` or equivalent.

### Phase 5: Capture Sequence Filter

Add an optional prefilter in the split-read adapter:

```text
capture_read=R1
capture_sequences=CS1|CS2
capture_max_hamming=0
```

First pass should support exact match anywhere in R1. Hamming support can be
added later if needed, but the config should not hardcode only CS1.

Track and report:

```text
capture_total
capture_cs1_hits
capture_cs2_hits
capture_either_hits
capture_filtered_out
```

### Phase 6: Recipe Integration

The wrapper in this repo is `scripts/run_star_multiome_lane_smoke.sh`; the
canonical recipe may live in `/mnt/pikachu/morphic-recipes`. Before editing an
external recipe checkout, read that repo's local instructions.

Add an opt-in guide profile/flags with defaults off, for example:

```text
--guide-r1 FASTQ
--guide-r2 FASTQ
--guide-r3 FASTQ
--guide-feature-ref CSV_OR_TSV
--guide-whitelist-map TSV
--guide-layout catatac_guide
--enable-guide-arm
```

The generated pf-multi config should include GEX plus one CRISPR Guide Capture
library row for the guide FASTQs. ATAC continues through existing chromap flags.

For CAT-ATAC benchmark runs, keep:

```text
--dynamicThreadInterface 1
--crAssignConsumerThreads -1
--crAssignSearchThreads 1
--soloFeatures GeneFull
--outSAMtype None
```

unless BAM output is explicitly needed.

Use `--dynamicThreadTelemetry 1` for smoke/debug gates that need to prove
Chromap exercised the STAR permit hooks. Production benchmark runs can leave
telemetry off unless permit diagnostics are part of the measurement.

## Validation Gates

### Gate 1: Barcode Window Preflight

Run against the fixture before guide assignment:

```bash
python3 - <<'PY'
import gzip
from pathlib import Path

r2 = Path('/mnt/pikachu/catatac_gse288996/guide_redump/fixture/guide_R2.fastq.gz')
wl = Path('/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt')
trans = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
whitelist = {line.split()[0].upper() for line in wl.open() if line.split()}
checks = {(8, 16): 0, (9, 16): 0, (8, 24): 0}
n = 0
limit = 200000
with gzip.open(r2, 'rt') as f:
    while n < limit:
        h = f.readline()
        if not h:
            break
        seq = f.readline().strip()
        f.readline()
        f.readline()
        n += 1
        for off, length in checks:
            s = seq[off:off + length]
            if len(s) == length:
                checks[(off, length)] += s.translate(trans)[::-1].upper() in whitelist
print(f'n={n}')
for (off, length), hits in checks.items():
    print(f'offset={off} length={length} rc hits={hits} rate={hits / n if n else 0:.6f}')
PY
```

Expected on the staged fixture:

```text
offset=8 length=16: high, about 0.94 exact whitelist hits on first 200k reads
offset=9 length=16: zero or near zero
offset=8 length=24: zero or near zero
```

### Gate 2: Capture Filter Preflight

```bash
python3 - <<'PY'
import gzip
from pathlib import Path

r1 = Path('/mnt/pikachu/catatac_gse288996/guide_redump/fixture/guide_R1.fastq.gz')
cs1 = 'CAAGTTGATAACGGACTAGCC'
cs2 = 'CAAGTTGTAAACGGACTAGCC'
n = cs1_hits = cs2_hits = either_hits = 0
limit = 200000
with gzip.open(r1, 'rt') as f:
    while n < limit:
        h = f.readline()
        if not h:
            break
        seq = f.readline().strip()
        f.readline()
        f.readline()
        n += 1
        a = cs1 in seq
        b = cs2 in seq
        cs1_hits += a
        cs2_hits += b
        either_hits += a or b
print(f'n={n}')
print(f'cs1={cs1_hits} rate={cs1_hits / n if n else 0:.6f}')
print(f'cs2={cs2_hits} rate={cs2_hits / n if n else 0:.6f}')
print(f'either={either_hits} rate={either_hits / n if n else 0:.6f}')
PY
```

Expected on the staged fixture: CS2 contributes most of the captured reads, so
CS1-only filtering is invalid.

### Gate 3: Native Split Layout Equals Oracle

Create a tiny synthetic fixture and compare:

- materialized oracle `CB+UMI` barcode read plus R3 feature read;
- native split layout adapter with R2 CB, R1 UMI, R3 feature.

Required equality:

- `matrix.mtx`
- `features.tsv`
- `barcodes.tsv`
- per-run barcode/capture metrics, allowing only ordering differences where the
  existing process-features output is not order-stable.

### Gate 4: Real Fixture Smoke (process-features)

```bash
CATATAC_GUIDE_MAX_READS=1000 tests/test_catatac_guide_arm_smoke.sh
# 2M-fixture gate:
CATATAC_GUIDE_MAX_READS=200000 tests/test_catatac_guide_arm_smoke.sh
```

Oracle/native split-read smoke checks:

- duplicate feature warnings for `HIC2_1` and `HIC2_2`;
- 54 unique guide rows in output `features.tsv`;
- oracle/native `matrix.mtx` parity for the same read cap;
- GEX-space `barcodes.tsv` stub remap (`gex_tsv_whitelist_hits == barcode_count`).

Run the 2M-read fixture. Output must include:

- non-empty guide raw MEX;
- guide names from `ps_ref_crispri.tsv`;
- GEX-space `barcodes.tsv` after ATAC-to-GEX output mapping;
- non-zero UMI counts;
- run summary with barcode window, capture, guide-match, and namespace metrics.

Keep outputs under `tests/*_output*/` or `/tmp/`. Add any new persistent artifact
location to `tests/ARTIFACTS.md`.

### Gate 4c: Downsampled Trimodal STAR Smoke

One STAR invocation with GEX Solo + Chromap ATAC + guide pf-multi:

```bash
tests/test_catatac_trimodal_downsample_smoke.sh
# default: 100K read cap on GEX/guide, 100K staged ATAC triplets, guide fixture

CATATAC_TRIMODAL_MAX_READS=2000000 \
CATATAC_TRIMODAL_INLINE_ATAC_PEAK_MEX=yes \
tests/test_catatac_trimodal_downsample_smoke.sh
```

Outputs: `tests/catatac_trimodal_downsample_smoke_output/` (`trimodal_verify.json`).
The harness writes `RUN_STAR_TRIMODAL_SMOKE.sh`, `pf_multi_config.csv`, STAR logs,
and a compact verification report. Required checks cover GEX GeneFull output,
Chromap ATAC BAM/sidecar/summary, guide MEX in GEX barcode space, and raw merged
feature MEX. Inline ATAC peak MEX is checked when enabled. Cross-modality Jaccard
is **not** required on the independent downsamples. The smoke enables
`--dynamicThreadTelemetry 1`; with telemetry on, STAR treats zero ATAC permit
acquires after a successful Chromap run as a hook-integration failure.

Validation record:

- 2026-06-14: 100K trimodal smoke passed in
  `/mnt/pikachu/catatac_gse288996/full_bench/catatac_trimodal_smoke_100k_20260614T023346Z`
  with `CATATAC_TRIMODAL_MAX_READS=100000`, `CATATAC_TRIMODAL_THREADS=16`, and
  `CATATAC_TRIMODAL_INLINE_ATAC_PEAK_MEX=no`.
- Verifier summary: GEX filtered barcodes `11730`; guide MEX `54` features,
  `38738` barcodes, `55625` nnz; guide/GEX overlap `8951`; raw merged MEX
  `36655` features, `21394` barcodes, `92780` nnz; filtered merged MEX
  `36655` features, `11730` barcodes, `77749` nnz; Chromap ATAC BAM,
  binary sidecar, and summary were present and non-empty.
- Expected duplicate-name warnings for `HIC2_1` and `HIC2_2` were emitted during
  feature loading; the verifier confirmed `54` unique feature rows in the
  written guide MEX.

### Gate 4b: STAR pf-multi Smoke (guide arm)

Exercises `star_layout=catatac_guide` plus `star_barcode_output_map` through a
real STAR `alignReads` + pf-multi assign path (guide-only library row; no GEX
Solo merge required for this gate). The smoke uses `--soloType None` and omits
`--defaultCrCompat` / `--dynamicThreadInterface` because those bundles enable
Solo Droplet mode or can segfault on guide-only runs without concurrent GEX
mapping; full trimodal runs should still use the benchmark flag bundle from
Phase 6.

```bash
CATATAC_GUIDE_MAX_READS=1000 tests/test_catatac_guide_arm_pf_multi_smoke.sh
```

Expected checks (`tests/catatac_guide_pf_multi_verify.py`):

- duplicate feature warnings for `HIC2_1` and `HIC2_2` in STAR/native logs;
- pf-multi `features.tsv` row count is 54 with no duplicate feature names;
- pf-multi `barcodes.tsv` is GEX-space (`gex_whitelist_hits == barcode_count`,
  `atac_direct_gex_hits` near zero, `mapped_diff_from_atac > 0`);
- pf-multi guide count matrix matches the native split-read baseline for the
  same `CATATAC_GUIDE_MAX_READS` cap.

Outputs: `tests/catatac_guide_arm_pf_multi_smoke_output/` (`pf_multi_verify.json`,
`star/cr_assign/CRISPR_Guide_Capture/catatac_guide_fixture/`).

### Gate 5: Cross-Modality Barcode Overlap

For a trimodal run, compare guide cells against chromap ATAC cells after both
are in GEX barcode namespace:

```text
guide_cells_in_atac_cells
guide_cells_in_gex_cells
jaccard_guide_atac
jaccard_guide_gex
```

This is a cheap way to catch barcode offset, reverse-complement, and output-map
mistakes.

### Gate 6: Full Count Concordance

On the full guide library:

1. Compare guide-by-cell UMI count matrix against the reference stitched count
   surface if available.
2. Then compare per-cell guide assignments against `DMSO1_calls.rds`.

Do not use assignment disagreement alone to judge extraction correctness because
the deposited calls and STAR-suite `call_features` may use different caller
models and thresholds.

## Build And Test Hygiene

Before debugging any crash or regression:

```bash
make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR
```

Useful focused builds/tests after implementation:

```bash
make -C core/features/process_features tests/test_feature_duplicate_names tests/test_catatac_split_read assignBarcodes
CATATAC_GUIDE_MAX_READS=1000 tests/test_catatac_guide_arm_smoke.sh
CATATAC_GUIDE_MAX_READS=1000 tests/test_catatac_guide_arm_pf_multi_smoke.sh
make core
make feature-barcodes-tools
tests/run_cbub_regression_test.sh
tests/test_cr_compat_crispr_calling.sh
```

Add focused CAT-ATAC guide-arm smoke tests rather than relying on broad benchmark
scripts. Do not run benchmark jobs in parallel from this checkout.

## Common Failure Modes

- Off-by-one barcode extraction: near-zero ATAC whitelist hit rate and empty
  guide matrix.
- Literal 24 bp `substr(...,9,24)` interpretation: zero hits against 16 bp ARC
  whitelist/map.
- CS1-only filter: drops most captured fixture reads.
- NXT/TRU translation applied to ATAC guide barcodes: namespace corruption.
- Translating to GEX before correction: observed barcode cannot be corrected in
  the intended ATAC namespace.
- Treating all spacers as fixed 20 bp: loses 19 bp guides or creates biased
  matching.
- Forcing anchor search on an unanchored R3 spacer library: silent feature
  undercount.
- Comparing only final guide calls: caller variance can hide correct counting or
  make correct counting look wrong.
