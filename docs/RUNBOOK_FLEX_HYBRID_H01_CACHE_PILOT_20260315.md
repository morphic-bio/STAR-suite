# Runbook: Flex Hybrid H0/H1 Cache Pilot

Date: 2026-03-15
Status: design ready, not yet implemented
Scope: 100-probe pilot on the 100K `SC2300771` Flex fixture

## Objective

Test whether a precomputed hybrid-reference cache can short-circuit Flex probe
assignment with zero output drift.

The pilot is intentionally narrow:

1. choose 100 pilot probes
2. precompute all `H0/H1` 50-mers for those probes
3. align those synthetic 50-mers to the current hybrid Flex reference
4. resolve them with the current Flex decision logic
5. benchmark on the 100K fixture with this runtime policy:
   - cache value `0` -> keep
   - cache value `1` -> keep
   - cache value `2` -> deny
   - cache miss -> fall through to normal alignment
6. compare outputs

Expected result for the pilot: exact byte match. If this fails, the first
question is whether `value=2 -> deny` is too aggressive.

## Core Hypothesis

If a 50-mer is within `H0/H1` of one of the pilot probes, and if its resolved
outcome against the **same hybrid STAR reference** is precomputed with the
**same current Flex resolver semantics**, then the runtime path can reuse that
decision instead of aligning the read again.

This is not probe-space uniqueness. This is a cached answer from the current
hybrid-reference decision surface.

## Non-Negotiable Rules

1. The cache builder must use the same hybrid reference currently used by Flex:
   - `/storage/flex_filtered_reference/star_index`
2. The cache builder must use the same current Flex decision rules as:
   - `core/legacy/source/flex/SoloReadFeature_record_flex.cpp`
   - `flex/source/GeneResolver.cpp`
3. The pilot must not silently change behavior for reads outside the 100-probe
   cache surface.
4. Any parameter that changes alignment or resolver behavior invalidates the
   cache and requires a rebuild.
5. Before any pilot benchmark, do a clean rebuild:

```bash
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
```

## Why Start With 100 Probes

The goal is not speed at this stage. The goal is falsifiability.

A 100-probe pilot is enough to answer the important questions:

1. Can `value 0/1` be treated as a strict fast-path keep?
2. Can `value 2` be treated as a strict deny without later rescue?
3. Does the pilot produce exact output parity on a real Flex dataset?

If the logic is wrong for 100 probes, it is wrong for the full probe set.

## Inputs

### Reference inputs

- Hybrid STAR index:
  `/storage/flex_filtered_reference/star_index`
- Filtered probe CSV:
  `/storage/flex_filtered_reference/filtered_reference/filtered_probe_set.csv`
- Probe gene list:
  `/storage/flex_filtered_reference/filtered_reference/probe_list.txt`

### 100K benchmark fixture

- R2 FASTQs:
  `/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L00{1..8}_R2_001.fastq.gz`
- R1 FASTQs:
  `/storage/downsampled_100K/SC2300771/SC2300771_GT23-14630_GATAATACCG-TTTACGTGGT_S5_L00{1..8}_R1_001.fastq.gz`
- Sample whitelist:
  `/storage/SC2300771_filtered_2M/sample_whitelist.tsv`
- Sample probe variants:
  `/mnt/pikachu/JAX_scRNAseq01_processed/probe-barcodes-fixed-rna-profiling-rna.txt`

### Existing reference benchmark surface

March 11 direct-scan artifacts and the matching STAR benchmark already exist
for this fixture. Reuse them for context, but do not treat them as proof of
parity for this cache experiment.

## Artifact Layout

Use a new untracked artifact root for each pilot run:

```text
/storage/100K/tmp/flex_h01_cache_pilot_YYYYMMDD_HHMMSS/
  baseline/
  pilot/
  selection/
  cache/
  comparison/
```

Recommended contents:

- `selection/pilot_probes.csv`
- `selection/pilot_genes.txt`
- `selection/selection_manifest.tsv`
- `cache/h01_cache.tsv`
- `cache/cache_summary.tsv`
- `pilot/cache_hit_summary.tsv`
- `comparison/whole_mex_diff.txt`
- `comparison/pilot_gene_subset_diff.txt`

If this experiment becomes repeatable infrastructure, add the artifact root to
`tests/ARTIFACTS.md`.

## Phase 0: Baseline Capture

Use the current Flex command line for the 100K fixture and keep the exact
baseline artifacts.

Minimum baseline outputs:

1. `Log.out`
2. `Log.final.out`
3. `Solo.out/Gene/raw/{matrix.mtx,barcodes.tsv,features.tsv}`
4. `per_sample/*/Gene/filtered/{matrix.mtx,barcodes.tsv,features.tsv}`
5. wall time and RSS

Whole-output comparison can reuse:

```bash
tests/compare_mex_outputs.sh <baseline_dir> <pilot_dir>
```

## Phase 1: Choose the 100-Probe Pilot Set

Do not choose an arbitrary silent block of probes. The pilot needs to be
exercised by the 100K fixture.

Recommended selection procedure:

1. Use the baseline 100K raw Gene MEX to rank genes by total count.
2. Walk that ranked gene list in descending count order.
3. For each gene with at least one probe in `filtered_probe_set.csv`, select
   one representative probe record.
4. Stop at 100 probe records.

Selection policy for reproducibility:

1. stable sort candidate probes by `gene_id`, then `probe_id`
2. within each chosen gene, take the first probe record
3. write the final probe records to `selection/pilot_probes.csv`
4. write the touched gene IDs to `selection/pilot_genes.txt`
5. write `selection/selection_manifest.tsv` with:
   - `gene_id`
   - `probe_id`
   - `probe_seq`
   - baseline gene count
   - selection rank

This gives exactly 100 probes and a well-defined affected gene set.

## Phase 2: Offline Cache Build

Build a pilot cache for only these 100 probes.

### 2A. Enumerate `H0/H1`

For each 50 bp pilot probe:

1. emit the exact probe sequence as `H0`
2. emit all single-substitution neighbors as `H1`
3. deduplicate identical 50-mers across probes

Expected upper bound before dedup:

- `100 * (1 + 50 * 3) = 15,100` synthetic 50-mers

### 2B. Align Against the Hybrid Reference

For each unique synthetic 50-mer:

1. align against `/storage/flex_filtered_reference/star_index`
2. use the same alignment settings that matter for current Flex probe
   resolution
3. feed the resulting candidates through the same current resolver logic

Important:

- this is against the **hybrid** reference, not probe-only space
- this must preserve current probe/genomic competition
- do not replace the resolver with a simplified "best probe wins" rule

### 2C. Encode the Pilot Cache

Write one row per synthetic 50-mer to `cache/h01_cache.tsv` with at least:

- `sequence_50mer`
- `pilot_probe_id`
- `pilot_gene_id`
- `best_distance`
- `resolved_gene_id`
- `cache_value`
- `reason`

Proposed encoding for the pilot:

- `0` = unique keep from exact (`H0`)
- `1` = unique keep from one-mismatch (`H1`)
- `2` = ambiguous or problematic under the hybrid/reference resolver surface

For the pilot, `value=2` intentionally means "deny" at runtime. This is the
hypothesis under test.

`reason` should distinguish at least:

- `UNIQUE_H0`
- `UNIQUE_H1`
- `PROBE_AMBIG`
- `GENOMIC_AMBIG`
- `MIXED_CONFLICT`
- `NONCANONICAL`
- `NO_CANDIDATES`
- `OTHER`

Also write `cache/cache_summary.tsv` with counts by:

- distance
- cache value
- reason
- probe
- gene

## Phase 3: Runtime Pilot Policy

Add a runtime pilot path that only consults the cache for the selected 100
probes.

### Pilot lookup policy

At the current expected Flex probe window:

1. extract the runtime 50-mer candidate
2. look up that 50-mer in `h01_cache.tsv`
3. if `cache_value=0` or `1`, keep the cached resolved gene and skip
   alignment for that read
4. if `cache_value=2`, deny this read for the pilot policy
5. if no cache entry exists, run the normal alignment path unchanged

This is intentionally strict because the experiment is about whether `2`
really propagates as a terminal negative.

### Required runtime telemetry

Write `pilot/cache_hit_summary.tsv` with at least:

- `cache_hit_0_reads`
- `cache_hit_1_reads`
- `cache_hit_2_reads`
- `cache_miss_reads`
- `cache_hit_unique_genes`
- `cache_hit_unique_probes`

Also write optional read-level traces for debugging:

- qname
- cache value
- resolved gene
- fallback used yes/no

## Phase 4: 100K Benchmark Run

Run the pilot on the same 100K fixture and with the same command shape as the
baseline, changing only the pilot-cache flags.

Required outputs:

1. full raw Gene MEX
2. full per-sample filtered MEX
3. wall time and RSS
4. cache hit telemetry

## Phase 5: Comparison

### Primary acceptance criterion

The correct result is exact whole-output parity:

1. byte-identical pooled raw MEX
2. byte-identical per-sample filtered MEX

Use:

```bash
tests/compare_mex_outputs.sh <baseline_dir> <pilot_dir>
```

### Secondary localization step

If whole-output parity fails, localize the diff to the touched pilot genes.

Minimum required subset comparison:

1. extract rows in `features.tsv` whose gene IDs are in `selection/pilot_genes.txt`
2. extract matching `matrix.mtx` rows
3. diff the extracted subset between baseline and pilot

This subset comparison is diagnostic only. It does not replace whole-output
parity.

## Pass / Fail Interpretation

### PASS

All of the following are true:

1. whole-output MEX parity is exact
2. the pilot touched a non-trivial number of reads
3. the pilot touched a non-trivial number of the 100 probes
4. runtime telemetry shows the cache actually exercised `value 0/1/2`

If this passes, the idea is strong enough to scale.

### FAIL MODE A: `value 0/1` causes drift

Interpretation:

- the cached positive path is not equivalent to current Flex
- full-read context still matters even after hybrid pre-alignment

Action:

- stop scale-up
- inspect read-level differences before changing policy

### FAIL MODE B: only `value 2` causes drift

Interpretation:

- `value=2 -> deny` is too aggressive
- negative propagation is not universally safe

Action:

1. rerun the pilot with this modified policy:
   - `0/1` keep
   - `2` fallback-align
   - miss fallback-align
2. if parity is restored, keep `2` as a fallback class rather than a deny class

This is the most likely failure mode and should be tested immediately if the
strict pilot drifts.

### FAIL MODE C: pilot has too little traffic

Interpretation:

- the selected probes were not exercised enough to say anything useful

Action:

- rebuild the 100-probe pilot set from more highly expressed baseline genes

## Success Criteria for Scaling Beyond 100 Probes

Do not scale to the full probe set until the 100-probe pilot shows:

1. exact whole-output parity on the 100K fixture
2. enough cache-hit traffic to make the result meaningful
3. no unresolved drift from `value 0/1`

If the strict `2 -> deny` policy fails but `2 -> fallback` passes, the design
can still scale. In that case the cache becomes:

- `0/1` = terminal fast-path keep
- `2` = non-terminal cached negative that forces fallback alignment
- miss = fallback alignment

That still delivers a useful filter and should remain compatible.

## Recommended First Implementation Notes

Keep the pilot implementation minimal:

1. no full-probe-set cache yet
2. no memory optimization yet
3. no clever compression yet
4. no attempt to optimize the builder before parity is established

The first version only needs to answer one question:

Can a 100-probe hybrid `H0/H1` cache preserve exact output parity on the 100K
Flex fixture?

If yes, then optimization and scale-up are worth doing. If no, the failure
mode will tell us whether the problem is:

- positive-cache equivalence
- negative-cache equivalence
- or just pilot selection quality
