# RUNBOOK: Velocyto Count Resolution on the Optimized Solo Path

Date: 2026-03-26
Branch context: `feature/flex-optimization-using-solo-20260325`

Related docs:
- `docs/HANDOFF_SCRNA_DOWNSTREAM_MEX_VELOCYTO_FINDINGS_20260313.md`
- `docs/HANDOFF_SOLO_OPTIMIZATION_20260324.md`

## Objective

Record how legacy STARsolo `Velocyto` counting actually works, why the current
non-Flex optimized bridge path does not emit `Velocyto`, and what the simplest
exact bridge-compatible design should look like.

This is a design/runbook document, not an implementation claim.

## Short Answer

Legacy `Velocyto` is **not** computed from final `Gene` or `GeneFull` matrices.

It is computed by:

1. annotating each read against transcript models
2. storing transcript-level compatibility evidence
3. deduping by `(CB, UMI)` while intersecting transcript support across reads
4. resolving the surviving transcript set to one gene if possible
5. collapsing the UMI to one of:
   - `spliced`
   - `unspliced`
   - `ambiguous`

Therefore the optimized bridge path can support exact `Velocyto`, but it must
preserve a separate transcript-level aggregate. `Gene` / `GeneFull` counts alone
are insufficient for exact parity.

## Legacy Data Flow

### 1. Per-read transcript annotation

`ReadAnnotations` is the per-read annotation container in:

- [ReadAnnotations.h](/mnt/pikachu/STAR-suite/core/legacy/source/ReadAnnotations.h)

Relevant fields:

- `trVelocytoType`
  - `vector<trTypeStruct>`
- `geneVelocytoSimple`
  - compact simplified velocyto summary

`trTypeStruct` is defined in:

- [SoloCommon.h](/mnt/pikachu/STAR-suite/core/legacy/source/SoloCommon.h)

```cpp
typedef struct{
    uint32 tr;
    uint8  type;
} trTypeStruct;
```

So each entry stores:

- `tr`: transcript id
- `type`: transcript-compatibility bitmask

The 4 transcript-compatibility classes are:

- [AlignVsTranscript.h](/mnt/pikachu/STAR-suite/core/legacy/source/AlignVsTranscript.h)

```cpp
enum {Intron=0, ExonIntron=1, ExonIntronSpan=2, Concordant=3, N=4};
```

This transcript-level evidence is produced in:

- [Transcriptome_classifyAlign.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/Transcriptome_classifyAlign.cpp)

Key lines:

- [Transcriptome_classifyAlign.cpp:254](/mnt/pikachu/STAR-suite/core/legacy/source/Transcriptome_classifyAlign.cpp#L254)
- [Transcriptome_classifyAlign.cpp:284](/mnt/pikachu/STAR-suite/core/legacy/source/Transcriptome_classifyAlign.cpp#L284)

Important behavior:

- for `ExonIntronSpan`, STAR also sets `Intron` and `Concordant` bits
- what survives is a categorical transcript-compatibility mask, not raw exon or
  intron base counts

### 2. Serialization into the Velocyto temp stream

For `--soloFeatures Velocyto`, STAR writes a dedicated per-read stream record in:

- [SoloReadFeature_record_base.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature_record_base.cpp)

Key lines:

- [SoloReadFeature_record_base.cpp:567](/mnt/pikachu/STAR-suite/core/legacy/source/SoloReadFeature_record_base.cpp#L567)

Serialized format:

- `iRead nTr tr1 type1 tr2 type2 ...`

This is transcript-level evidence, not gene-level output.

### 3. Velocyto UMI collapse

Legacy `Velocyto` counting is implemented in:

- [SoloFeature_countVelocyto.cpp](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countVelocyto.cpp)

Key structure:

- `vector<unordered_map<uintUMI, vector<trTypeStruct>>> cuTrTypes`
  at [SoloFeature_countVelocyto.cpp:19](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countVelocyto.cpp#L19)

The actual dedupe key is:

- `(CB, UMI)`

The payload is:

- the currently surviving transcript set for that UMI
- plus the merged 4-bit compatibility mask for each surviving transcript

### 4. Intersection semantics

When another read with the same `(CB, UMI)` is seen:

- transcript ids are intersected
- only transcripts present in both old and new support survive
- for surviving transcripts, type bits are ORed

This happens at:

- [SoloFeature_countVelocyto.cpp:68](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countVelocyto.cpp#L68)

This is the central point:

- the logic is **not** "sum spliced/unspliced/ambiguous per read"
- it is "intersect transcript support per `(CB, UMI)`, then classify once"

### 5. Resolution to gene and final class

After all reads for one `(CB, UMI)` are merged:

1. STAR checks whether all surviving transcripts map to the same gene
   - [SoloFeature_countVelocyto.cpp:110](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countVelocyto.cpp#L110)
2. if surviving transcripts map to multiple genes:
   - the UMI is dropped
   - it is **not** split across genes
   - [SoloFeature_countVelocyto.cpp:124](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countVelocyto.cpp#L124)
3. otherwise the merged transcript evidence is classified as:
   - `spliced`
   - `unspliced`
   - `ambiguous`
   - [SoloFeature_countVelocyto.cpp:127](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countVelocyto.cpp#L127)

## What the Current Optimized Bridge Path Does

The optimized non-Flex bridge path applies only to:

- `Gene`
- `GeneFull`
- `GeneFull_Ex50pAS`
- `GeneFull_ExonOverIntron`

See:

- [SoloFeature_countCBgeneUMI.cpp:24](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_countCBgeneUMI.cpp#L24)

`Velocyto` still uses the legacy stream-based path via:

- [SoloFeature_processRecords.cpp:106](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_processRecords.cpp#L106)

The bridge path does not carry forward `trVelocytoType` or an equivalent
transcript-level UMI aggregate.

As a result:

- optimized bridge `Gene*` counting works
- exact `Velocyto` cannot currently be emitted from the optimized bridge state

## Exact Bridge-Compatible Design

### Required information content

For exact legacy-style `Velocyto`, the minimal faithful representation is:

- key: `(CB, UMI, transcript)`
- value: `4-bit compatibility mask`

This is the simplest exact design because transcript identity is required for
the intersection step. A payload that stores only a 4-bit mask without
transcript identity is insufficient.

Equivalent alternative:

- key: `(CB, UMI)`
- value: compact aggregate of surviving transcript ids plus their 4-bit masks

This may be better for implementation efficiency, but the information content
must be equivalent.

### Online update rule

For each read:

1. derive `trVelocytoType` entries as legacy STAR does now
2. update the per-UMI transcript state
3. for repeated evidence on the same transcript, OR the 4-bit mask
4. for repeated reads with the same `(CB, UMI)`, preserve only transcripts that
   remain supported after intersection

### Finalization rule

For each `(CB, UMI)`:

1. gather surviving transcript entries
2. ensure they all map to one gene
3. if multigene, drop the UMI
4. collapse the merged masks to one final class:
   - `spliced`
   - `unspliced`
   - `ambiguous`
5. add `+1` UMI to that gene in the corresponding matrix

This reproduces the current legacy counting semantics.

## What Is Not Exact

The following approximation is simple but not exact:

- derive `spliced` from `Gene`
- derive `unspliced` from `GeneFull - Gene`
- infer or zero `ambiguous`

Why it is not exact:

- legacy `Velocyto` is driven by transcript-model compatibility, not final gene
  counts
- multigene UMIs are dropped upstream
- ambiguous classification depends on transcript evidence that is lost once only
  final gene matrices remain

This approximation may still be useful as an optional fast mode, but it should
not be presented as legacy-equivalent output.

## Cell Filtering / Output Notes

Legacy STAR uses `Gene` cell filtering to gate `Velocyto` kept barcodes:

- [SoloFeature_cellFiltering.cpp:16](/mnt/pikachu/STAR-suite/core/legacy/source/SoloFeature_cellFiltering.cpp#L16)

So the bridge-compatible implementation should preserve this behavior:

- perform `Velocyto` UMI classification independently
- apply the normal `Gene` cell filter to decide filtered-barcode output

No extra normalization to total gene counts is required.

## Recommended Implementation Order

1. Add a separate bridge-side velocyto aggregate with transcript identity.
2. Keep the existing optimized `Gene*` bridge path untouched.
3. Implement exact `(CB, UMI)`-level transcript intersection and final class
   assignment.
4. Reuse the existing `Velocyto` MEX/output conventions.
5. Validate against legacy `countVelocyto()` on:
   - 100K fixture
   - UCSF 2M fixture
   - full sample if needed

## Validation Targets

Compare new bridge-compatible `Velocyto` against legacy STARsolo `Velocyto` for:

- raw `spliced.mtx`
- raw `unspliced.mtx`
- raw `ambiguous.mtx`
- filtered counterparts
- filtered barcode set
- top-line `Summary.csv` / `Features.stats` where relevant

Initial acceptance criteria:

- exact or near-exact raw matrix parity on smoke fixtures
- exact filtered barcode parity after Gene filtering
- no regression to `Gene` / `GeneFull` optimized path wall time outside the new
  feature’s enabled surface

## Conclusion

Exact `Velocyto` support on the optimized Solo path is feasible.

The key design requirement is:

- preserve transcript identity plus 4-bit compatibility state through the UMI
  merge

The key non-requirement is:

- no post hoc normalization to `Gene` or `GeneFull` totals

The existing bridge path already has the right execution point to add this; it
simply does not yet preserve the needed transcript-level aggregate.
