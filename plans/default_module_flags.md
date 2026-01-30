# Default Module Flag Groups

Goal: define `--default*` flag bundles that simplify common workflows.
These bundles expand to predefined parameter sets while allowing explicit overrides.

**Implemented CLI flags** (CamelCase, as implemented in STAR):
- `--defaultBulk Yes`
- `--defaultCoreScrna Yes`
- `--defaultCrCompat Yes`
- `--defaultFlex Yes`
- `--defaultSlam Yes`
- `--defaultVbem Yes`
- `--defaultA375Parity Yes`

Policy: include **all non‑default** parameters that are **not** dataset‑ or
hardware‑specific (e.g., no FASTQs, references, prefixes, threads, whitelists).
User-specified parameters always override default group values.

## Smoke tests for default bundles (proposed)

Goal: ensure each `--default-*` group runs end‑to‑end on a small fixture.

Suggested coverage (existing scripts to wire to `--default-*`):
- **Bulk**: `tests/run_baseline_test.sh`
- **scRNA**: `tests/run_solo_smoke.sh`
- **CR‑compat**: `tests/run_cr_compat_integration_smoke.sh` (or `tests/test_cr_compat_crispr_calling.sh`)
- **Flex**: `tests/flex_smoke/run_flex_smoke.sh` (if fixtures available)
- **CB/UB**: `tests/run_unsorted_cbub_smoke_test.sh`
- **SLAM**: `tests/run_slam_*` *(choose the smallest fixture)*
- **VBEM**: `tests/transcriptvb/*` *(vbem smoke/parity)*
- **Feature tools**: `core/features/process_features/tests/*`

TODO:
- Add a small wrapper that injects `--default-*` into each smoke script and
  confirms overrides for dataset‑specific parameters still work.
- Record fixture locations in `tests/ARTIFACTS.md` if not already listed.

## Bulk (permissive, VB‑quant friendly)

### `--defaultBulk Yes`
Use for bulk RNA‑seq when you want **permissive** alignment filters (appropriate
for Salmon/VB quant and parity runs).

Suggested bundle (edit):
- `--runThreadN <N>` *(auto‑detect in wrapper; hardware‑specific)*
- `--outSAMtype BAM SortedByCoordinate` *(binary output; standard downstream input)*
- `--outSAMattributes NH` *(trim tags; smaller outputs for parity)*
- `--outSAMprimaryFlag AllBestScore` *(keep all best alignments)*
- `--outSAMmultNmax 10000` *(cap multi‑mapper outputs)*
- `--outFilterMultimapNmax 10000` *(retain multimappers for VB/quant)*
- `--outFilterMultimapScoreRange 4`
- `--outFilterMismatchNmax 6` *(tighter than default 10; still permissive)*
- `--outFilterMatchNmin 25` *(minimum matched length for stability)*
- `--outFilterMatchNminOverLread 0` *(no fraction cutoff)*
- `--outFilterScoreMinOverLread 0` *(no fraction cutoff)*
- `--alignIntronMax 500000` *(allow spliced alignments)*
- `--alignMatesGapMax 1000` *(PE insert sizes; omit for SE‑only workflows)*
- `--winAnchorMultimapNmax 200` *(more permissive anchor multimapping)*
- `--outBAMcompression 6` *(smaller BAMs)*

Optional if standardized:
- `--outBAMsortMethod samtools` *(only if we standardize on samtools sort)*

Related scripts (bulk):
- `tests/run_baseline_test.sh`
- `tests/run_pe_preprod_validation.sh`
- `tests/run_ychrom_*`
- `tests/run_y_removal_comprehensive_test.sh`

Notes:
- Some no‑intron tests use `--alignIntronMax 1` (not part of this permissive set).
- `--alignEndsType Local` is STAR default; keep unless you explicitly want EndToEnd.

Close calls / mixed usage (review before locking):
- `--outFilterMultimapNmax 10000` vs `10` (strict vs permissive)
- `--outFilterMismatchNmax 6` vs `0` (stringency)
- `--outFilterMatchNmin 25` vs `1` (stringency)
- `--alignEndsType EndToEnd` vs `Local` (default)
- `--alignIntronMax 500000` vs `1` (spliced vs no‑intron tests)

## scRNA (non‑Flex)

### `--defaultCoreScrna Yes`
Use for standard scRNA‑seq (GEX) without Flex.

Bundle includes:
- `--soloType Droplet` *(enable droplet pipeline)*
- `--soloFeatures Gene` *(explicit GEX)*
- `--soloCellFilter EmptyDrops_CR` *(modern cell calling; CR defaults)*
- `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts` *(robust CB matching)*
- `--soloUMIdedup 1MM_CR` *(match perturb defaults)*
- `--soloUMIfiltering MultiGeneUMI_CR` *(match perturb defaults)*
- `--soloMultiMappers Rescue` *(match perturb defaults)*
- `--soloCbUbRequireTogether no` *(match perturb defaults; CB emitted if UB invalid)*

Protocol‑specific (still required, but dataset‑specific):
- `--soloCBstart/--soloCBlen`
- `--soloUMIstart/--soloUMIlen`
- `--soloCBwhitelist <path>`

Related scripts:
- `tests/run_solo_smoke.sh`
- `tests/run_a375_cr_like_gex.sh`
- `tests/run_cr_parity_100k.sh`

## CR‑compat / CRISPR

### `--defaultCrCompat Yes`
CR‑compat mode with CRISPR calling.

Bundle includes:
- `--soloType Droplet`
- `--soloFeatures Gene GeneFull` *(CR parity; introns included)*
- `--soloCellFilter EmptyDrops_CR`
- `--soloCrMode CR`
- `--soloCrGexFeature GeneFull`
- `--crMinUmi 10` *(configurable)*
- `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`
- `--soloUMIdedup 1MM_CR`
- `--soloUMIfiltering MultiGeneUMI_CR`
- `--soloMultiMappers Rescue`
- `--soloCbUbRequireTogether no`

Note: `--soloKeysCompat cr` is NOT set by default (requires `--soloProbeList`).

Protocol‑specific (still required, but dataset‑specific):
- `--soloCBstart/--soloCBlen`
- `--soloUMIstart/--soloUMIlen`
- `--soloCBwhitelist <path>`

Related scripts:
- `tests/test_cr_compat_crispr_calling.sh`
- `tests/run_a375_gex_features_cr_parity.sh`
- `tests/run_a375_gex_features_cr_parity_genefull.sh`
- `tests/run_cr_compat_integration_smoke.sh`

## Flex (separate pipeline)

### `--defaultFlex Yes`
STAR‑Flex standard defaults. **Do not reuse for non‑Flex pipelines.**

Bundle includes:
- `--flex Yes` *(enable Flex mode)*
- `--soloFeatures Gene` *(plus Flex probe features)*
- `--soloCellFilter EmptyDrops_CR`
- `--soloRunFlexFilter yes` *(enable Flex filter pipeline)*
- `--soloRemoveDeprecated Yes`
- `--removeDeprecated Yes` *(index build)*
- `--soloSampleSearchNearby no` *(strict demux)*
- `--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts`
- `--soloUMIdedup 1MM_CR`
- `--soloUMIfiltering MultiGeneUMI_CR`
- `--soloMultiMappers Rescue`

Flex alignment defaults (part of the bundle):
- `--outFilterMultimapNmax 10000`
- `--outFilterMultimapScoreRange 4`
- `--outFilterMismatchNmax 6`
- `--outFilterMatchNmin 25`
- `--outFilterMatchNminOverLread 0`
- `--outFilterScoreMinOverLread 0`
- `--alignIntronMax 500000`
- `--winAnchorMultimapNmax 200`
- `--outSAMtype BAM Unsorted` *(unsorted by default)*
- `--outSAMattributes NH`
- `--outSAMprimaryFlag AllBestScore`
- `--outSAMmultNmax 10000`
- `--outBAMcompression 6`
- `--outBAMsortMethod samtools` *(only if SortedByCoordinate is requested)*

Protocol‑specific (still required, but dataset‑specific):
- `--soloCBstart/--soloCBlen`
- `--soloUMIstart/--soloUMIlen`
- `--soloCBwhitelist <path>`
- `--soloSampleWhitelist <path>`
- `--soloSampleProbes <path>`
- `--soloProbeList <path>`

Related scripts:
- `tests/flex_smoke/run_flex_smoke.sh`
- `tests/run_flex_with_tags_test.sh`
- `tests/run_flex_inline_test.sh`
- `tests/run_flex_multisample_test.sh`
- `tests/run_flex_cbub_validation_test.sh`
- `flex/run_all.sh`
- `flex/run_SC2300771.sh`
- `flex/run_SC2300772.sh`

## CB/UB Tagging and BAM injection (non‑Flex)

CB/UB is **not** a separate default group. It should follow the **parent
pipeline** defaults (scRNA or Flex), plus the Solo CR‑style UMI settings
when needed.

CB/UB related scripts:
- `tests/run_cbub_regression_test.sh`
- `tests/run_sorted_bam_cbub_test.sh`
- `tests/run_sorted_bam_cbub_nonflex_test.sh`
- `tests/run_unsorted_bam_cbub_test.sh`
- `tests/run_unsorted_cbub_smoke_test.sh`
- `tests/run_unsorted_cbub_a375_nonflex.sh`
- `tests/run_cbub_tag_validation_test.sh`
- `tests/compare_cb_ub_tags.sh`

## SLAM‑seq (bulk superset)

### `--defaultSlam Yes`
SLAM is a **superset of bulk**: applies all `--defaultBulk` parameters first,
then adds SLAM‑specific flags.

SLAM-specific add‑ons:
- `--slamQuantMode 1`
- `--slamSnpMaskModel em`
- `--slamSnpMaskOnly 1`
- `--slamWeightMode dump`

Related scripts:
- `tests/run_slam_*`
- `tests/slam/*`

## VBEM / TranscriptVB

### `--defaultVbem Yes`
VBEM quant/compat workflows.

Bundle includes:
- `--quantMode TranscriptVB` *(enable VBEM path)*
- `--quantVBem 1` *(enable EM)*
- `--quantVBgcBias 1` *(enable GC bias)*
- `--quantVBgenes 1` *(keep gene‑level VB; also keep transcript‑level outputs)*
- `--quantVBgenesMode Tximport` *(compat output)*
- `--quantVBprior 0.01` *(default; prior sweep shows 0.0001–0.1 tie, 0.01 acceptable)*

Related scripts:
- `core/features/vbem/test/*`
- `core/features/vbem/tools/*`
- `tests/transcriptvb/*`

## Feature Tools (assignBarcodes / process_features)

Feature tools have their own CLI (not STAR flags). No separate `--default-*`
flag was implemented because the tools already have sensible defaults:
- Auto-detect offset from pattern column
- Warn (not error) on heterogeneous offsets
- Use `--strict-offset-check` to hard-fail on heterogeneity

Related scripts:
- `core/features/process_features/tests/*`
- `core/features/process_features/scripts/*`
- `core/features/feature_barcodes/*`

## GEX + Feature Parity (A375)

### `--defaultA375Parity Yes`
Defaults used for A375 GEX + feature parity runs.

Bundle includes:
- `--soloCrMode CR`
- `--soloCrGexFeature GeneFull`
- `--crMinUmi 10`

Related scripts:
- `tests/run_a375_gex_features_cr_parity.sh`
- `tests/run_a375_gex_features_cr_parity_genefull.sh`
