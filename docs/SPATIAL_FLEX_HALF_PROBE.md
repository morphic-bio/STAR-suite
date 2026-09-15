# Spatial Flex half-probe processing

The spatial half-probe migration replaces the alignment-backed spatial Flex
adapter. `--soloSpatialFlexIntegrated yes --flex yes` now uses the fully fused
FASTQ reader and shared H0/H1X2 classifier, and completes the existing spatial
molecule pipeline without loading a genome. Visium HD GEX retains its normal
GeneFull alignment path. Ordinary Chromium Flex legacy reproduction is unchanged.

Supply the matched half-probe cache and ordered gene list explicitly. Generate
the cache once with `--runMode hashCacheGenerate --hashCacheTiers H0,H1X2`,
then optionally pack it as described in [FLEX_KHASH_CACHE.md](FLEX_KHASH_CACHE.md).
Generation may use the probe reference; production processing does not.
Spatial input is paired FASTQ ordered R2, raw R1. Plain text, gzip and native
BGZF are supported. Packed Chromium CBQ cannot represent this raw-R1 contract
and is rejected. An external decompressor such as rapidgzip can be selected
for ordinary gzip with `--readFilesCommand`.

Keep the existing spatial barcode contract, oligo files and capacity settings.
Use `--soloType CB_UMI_Complex --soloFeatures Gene --soloSkipProcessing yes
--soloRunFlexFilter no --soloStrand Unstranded --soloMultiMappers Rescue`,
`--outSAMtype None --outSAMattributes None --outSJtype None --chimSegmentMin 0`,
and the documented spatial UMI/filter options. The normal defaults select
`--flexPipeline yes --flexPipelineNTriage 0 --flexPipelineNSolo 0 --flexNoAlign 1`.
A failed no-genome activation guard is fatal for spatial Flex; it never falls
back to reference loading. Spatial `--flexLegacy yes` is no longer supported.
No BAM, feature sidecar, candidate TSV or alignment-reader tap is produced.

Each fused worker owns a raw-R1 decode and exactly one terminal feature
completion. All valid coordinates, their likelihoods and the raw UMI reach the
spatial engine; ordinary single-barcode correction and Solo count stores are
bypassed. Candidate-specific UMI correction, spatial family resolution, bounded
spill and all four assignment products at 2/8/16 um are reused.

Source ordinals are `record_in_lane * number_of_lanes + lane_index`, all
zero-based. This identity is stable across thread schedules, gzip and BGZF.
The compact uint32 limit is checked before multiplication or truncation. For
uneven lanes this can be reached before the total-read capacity, and the run
then fails explicitly. It is not a contiguous global input counter.

The run summary reports H0, H1 (including conservative single-N decisions),
H1X2, deny and terminal miss counts. Assigned reads equal the sum of the three
keep classes; unassigned reads equal deny plus miss. The feature/spatial
cross-tab covers every decoded R1. The historical alignment counters remain
zero and `feature_route=half_probe_no_alignment` identifies the new semantics.
Never treat these misses as alignment-fallback attempts.

`tests/test_spatial_flex_half_probe.py` exercises the production binary on
small real paired inputs with an empty genome directory. It compares every MEX
component across plain/gzip/BGZF, one/four workers and forced spill, and checks
rejection of malformed mates, unsupported legacy routing and disabled
no-genome routing. Use a fresh output directory and the matched spatial cache,
feature list and barcode inputs. It records commands, binary hashes and outcomes.

The campaign runbook, including CRC, both SPATCH reference arms, ovarian GEX
and both CODEX analyses, is owned by `visium-hd-processing-recipes`:
`docs/RUNBOOK_VISIUM_HD_FLEX_HALF_PROBE_MIGRATION_20260915.md`.
Implementation tests are development evidence; paper runs require the accepted
public release containing this migration.
