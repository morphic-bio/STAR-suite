# No-BAM bridge Velocyto read-information gap

Discovered while validating the allocation audit. No fix for this issue is included
in the allocation changes. The failing run uses the preserved pre-pooling binary
`allocation_audit/ambiguity_build/STAR`, SHA256
`4d78a256899e8f6bbf933e203f4bb828370fb1ec0a94ff8fdda30f97d9207f6a`.

All paths below are relative to
`/home/lhhung/pf_larry_regression_20260911/`.

## Evidence

- `allocation_audit/velocyto_before`: 100,000 A375 GEX reads plus the full guide
  library; `GeneFull Velocyto`, `--soloInlineHashMode yes`,
  `STAR_SOLO_NONFLEX_HASH_BRIDGE=1`, `--outSAMtype None`. Exit status zero,
  GeneFull has 47,458 stored entries, but all three raw Velocyto matrices have
  zero entries. Exact argv, environment and input identities are in `execution.json`.
- `validation.json` marks that run unusable for allocation count parity.
- `allocation_audit/velocyto_before_legacy`: same input cap, no BAM, with
  `--soloInlineHashMode no` and the bridge environment variable **unset**.
  This yields 34,216 spliced, 6,813 unspliced and 5,196 ambiguous UMIs.
  The bridge switch tests variable presence, so setting it to `0` does not disable it.

The allocation comparison uses matched legacy-counting controls. Do not treat the
empty bridge matrices as a successful correctness control or compare these two
different counting configurations as a memory-performance experiment.

## Source trace

- `ParametersSolo.cpp`: the Velocyto configuration sets the selected gene-like
  feature's `readInfoYes`, but `trackReadIdsForTags` is enabled by BAM CB/UB tags.
- `SoloFeature_countCBgeneUMI.cpp`: `needPackedReadInfo` requires
  `readIndexYes && (trackReadIdsForTags || (readInfoYes && !nonFlexBridgePath))`.
  The no-BAM direct bridge therefore skips this allocation despite Velocyto's request.
- `SoloReadFeature.cpp`: the inline read-ID tracker is created only when
  `trackReadIdsForTags` is enabled.
- `SoloFeature_countVelocytoBridge.cpp`: the selected gene-like feature supplies
  packed CB/UMI/status. Records with status other than one are skipped before UMI merging.

## Next isolated correction

Make Velocyto's per-read information requirement reach the direct bridge even
without BAM tag output. Preserve the per-feature `readIndexYes`/stride guards:
forcing tracking on a stride-two stream previously caused invalid read IDs.
Validate a nonempty, no-BAM bridge result and its GEX/guide outputs on a changed
binary before full-scale benchmarking. Keep this correction separate from the
allocation-only changes. Existing executions are preserved; this handoff does
not authorize identical reruns.
