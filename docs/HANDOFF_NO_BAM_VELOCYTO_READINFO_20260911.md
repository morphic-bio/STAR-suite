# No-BAM bridge Velocyto read-information correction

Implemented as a separate correctness change after the allocation pass (`5a194e7`).
The direct bridge now provides corrected read identities to Velocyto with
`--outSAMtype None`; BAM and tag-table output are not required.

## Defect and correction

Velocyto requests `readInfoYes` on its selected gene-like feature, but the direct
bridge previously skipped packed read storage and discarded individual read IDs.
The existing BAM tracker cannot simply be enabled: direct bridge collapse explicitly
rejects that configuration and does not populate that tracker.

The correction retains flat read-identity records only on the selected source
feature. It uses the same resolved ambiguous-CB decisions as the count hash, applies
its one-exact-match gate, and replays the final gene-specific corrected UMI or
rejected-UMI status after MultiGeneUMI_CR. Other feature streams retain their existing
read-index/stride behavior; GEX-only runs do not allocate the new read buffers.
The gate uses the selected Velocyto feature explicitly: the existing default
tag-table setup also sets `readInfoYes` when no tag file is emitted, so testing
`readInfoYes` alone would retain unnecessary GEX-only read buffers.

Read eligibility follows the current `CountingSink::onRecord`: only reads assigned
to the selected gene feature enter replay. A transcript overlap alone does not
admit a no-feature or multi-gene read. This preserves existing counting policy.
The bridge's established key-level ambiguous barcode resolution remains unchanged.

Velocyto now reports an error if its source lacks packed read storage instead of
silently producing empty layers. Hash-snapshot replay rejects requests for per-read
information because those snapshots do not contain read identities.

## Validation

Artifact root: `/home/lhhung/pf_larry_regression_20260911/velocyto_readinfo_fix/`.
Each run uses 100,000 A375 GEX reads plus the full guide library, 32 threads,
auto BGZF detection, and no BAM or tag sidecars. These are correctness fixtures,
not full-set paper benchmarks.

Final clean build: `build4/STAR`, SHA256
`6c5bd5c138482ceb6b153dada5e9cd52e9c4f44a1c132f952535a342ffff9994`.
The build manifest records base revision `5a194e7`, the complete source patch,
clean/build commands and completion. Test wrappers record inputs, argv, environment
and binary/driver hashes.

| Configuration | Spliced UMIs | Unspliced UMIs | Ambiguous UMIs |
| --- | ---: | ---: | ---: |
| Preserved defective standard bridge | 0 | 0 | 0 |
| Fixed standard bridge (`bridge_final`) | 34,091 | 6,784 | 5,177 |
| Exact-barcode legacy (`exact_legacy`) | 33,617 | 6,668 | 5,096 |
| Exact-barcode bridge (`exact_bridge_fixed`) | 33,617 | 6,668 | 5,096 |

- `restored_fixed_comparison.json`: the fixed standard bridge preserves all six
  GEX/guide MEX surfaces and three guide-call tables against the defective run.
  All raw and filtered Velocyto layers are now nonempty. Filtered layer sums are
  31,897 / 6,442 / 4,815.
- `exact_fixed_comparison.json`: all eight matrix surfaces, three guide tables,
  six raw/filtered layer matrices and four barcode/gene axis files match the
  existing legacy counting path. Exact matching isolates UMI replay from the
  two paths' different ambiguous-barcode resolution policies.
- `bucket_comparison.json`: standard stream replay and CB-bucket replay match
  on all eight matrix surfaces, three guide tables, six layer matrices and four
  axis files. Both use `build3`, before the final allocation gate; the final-build
  comparison below verifies unchanged layer results.
- `final_bridge_comparison.json`: the final allocation gate preserves all eight
  matrix surfaces, three guide tables, six layers and four axes from `build3`.
- `gex_only_comparison.json`: the final GEX-only control retains no replay state
  or packed read array and preserves all six GEX/guide matrices and three guide
  tables against `bridge_final`.
- `multi_feature_final_comparison.json` and `multi_feature_allocation_check.json`:
  `Gene GeneFull Velocyto` retains identities only for Gene, allocates one packed
  array, and preserves all ten matrices, three guide tables, six layers and four
  axes against the earlier mixed-feature control. GeneFull remains a stride-two
  stream in this final configuration.
- The standard bridge retained 61,040 assigned gene reads after CB resolution:
  60,963 receive valid corrected identities and 77 retain rejected-UMI status;
  264 accepted original UMI keys underwent correction.
- `unit_fixed`: AddressSanitizer/UndefinedBehaviorSanitizer test passes for repeated
  reads, gene-specific UMI correction/rejection, resolved/unresolved ambiguous CBs,
  barcode gating and exclusion of reads without a selected gene assignment.

Regression tools: `tests/run_velocyto_allocation_probe.py`,
`tests/compare_velocyto_readinfo_outputs.py` and `tests/test_bridge_read_info.cpp`.
The wrapper checks actual layer UMI sums, because stored Matrix Market rows may
include explicit zeros.

## Preserved investigation evidence

The original empty-layer run is
`../allocation_audit/velocyto_before`; its pre-pooling binary SHA256 is
`4d78a256899e8f6bbf933e203f4bb828370fb1ec0a94ff8fdda30f97d9207f6a`.
It is excluded as an allocation correctness/performance control by `validation.json`.
Its GEX/guide outputs remain useful unchanged-count controls for this repair.

The first repair draft admitted no-feature reads. Although it preserved GEX counts,
`exact_comparison.json` caught a 54-UMI Velocyto difference against the legacy path.
That unintended policy change was removed. `bridge` and `exact_bridge` are marked
superseded in their validation records; use the `_fixed` runs above.

Full-scale Velocyto runtime and memory remain unmeasured. Per-read state is required
when requesting Velocyto and adds storage relative to GEX-only direct hashing.

## Existing source-feature policy

If both Gene and GeneFull are requested, Velocyto uses Gene by the existing
priority rule. On this fixture that produces 33,635 spliced, zero unspliced and
5,076 ambiguous UMIs: the current CountingSink eligibility excludes reads with
no Gene assignment, including purely intronic reads. This correction preserves
that policy; it does not repair or redefine Gene-only intronic eligibility.
The validated intron-inclusive configuration here is `GeneFull Velocyto`, with
6,784 unspliced UMIs. Any change to Gene-only eligibility needs a separate
correctness decision and control.
