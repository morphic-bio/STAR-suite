# OCM count routing without intermediate MEX reads

Date: 2026-09-10. Baseline OCM implementation: `8ec9945`.
Follow-up to [the cross-module performance runbook](../RUNBOOK_CROSS_MODULE_PERFORMANCE_20260910.md).

## Implementation

Native per-sample EmptyDrops now reads the pooled raw matrix once. That pass
writes the required pooled raw MEX and routes compact `(gene, cell, count)`
records to each configured sample, including tag unions. The caller prepares
its sparse arrays directly from those records. Required sample raw and filtered
MEX files are written after calling; no generated sample MEX is read back.

The records preserve source order. Sparse preparation preserves the original
within-column order, detects nonzero columns before allocating the count arrays,
and keeps the existing barcode sort, feature universe, OrdMag/ambient/MC settings,
seeds and decisions. Export rounding and zero-entry behavior are preserved for
each existing path. Invalid dimensions, coordinates, count/UMI overflow,
truncation and gzip errors now fail explicitly.

`--ocmCellCallMaxMemory` also bounds retained count blocks: at most one quarter
of its value is available for those blocks. Stores that outgrow this allowance
spill in 12-byte binary records using bulk buffered I/O. Spill files are created
in the existing OCM temporary directory and immediately unlinked; handles close
on success or exception. This is process-local scratch, not a new persistent
cache format. Retained blocks and fixed 48 KiB per-store I/O buffers are reserved
before admitting caller matrices. CPU permits and the oversized-sample-runs-alone
rule are preserved. Output directories are created by the coordinator before
parallel tasks start, so workers never write directory messages to STAR's shared
log stream; sample messages are buffered and emitted in configuration order. As before, this is an estimated matrix admission budget,
not a hard process RSS limit; axes, routing maps, compression and other caller
allocations also consume memory.

When pooled filtered calls already exist, the materializer routes the pooled
raw output plus all sample raw outputs in one source scan, and all sample
filtered outputs in one filtered-source scan. It no longer scans each source
matrix twice per sample. This path still uses temporary text bodies to assemble
MEX headers; they are checked for write failures and removed on exceptions.
Velocyto routing and downstream directory/mirror contracts are preserved.

## Larger regression and measurements

All measurements were serialized on the same i9-13900KF host, affinity CPUs 0–2,
three workers, no alignment, reference index, BAM or sidecar. These are single
measurements, not replicated estimates. Controls already completed during the
cross-module work were reused; the new existing-call input had one control run.

Input: real A375 counts, 38,606 features, 290,020 detected barcodes and 7,533,915
nonzeros. Synthetic OCM suffixes divide the counts into 116,008 / 87,006 / 87,006
/ 0 columns plus a 203,014-column union. This tests routing and execution on real
counts; it is not a biological OCM benchmark. Native calling uses the same 1,000
simulations as its saved control. The existing-call arm uses the saved global
caller's 1,188 barcodes, preserving 6,274,080 filtered nonzeros.

| Materialization | Control wall / peak RSS | Final wall / peak RSS | Result |
| --- | --- | --- | --- |
| Native calling, 100,000,000-byte matrix budget | 28.7536 s / 263,416 KiB | **8.5893 s / 278,992 KiB** | **70.1% less wall time**, 95 files exact; four nonempty samples spill |
| Existing filtered calls, same 100,000,000-byte budget | 27.5632 s / 158,668 KiB | **9.4468 s / 184,340 KiB** | **65.7% less wall time**, 75 files exact |
| Native calling, default 1 GiB budget | No matching-budget timing control | **7.9689 s / 467,812 KiB** | 95 files exact against saved calls; all records retained, no spill |

The native sample callsets remain **464 / 342 / 352 / 0**, with **826** for the
union. Exact comparison includes counts, axes, caller diagnostics, routing JSON,
CSV and downstream mirrors. Gzip payloads are compared after decompression;
only the copied fixture's config path is normalized in the summary JSON.

This is a runtime optimization, not a memory-reduction claim. At the matched
100 MB admission setting, native peak RSS rises 15.2 MiB and existing-call peak
RSS rises 25.1 MiB. The higher-budget native arm illustrates the memory/time
tradeoff and is not a like-for-like speed comparison with the 100 MB control.
The existing-call arm's temporary bodies also increase recorded filesystem
writes (418,560 to 1,216,344 units of 512 bytes); avoiding repeated parsing still
reduces its CPU time from 27.03 to 9.07 seconds. Native spill-mode filesystem
writes decrease from 1,003,232 to 718,440 such units.

## Correctness and failure gates

- Tiny native retained/spilled fixtures: 155 files exact each, including caller
  diagnostics, empty groups, union samples and downstream Velocyto layers.
- Existing-call tiny fixture: 135 files exact.
- Interleaved coordinates, duplicate entries, zero-UMI known tags, fractional
  export/caller rounding, plain/gzip input and one/multiple worker budgets:
  exact outputs against frozen controls.
- Invalid rows/columns/axes, truncated records, malformed entries, count and
  per-barcode UMI overflow, corrupt CRC and truncated gzip: rejected before
  publishing sample matrices. Existing-call axis failures also clean up scratch.
- The independent count-store test passes AddressSanitizer and UBSan: empty,
  exact/partial block boundaries, shared retention limits, immediate/late spill,
  replay following exceptions, invalid spill directory and allocation accounting.
- Full STAR and the native OCM materializer harness build successfully with
  Chromap enabled. The initial sanitizer link exposed a C++11 constant-definition
  issue; enum constants fixed it before final validation. No output-parity failures
  occurred during implementation.

Reusable regression sources are `tests/run_ocm_direct_counts_regression.py` and
`tests/test_ocm_count_store.cpp`. The Python runner takes `--candidate`,
`--baseline`, `--output`, and optional `--saved-native`, `--saved-existing`,
`--saved-ordered` controls. Output directories must be fresh. On this host,
reuse saved controls and follow the identical-execution rule in AGENTS.local.md.
Do not rerun the whole prior cross-module suite to test this isolated change.

## Evidence

Root: `/mnt/pikachu/star_suite_paper/analysis/ocm_direct_counts_20260910`.

- `bin/`: frozen final STAR and OCM harness.
- `regression_accepted/PARITY.json`: exact-output and malformed-input results.
- `STORE_COMPLETE.json`, `store_sanitizer.log`: sanitizer test completion.
- `native_accepted`, `native_default_accepted`, `existing_accepted`: final commands,
  binary hashes, wall/CPU/RSS/I/O data, logs and all outputs.
- `existing_control`: the new existing-call comparator execution.
- `ACCEPTED_LARGE_RESULTS.json`: final measurements and exact file counts.
- `large_inventory.json`, `existing_inventory.json`: input identities/provenance.
- Initial build/gate/measurement attempts are retained. The final large arms
  follow added axis, temporary-file and shared-log checks; frozen controls were not repeated.

The saved native timing/output control remains in
`../cross_module_performance_20260910/ocm_large/parallel`.
No Cell Ranger source was inspected, and no caller model or concordance target
was changed.
