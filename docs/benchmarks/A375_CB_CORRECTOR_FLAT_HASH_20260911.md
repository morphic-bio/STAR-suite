# A375 barcode-corrector storage benchmark, 2026-09-11

Status: implemented and validated on full A375 in the development worktree. No release.

Full runtime falls by 44–47 seconds (23–24%), peak RSS falls by about 4.4 GiB, and all count matrices
and guide-call tables remain exactly equal. The measured success-message-to-exit interval falls from
about 28.4 seconds to roughly one second. Most of the other savings come from table construction.

The owner authorized two changes together: replace the barcode-correction tables with khash and
flatten the ambiguity candidates. The control is the corrected A375 build from
`HANDOFF_A375_HIERARCHY_CORRECTNESS_20260911.md`, including the zero-N and deterministic-learning fixes.
Only the storage representation changes in this comparison. Scheduler settings match within each pair.

## Implementation

`flex/source/solo/CbCorrector.{h,cpp}` now uses khash for exact matches, H1 matches and ambiguity ranges.
The keys remain packed 32-bit barcodes, hashed with the existing Wang integer mixer. Exact values
remain zero-based whitelist indices; H1 values remain one-based indices with zero denoting ambiguity.
The tables have RAII ownership and check allocation failures.

Ambiguity candidates occupy one contiguous `vector<uint32_t>`, with an offset/count per ambiguous key.
The first construction pass counts candidates; the second repeats the original whitelist/position/base
traversal to fill each range in its original order. No per-key vectors are allocated during construction
or retained afterward. Duplicate multiplicity and candidates above the correction limit are preserved.
Read-only views replace the former map accessor. The opt-in legacy ambiguity export in
`core/legacy/source/ParametersSolo.cpp` consumes those views and retains its existing output convention.

Exact-first priority, the five-candidate ambiguity limit, N expansion, FASTQ/CBQ packing order and
Bayesian correction policy are unchanged. The whitelist string storage also remains unchanged.

## Build and correctness

Clean Chromap-enabled core build:

```sh
make -C core/legacy/source clean
make -C core/legacy/source -j8 STAR
```

Frozen binaries (full SHA256):

- Control: `ce9d8142f6e8f0f52d4755b669265166f4a08e673c65650be362e8ffe2577a3c`
- New storage: `9fcde3060a5d5102831c1450a0d0551bf504933a31d04cae14d443c0c2e156b0`

`tests/run_cb_corrector_storage_test.py` compiles the current implementation alongside a renamed,
preserved copy of the old one. Both ordinary and ASan/UBSan runs pass. Coverage includes exhaustive
six-base A/C/G/T/N queries, 16-base mutations, exact/H1/N correction, both packing conventions,
duplicate input barcodes, ordered candidate lists, the five/six-candidate boundary, empty inputs and
concurrent read-only queries. All result fields and candidate ordering must match.

Full A375 uses 32 permits, native BGZF, four inflater workers per parent, MAP/FEATURE floors 24/8,
GeneFull, CR-compatible multimapper rescue and no BAM or sidecars. Input, reference, environment and
configuration identities are retained in each manifest. Both fixed-four and balanced arms pass exact
equality for all six raw/filtered MEX matrices and all three guide-call CSVs against their respective
controls. The two new arms also agree with each other. Feature totals remain 9,275,612 assigned counts,
3,221,846 deduplicated counts, 109,746 raw barcodes and 348,305 unmatched reads. GEX retains 1,170 cells.

## Component measurement

The real A375 whitelist contains 3,686,400 entries, yielding 88,270,080 H1 keys and 42,938,880 ambiguous
keys. The new candidate array contains 131,616,000 indices. Both component probes create and join a
thread before destruction so glibc's single-thread allocator shortcut is disabled.

| Measurement | Before | After |
| --- | ---: | ---: |
| Construct correction tables | 43.91 s | 26.44 s |
| Destroy ambiguity table and candidate storage | 4.411 s | 0.058 s |
| Destroy H1 table | 7.473 s | 0.039 s |
| Destroy exact table | 0.056 s | 0.004 s |
| Destroy both whitelist string copies | 0.061 s | 0.062 s |
| Total measured destruction | 12.000 s | 0.163 s |
| Probe peak RSS | 7.61 GiB | 3.31 GiB |

These are isolated component measurements. The full-process exit interval is measured separately.

## Full A375 measurements

| Full-run measurement | Fixed-four before | Fixed-four after | Balanced before | Balanced after |
| --- | ---: | ---: | ---: | ---: |
| Total wrapper wall | 192.10 s | 145.47 s | 190.60 s | 146.21 s |
| Peak RSS | 44.20 GiB | 39.82 GiB | 44.20 GiB | 39.80 GiB |
| Startup to mapping | 54 s | 36 s | 55 s | 36 s |
| Of which explicit genome load to mapping | 8 s | 8 s | 8 s | 8 s |
| Mapping | 91 s | 89 s | 89 s | 90 s |
| Feature API, including gather/cleanup | 24.48 s | 26.74 s | 24.69 s | 26.58 s |
| Solo counting | 5 s | 5 s | 5 s | 5 s |
| Joint merge/finalization and guide calling | 13 s | 13 s | 13 s | 12 s |
| Feature completion to Solo completion | 66 s | 63 s | 67 s | 67 s |
| Success message to process exit | ~28.39 s | ~1.10 s | ~28.45 s | ~0.61 s |
| Whole-run average CPU cores | 14.33 | 18.49 | 14.43 | 18.50 |

These are one completed execution per configuration, compared with preserved controls. Log-derived
phase boundaries have one-second resolution; the wrapper and feature API durations have finer timers.
The exit intervals include timestamp rounding and wrapper polling, so subsecond precision should not
be interpreted as an exact destructor measurement.

The main gains are startup and exit. Mapping remains near its previous duration. Feature API time
increases by about two seconds while overlapping GEX, without extending overall runtime. The gathers
remain included in the phase metrics. These runs do not establish a scheduler improvement or a faster
read-lookup loop; they establish the benefit of the storage changes together. The original isolated
12-second teardown did not predict the full-process cost: almost the entire observed 28-second exit
gap disappears after this replacement. No cleanup was skipped.

## Provenance

Artifact root: `/home/lhhung/pf_larry_regression_20260911/hierarchy/`.

- `cbcorrector_flat_hash/baseline/`: preserved source/object, pre-change patch and baseline identities.
- `cbcorrector_flat_hash/build/`: frozen STAR/object, source snapshot, clean-build commands and logs.
- `cbcorrector_flat_hash/{unit,sanitizer,probe}/`: source identities, commands, execution records and logs.
- `exit_cleanup_probe/threaded_allocator/`: preserved component control.
- `a375_corrected_{fixed4,balanced4}/`: preserved full-run controls.
- `a375_flat_hash_{fixed4,balanced4}/`: new full-run outputs, argv, environment, CPU samples and GNU time.
- `cbcorrector_flat_hash/analyze.py`: phase, memory, CPU and strict output comparisons.
- `cbcorrector_flat_hash/results.json`: combined before/after component and full-run measurements.
- `cbcorrector_flat_hash/*comparison.json`: exact matrix and guide-table equality results.

All timed runs are serialized and use fresh output directories. Each configuration is executed once.
The comparator requires an explicit `--allow-different-binaries` for the storage comparison, records
both binary identities, and retains strict biological equality. This is development validation, not a
new paper release or a repeated estimate of run-to-run variance.
