# Molecule-first Visium HD materializer optimization and validation

Date: 2026-07-22

Branch: `performance/molecule-first-materializer-20260722`

Base: `ad8a20c41f60a404d40d5c6a56386ff3c33a6ccb`

Implementation commit: `be0ab3059aaf6d7a1fab20ffd8433f02d2d64d4b`

Test commit: `0e2b62f5caf75ef0bc8ad4f190bd3152608ac8d3`

## Scope and scientific boundary

This change optimizes only `molecule_first_materialize`, after the resolver has
completed candidate preservation, clique formation, posterior assignment,
candidate-specific UMI correction, and molecule collapse. It does not change
candidate enumeration, clique membership, posterior policy, hard/gated calls,
UMI correction, occupancy, molecule identifiers, or resolver outputs.

The full-slide test consumed the existing resolved directory without rerunning
the resolver:

```text
/mnt/pikachu/star-spatial/runs/20260722_human_crc_flex_full_native_policy_comparison_v1/resolved
```

No Space Ranger assignment, corrected barcode/UMI, called-bin field, or vendor
molecule definition is read by the materializer.

## Architecture and integration point

The released resolved-directory CLI remains compatible. Visium HD now uses a
bounded post-collapse path by default; non-HD assays retain the generic path,
but load and write one policy at a time and no longer construct complete Matrix
Market files in an `ostringstream`.

The HD path is:

1. Scan `read_cliques.tsv` once to build the feature axis and compact observed
   coordinate sets. Fine and parent coordinate axes retain the exact lexical
   ordering of the historical `std::set<string>` implementation.
2. Parse each post-collapse policy TSV exactly once and process one policy at a
   time. For feature-sorted resolver output, aggregate only the current
   completed feature chunk before emitting counts. Older resolved directories
   without the streaming marker use the same bounded external transpose without
   assuming feature order.
3. External-sort deterministic 24-byte binary records under
   `--sort-memory-mb` (default 1024 MiB). Records contain matrix key, stable
   source order, and value. Spill directories are mode 0700, private to the
   process, hidden under `--tmp-dir` or the output root, and removed after use.
4. During the fine-scale merge, emit both 8 and 16 um parent records. Thus all
   scales share one policy parse and one fine aggregation; coarse summation
   retains fine-barcode lexical order.
5. Merge once to determine NNZ/mass and again to write Matrix Market entries
   incrementally into an atomic `matrix.mtx.tmp`. Complete matrices are never
   retained as strings or row vectors.
6. Write one feature axis and one barcode axis per scale, then atomically
   hard-link identical axes into all policy directories. The output schema and
   file contents remain unchanged while duplicate physical storage is avoided.

The feature-sorted TSV boundary is the current resolver-to-materializer
integration point: collapsed counts are flushed after each resolver feature
chunk, without reassigning barcodes or recollapsing UMIs. An in-process fused
resolver sink was intentionally not required for this tranche because the
production equivalence run had to reuse the sealed resolved directory.

New optional CLI controls:

```text
--sort-memory-mb INT   combined in-memory spill-record budget; default 1024
--tmp-dir DIR          root for the private execution-only spill directory
```

## Build and binary provenance

The benchmark host was x86-64, little-endian, with a 13th Gen Intel Core
i9-13900KF (24 cores / 32 logical CPUs). The production binary was built with
the repository Makefile (`g++ -O3 -std=c++11`) after a clean tool build.

```text
Binary: flex/tools/molecule_first_resolver/molecule_first_materialize
SHA-256: 8ccaaea8841886da10215559bf68e3b44a1f18847538181a365bc40f839475ce
STAR Suite version string: 1.5.0
Source commit: be0ab3059aaf6d7a1fab20ffd8433f02d2d64d4b
```

## Tests executed

All commands exited zero unless explicitly described as a comparison that
reported the documented summary-only float differences.

```bash
make -C flex/tools/molecule_first_resolver clean
make -C flex/tools/molecule_first_resolver -j8 all
make -C flex/tools/molecule_first_resolver test
tests/run_molecule_first_production_adapter_test.sh
MOLECULE_FIRST_RESOLVER="$PWD/flex/tools/molecule_first_resolver/molecule_first_resolver" \
  tests/run_molecule_first_native_smoke.sh
tests/run_molecule_first_materializer_bounded_test.sh
```

The bounded test runs the feature-streaming materializer twice, checks exact
directory determinism, exercises HD row/column lexical ordering, verifies all
12 matrices, checks scale mass conservation, verifies shared-axis hard links,
and checks private-spill cleanup.

ASan/UBSan also passed:

```bash
g++ -Iflex/source -Icore/legacy/source/htslib \
  -DSTAR_SUITE_VERSION=\"1.5.0\" -O1 -g -std=c++11 -Wall -Wextra \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -o /tmp/molecule_first_materialize_asan_20260722 \
  flex/tools/molecule_first_resolver/molecule_first_materialize.cpp
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  MOLECULE_FIRST_RESOLVER="$PWD/flex/tools/molecule_first_resolver/molecule_first_resolver" \
  MOLECULE_FIRST_MATERIALIZER=/tmp/molecule_first_materialize_asan_20260722 \
  tests/run_molecule_first_materializer_bounded_test.sh
```

The existing Visium HD Flex 100K resolved fixtures tested both compatibility
routes with `--sort-memory-mb 1`, forcing multiple spill chunks. Repeated new
runs were byte-identical. Against the original `ad8a20c` materializer on the
same feature-streaming resolved directory, all 36 matrix/axis files were
byte-identical. The original took 2.06 s / 108,620 KiB RSS; the bounded path
took 0.53 s / 18,048 KiB RSS. A small generic exact-mode run was also
byte-identical to the original materializer. No full-slide exact-mode run was
performed.

## Full-slide command and outputs

```bash
/usr/bin/time -v \
  flex/tools/molecule_first_resolver/molecule_first_materialize \
  --resolved-dir /mnt/pikachu/star-spatial/runs/20260722_human_crc_flex_full_native_policy_comparison_v1/resolved \
  --out-dir /mnt/pikachu/molecule_first_materializer_fullslide_1mm_cr_20260722_v1 \
  --assay visium-hd \
  --umi-mode 1mm_cr \
  --sort-memory-mb 1024
```

Timing record:

```text
/mnt/pikachu/molecule_first_materializer_fullslide_1mm_cr_20260722_v1.time.txt
```

Reference:

```text
/mnt/pikachu/star-spatial/runs/20260722_human_crc_flex_full_native_policy_comparison_v1/policy_mex_1mm_cr
```

## Full-slide equivalence

Every one of the 12 `matrix.mtx` files, 12 `features.tsv` paths, and 12
`barcodes.tsv` paths is byte-identical to the reference (`36/36`). This proves
entry-wise equality, axis equality, matrix dimensions/order, and barcode
hierarchy equality. In particular, soft expected matrix entries have zero
entry-wise difference; only summary accumulation order differs as described
below.

| Policy | Scale | Dimensions (features x barcodes) | NNZ | Reference mass | Optimized summary mass | Matrix/axes |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| strict | 2 um | 18,072 x 8,658,043 | 64,586,042 | 65,842,882 | 65,842,882 | byte-identical |
| strict | 8 um | 18,072 x 678,773 | 59,726,279 | 65,842,882 | 65,842,882 | byte-identical |
| strict | 16 um | 18,072 x 175,292 | 54,087,967 | 65,842,882 | 65,842,882 | byte-identical |
| soft_expected | 2 um | 18,072 x 8,658,043 | 92,065,457 | 78,197,638.581787929 | 78,197,638.581749737 | byte-identical |
| soft_expected | 8 um | 18,072 x 678,773 | 75,667,696 | 78,197,638.581868380 | 78,197,638.581868812 | byte-identical |
| soft_expected | 16 um | 18,072 x 175,292 | 65,712,241 | 78,197,638.581862047 | 78,197,638.581848070 | byte-identical |
| hard | 2 um | 18,072 x 8,658,043 | 76,487,365 | 78,142,032 | 78,142,032 | byte-identical |
| hard | 8 um | 18,072 x 678,773 | 70,237,214 | 78,142,032 | 78,142,032 | byte-identical |
| hard | 16 um | 18,072 x 175,292 | 63,065,350 | 78,142,032 | 78,142,032 | byte-identical |
| gated_hard | 2 um | 18,072 x 8,658,043 | 69,647,101 | 71,066,741 | 71,066,741 | byte-identical |
| gated_hard | 8 um | 18,072 x 678,773 | 64,211,057 | 71,066,741 | 71,066,741 | byte-identical |
| gated_hard | 16 um | 18,072 x 175,292 | 57,934,607 | 71,066,741 | 71,066,741 | byte-identical |

Strict, hard, and gated-hard mass is exact at every scale. The soft summary
deltas (optimized minus reference) are:

```text
2 um:  -3.819167613983154e-05
8 um:   4.321336746215820e-07
16 um: -1.397728919982910e-05
```

The largest absolute summary delta is `3.8192e-05`. The declared mass tolerance
is `1e-9 * mass = 0.0781976386`, so the delta is about 2,047 times smaller than
the tolerance. The optimized soft mass range across scales is `0.0001190752`,
also far below the same tolerance. This is only the order of the final summary
reduction (matrix-order in the bounded merge versus feature-order in the old
map); the matrix values themselves are byte-identical and therefore conserve
the same posterior mass entry by entry.

## Performance and disk

| Metric | Original | Bounded | Change |
| --- | ---: | ---: | ---: |
| Wall time | 58:46.95 | 11:09.44 | 5.27x faster (81.0% lower) |
| User CPU | 2,394.99 s | 296.88 s | 8.07x lower |
| System CPU | 68.80 s | 53.90 s | 1.28x lower |
| Total CPU | 2,463.79 s | 350.78 s | 7.02x lower (85.8% lower) |
| Maximum RSS | 63,381,476 KiB (60.45 GiB) | 1,221,652 KiB (1.17 GiB) | 51.9x lower (98.1% lower) |
| Major page faults | 596,978 | 0 | eliminated in this run |
| Final allocated output | 13,319,405,568 B | 12,764,155,904 B | 555,249,664 B lower (4.17%) |
| Final apparent output | 13,319,287,227 B | 12,764,071,827 B | 555,215,400 B lower |

The optimized feature file has link count 12; each scale-specific barcode file
has link count 4. Matrix files remain independent atomic outputs.

## Failures and caveats

- The first 100K development comparison exposed row-component lexical ordering
  around `1`, `10`, and `1000`. The comparator was corrected to include the
  underscore delimiter used by the canonical candidate string, followed by a
  clean rebuild. The tracked test now seals this case. All subsequent 100K and
  full-slide axes are byte-identical.
- External sorting intentionally trades bounded memory for temporary I/O. The
  optimized `/usr/bin/time` record reports about 32.19 GiB of filesystem output
  versus 12.40 GiB for the original; private spill files were removed after
  success. Final output storage is smaller because axes are shared.
- The three soft summary strings are not byte-identical because of reduction
  order. All soft matrix entries are byte-identical, and the summary deltas are
  well inside the frozen floating-point tolerance.
- The current direct integration boundary is the feature-sorted post-collapse
  TSV stream, not an in-process resolver callback. Each large TSV is parsed once
  and no resolver stage is repeated. A future fused runner could avoid writing
  those TSVs, but is not needed to obtain the measured memory/performance result
  or to consume already-resolved directories.
- Only the required full-slide candidate-specific 1MM-CR mode was benchmarked.
  Exact mode remains supported and passed a small backward-compatibility check.
