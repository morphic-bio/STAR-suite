# Runbook: post-collapse molecule-first production integration

Date: 2026-07-17  
Branch: `agent/post-collapse-production-integration`  
Base: STAR Suite `master` at `b678a61` (released resolver: `v1.5.0`)  

## Objective

Move the opt-in `molecule_first_resolver` from a validated normalized-ledger
release surface into reproducible production adapters for Chromium 3-prime,
Chromium Flex, classic Visium, and Visium HD. The operation order remains:

```text
raw candidate-bearing reads
  -> globally candidate-consistent read cliques
  -> candidate-specific UMI correction and occupancy
  -> strict / soft expected / hard / gated-hard matrices
  -> spatial-bin inclusion or cell calling
```

The frozen STAR/STARsolo/Flex decoder remains the default compatibility path.
The production adapters are enabled only by an explicit command or manifest.

The mandatory invariant remains:

> Ambiguous molecules may refine a conservative broad field, but they may not
> manufacture unsupported fine-scale structure.

## Frozen evidence and prohibited priors

The production adapters may use only:

- the finite barcode candidates supported by the observed barcode sequence;
- observed bases and Phred qualities;
- PCR-inclusive exact-read frequencies as a read prior;
- STAR/STAR-Flex feature identity and raw UMI for clique construction;
- candidate-specific exact or directional one-mismatch UMI correction.

They may not use expression, called-cell membership, cell type, spatial
neighborhood, image, cluster, graph, GPU, Space Ranger assignment, or Cell
Ranger assignment as a prior. Vendor outputs remain final compatibility
diagnostics only.

## Tranche 1: portable release containers

The failed Docker publication run `29564366223` authenticated successfully but
failed while compiling the ARM64 target because
`core/legacy/source/Makefile` supplied the x86-only `-mavx2` flag to the SIMDe
OPAL source. The source already includes `simde_avx2.h`.

Implementation:

1. Detect the compilation target architecture.
2. Use `-mavx2` for x86-64 and no x86 ISA flag for ARM64, allowing SIMDe to
   select its ARM implementation.
3. Add a compile-time assertion/test recording which SIMD backend was built.
4. Build and smoke-test `suite-base` for `linux/amd64` and `linux/arm64`.
5. Publish Docker tags only after both targets pass.

No release workflow is permitted to silently reduce the declared platform
set.

## Tranche 2: spatial production adapters

### Classic Visium

Classic Visium uses the Chromium-family adapter with the slide spot whitelist.
The adapter reads STAR-produced BAM records, preserves exact or one-substitution
whitelist candidates, computes Phred likelihoods, counts exact reads in a
first pass, and writes the v1 normalized ledger. It never reads corrected CB
or UB tags.

### Visium HD

The Visium HD adapter consumes the existing raw R1 composite candidates joined
to STAR/STAR-Flex feature assignments. It converts each coordinate candidate
to `s_002um_<row>_<col>` and applies the frozen factorized exact-read prior.
The adapter writes the same normalized ledger consumed by the native resolver.

Both spatial paths materialize all four policy matrices. Visium HD matrices
must conserve mass exactly (within declared floating tolerance) when
aggregated from 2 um to 8 um and 16 um. Strict remains the conservative
fine-scale field.

## Tranche 3: single-cell production adapters

Chromium 3-prime and Flex use the same BAM-to-ledger implementation. Flex is
run on each already-demultiplexed sample surface so a clique cannot cross
sample tags. The adapter requires STAR provenance in the BAM header and uses
raw barcode (`CR`), raw barcode quality (`CY`), raw UMI (`UR`), and feature
(`GX` by default) only.

Exact whitelist observations remain single-candidate controls. Non-exact
observations retain every one-substitution whitelist candidate. The exact-read
count for a candidate is counted before deduplication and therefore retains
PCR amplification, but is applied by the resolver once per read clique.

The hard, gated-hard, and strict integer MEX products may be passed to
`STAR --runMode soloCellFiltering ... --soloCellFilter EmptyDrops_CR` after
materialization. Soft expected counts remain a real-valued audit/analysis
matrix; the current integer EmptyDrops loader must not silently round it and
call it equivalent.

## Tranche 4: regression matrix

Required tracked synthetic tests:

- architecture-aware SIMD selection for x86-64 and ARM64;
- BAM provenance and raw-tag enforcement;
- exact controls remain single-candidate;
- all supported one-substitution candidates are retained;
- corrected `CB`/`UB` and called-cell fields cannot enter the ledger;
- PCR-inclusive exact counts agree across every candidate row;
- normalized ledger order is deterministic;
- strict/soft/hard/gated products contain no unsupported barcode;
- policy MEX mass and dimensional checks;
- Flex sample boundaries cannot be crossed;
- Visium HD 2/8/16 um mass conservation;
- default-off STARsolo and Flex outputs remain unchanged.

Required fixture smokes:

| Assay | Fixture | Primary gate |
| --- | --- | --- |
| Chromium 3-prime | frozen 100K challenge/source fixture | native/reference equality plus pre-EmptyDrops materialization |
| Chromium Flex | frozen 100K challenge/source fixture | native/reference equality per sample |
| Classic Visium | current 100K Visium fixture | supported-candidate and matrix determinism |
| Visium HD 3-prime | frozen 100K H0 candidate ledger | sealed native/reference equality and 2/8/16 um mass |
| Visium HD Flex | frozen 100K candidate ledger | sealed native/reference equality and 2/8/16 um mass |

Benchmark runs are serialized and use fresh output directories. Missing
external fixtures are explicit skips, never passes.

## Tranche 5: frozen paper evidence

The implementation commit, executable version, input paths, SHA-256 digests,
resolved parameters, command lines, output summaries, and compatibility
diagnostics are copied as small immutable records into
`visium-hd-processing-paper`. Raw FASTQs, BAMs, matrices, and vendor outputs
remain external.

The frozen paper bundle must include:

- the fitted sequence-error generator and synthetic recovery summaries;
- the independent human and mouse biological oracle summaries;
- ambiguous-fraction contrasts in addition to whole-matrix metrics;
- strict-field subset checks and soft/hard/gated mass accounting;
- standard/Flex/Visium/Visium HD regression summaries;
- Space Ranger and Cell Ranger only as compatibility diagnostics;
- exact STAR Suite, processing, and paper repository commits.

Checksums are regenerated only after every included result passes its schema
and invariant checks.

## Planned implementation surfaces

- `core/legacy/source/Makefile`: architecture-aware SIMD default.
- `docker/Dockerfile`: target architecture propagation and build audit.
- `flex/tools/molecule_first_resolver/`: raw-tag candidate-ledger adapter and
  shared policy materializer/runner.
- `tests/molecule_first/`: synthetic BAM/ledger/MEX and four-assay integration
  harnesses.
- `visium-hd-processing`: Visium HD normalization and manifest-driven native
  execution adapter.
- `visium-hd-processing-paper`: sealed provenance and result records only.

## Execution ledger

| Gate | Command/artifact | Result |
| --- | --- | --- |
| Plan frozen before code | this runbook | complete |
| Existing native/reference unit conformance | `tests/run_molecule_first_reference_conformance.sh` | passed in v1.5.0 |
| Existing four 100K normalized-ledger conformance | `20260717_star_suite_molecule_first_integration_v1` | passed in v1.5.0 |
| Docker Hub authentication | release run `29564366223`, `Login to Docker Hub` | passed |
| ARM64 container build | release run `29564366223` | failed before this tranche: x86-only `-mavx2` |
| Architecture flag test | `tests/run_simd_arch_flags_test.sh` | passed: x86-64 selects `-mavx2`; AArch64 selects the portable SIMDe path |
| Clean ARM64 container build | `docker buildx build --platform linux/arm64 --target suite-base .` | passed after repository-wide removal of stale host objects and archives |
| Clean AMD64 container build | `docker buildx build --platform linux/amd64 --target suite-base .` | passed |
| Production adapter synthetic test | `tests/run_molecule_first_production_adapter_test.sh` | passed |
| Existing reference conformance after adapter changes | `tests/run_molecule_first_reference_conformance.sh` | passed |
| Installed native resolver smoke | `tests/run_molecule_first_native_smoke.sh` | passed with all three installed executables |
| Default-off standard STARsolo fixture | existing standard smoke harness | passed |
| Default-off Flex fixture | frozen Flex 100K smoke, resolver both disabled and enabled explicitly | passed; both runs reached `ALL DONE` |
| Four-assay fixture regression | pending | pending |
| Frozen paper bundle validation | pending | pending |

The container rebuilds exposed two additional portability hazards after the
compiler flag was corrected: a stale host-built `opal.o`, followed by a stale
host-built `libprocess_features.a`. The Docker context now excludes object and
archive files, and the image build cleans all nested build surfaces before
compiling. This makes the platform claim a property of the actual target
build, rather than of whichever objects happened to be present in the source
tree.

## Production commands

The tracked adapters are deliberately explicit rather than automatic:

```bash
molecule_first_bam_ledger \
  --input-bam Aligned.sortedByCoord.out.bam \
  --whitelist 3M-february-2018.txt \
  --output-ledger read_candidates.tsv \
  --output-summary ledger_summary.tsv

molecule_first_resolver \
  --ledger read_candidates.tsv \
  --output-dir resolved

molecule_first_materialize \
  --resolver-dir resolved \
  --output-dir matrices

scripts/run_molecule_first_cell_calling.sh \
  --star STAR \
  --policy-root matrices \
  --out-root calls \
  --policy hard
```

Flex ledger generation additionally requires `--sample-id` and must be run on
one already-demultiplexed sample at a time. Visium HD normalization and the
manifest/provenance wrapper live in `visium-hd-processing`; they produce the
same normalized ledger and invoke these native executables without changing
the resolver.

## Recovery after interruption

1. Read this runbook and
   `docs/MOLECULE_FIRST_BARCODE_ASSIGNMENT_IMPLEMENTATION_20260717.md`.
2. Confirm the worktree is on `agent/post-collapse-production-integration`
   and based on `b678a61` or a documented descendant.
3. Run a clean rebuild before diagnosing any native failure.
4. Never reuse a partially written adapter or benchmark output directory.
5. Preserve the existing decoder default and the frozen resolver parameters.
6. Do not copy large or vendor-controlled inputs into a Git repository.
