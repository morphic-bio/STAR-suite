# Runbook: Chromap ATAC CBQ In-Memory Integration

Date: 2026-05-28

Status: planned. This runbook defines the ATAC-only path for replacing the
current CBQ-to-Chromap FASTQ materialization adapter with production CBQ
support.

## Goal

Support Chromap ATAC from CBQ without temporary FASTQ files:

```text
shared CBQ reader -> decoded read-batch view -> libchromap ATAC provider ->
Chromap ATAC mapping/fragments/peaks
```

Production CBQ support means the consumer receives either an in-memory decoded
sequence view or a direct adapter into native consumer structs. FASTQ
materialization is allowed for tests and compatibility shims, but it is not the
production surface.

## Scope

In scope:

- ATAC paired-end reads with a barcode stream.
- STAR-suite multiome `chromapAtac*` integration.
- Chromap ATAC fragments, BAM/CRAM, binary fragment sidecar, MACS3 FRAG peak
  materialization, and ATAC peak-MEX/evidence flows that already exist.
- Native CBQ reader reuse; no separate Chromap CBQ parser.

Out of scope for this runbook:

- Hi-C / pairs.
- ChIP or bulk non-barcode Chromap modes.
- Generic Chromap single-end input.
- New BINSEQ variants beyond CBQ.
- Forcing existing FASTQ paths through the CBQ contract.

## Current State

STAR-suite currently has an ATAC-specific libchromap contract:

- `core/features/libchromap_contract/include/star_chromap_contract.h`
- `core/features/libchromap_contract/src/star_chromap_contract.cpp`
- `core/legacy/source/star_chromap_orchestration.cpp`

That contract is path based today: it carries `read1_fastqs`, `read2_fastqs`,
and `barcode_fastqs`, then converts them into Chromap
`MappingParameters.read_file1_paths`, `read_file2_paths`, and
`barcode_file_paths`.

Chromap-suite has `libchromap` entrypoints:

- `/mnt/pikachu/Chromap-suite/src/libchromap.h`
- `/mnt/pikachu/Chromap-suite/src/libchromap.cc`
- `/mnt/pikachu/Chromap-suite/src/chromap_lib_runner.cc`

However, verify the active Chromap-suite branch before starting work. In the
current local checkout inspected for this runbook, the main `chromap` CLI still
constructs `Chromap` and dispatches `MapSingleEndReads` /
`MapPairedEndReads` directly in `src/chromap_driver.cc`, while
`chromap_lib_runner` calls `chromap::RunMapping()`. The first implementation
phase is therefore to make `libchromap` the ATAC mapping source of truth for
the CLI, runner, and STAR-suite contract.

Chromap internally maps from `SequenceBatch` objects. That makes an in-memory
CBQ provider tractable: the adapter can populate or feed the same sequence,
quality, name, barcode, and read-id surfaces that the FASTQ loader currently
fills.

## Design Policy

- One shared STAR-suite native CBQ reader remains the only CBQ parser.
- Chromap receives decoded sequence state, not FASTQ text.
- FASTQ remains on its current optimized Chromap and STAR-suite paths.
- The new provider must be ATAC-specific until non-ATAC modality contracts are
  intentionally designed.
- CLI compatibility is preserved: existing FASTQ Chromap commands must keep
  working and must remain parity baselines.

## Target Architecture

Add an ATAC read-provider boundary in Chromap-suite:

```text
ChromapAtacReadProvider
  next_batch(read_batch1, read_batch2, barcode_batch) -> count | EOF | error
```

The concrete providers are:

- `FastqAtacReadProvider`: wraps the existing FASTQ-backed `SequenceBatch`
  loader and preserves current behavior.
- `ExternalAtacReadProvider` or `CbqAtacReadProvider`: receives decoded
  in-memory records from STAR-suite and fills the same Chromap sequence batch
  state without path-based FASTQ loading.

The provider must preserve:

- synchronized R1/R2/barcode records;
- read names and stable read IDs;
- sequence and quality strings/spans;
- barcode effective range and read-format behavior (`bc:start:end[:strand]`);
- barcode whitelist/correction behavior;
- read order as emitted by the selected provider;
- low-memory ATAC spill behavior;
- STAR dynamic permit hooks;
- fragment, BAM/CRAM, binary sidecar, and peak materialization behavior.

## Phase 1: Make libchromap The ATAC Dispatch Owner

Goal: ensure there is one ATAC mapping implementation path before adding CBQ.

Tasks:

- Move ATAC mapping dispatch from `chromap_driver.cc` into `libchromap` if it is
  still duplicated.
- Keep CLI argument parsing in `chromap_driver.cc`, but convert parsed values
  into `MappingParameters` and call `chromap::RunAtacMapping()` for ATAC
  presets.
- Keep `chromap_lib_runner` on the same `RunAtacMapping()` / `RunMapping()`
  path.
- Centralize ATAC validation in libchromap. CLI, runner, and STAR-suite should
  report the same missing-input and invalid-combination errors.
- Do not change non-ATAC dispatch except where needed to keep existing CLI
  behavior compiling and tested.

Exit criteria:

- Chromap CLI ATAC FASTQ output is byte- or tuple-equivalent to the pre-change
  CLI on the PBMC 3K 100K fixture.
- `chromap_lib_runner` ATAC output matches CLI output for the same fixture.
- STAR-suite `WITH_CHROMAP=1` ATAC fixture still passes.
- Negative validation cases return the same structured failures through CLI and
  lib runner.

## Phase 2: Add The ATAC Read-Provider Interface

Goal: separate ATAC mapping from FASTQ file loading without changing output.

Tasks:

- Introduce a small ATAC provider interface in Chromap-suite. Keep it narrow:
  paired-end ATAC plus optional barcode stream.
- Refactor `Chromap::LoadPairedEndReadsWithBarcodes()` or the surrounding loop
  so ATAC mapping can ask a provider for the next synchronized batch.
- Implement `FastqAtacReadProvider` first using existing `SequenceBatch`
  loading. This should be behavior-preserving.
- Keep provider-owned batch size and read-id accounting identical to current
  `SequenceBatch` semantics unless a change is explicitly justified.
- Preserve barcode effective range application. If this stays inside
  `SequenceBatch`, add a supported way to fill `SequenceBatch` records from
  memory and then apply the same range logic.

Exit criteria:

- ATAC FASTQ CLI parity is unchanged after the provider refactor.
- Existing low-memory and dual BAM/fragments paths still pass.
- Existing barcode correction and read-format tests still pass.

## Phase 3: Add STAR-suite CBQ To Chromap ATAC Adapter

Goal: replace CBQ-to-FASTQ materialization with a production in-memory path.

Tasks in STAR-suite:

- Extend the ATAC contract with an input-source choice:
  - FASTQ paths, current behavior;
  - CBQ decoded provider, new behavior.
- Reuse the shared CBQ reader and `CbqReadBatchView`.
- Define the ATAC CBQ source shape. For current multiome ATAC this is expected
  to be either:
  - one paired-read CBQ for genomic R1/R2 plus one barcode CBQ; or
  - a future three-stream/three-mate CBQ shape if real artifacts require it.
- Add an adapter that maps CBQ records into the Chromap ATAC provider without
  copying through FASTQ text.
- Preserve `--chromapAtacReadFormat`, barcode whitelist, barcode translation,
  output paths, sidecar, peak-MEX, and dynamic permit settings.

Tasks in Chromap-suite:

- Expose a libchromap ATAC entrypoint that accepts the provider or a provider
  factory.
- Ensure provider lifetime covers the full mapping run and that batches do not
  dangle while Chromap worker threads consume them.
- Return structured errors for mismatched R1/R2/barcode record counts and
  malformed records.

Exit criteria:

- The existing synthetic CBQ Chromap smoke no longer materializes FASTQ for the
  production path.
- FASTQ-vs-CBQ ATAC mapping parity passes on a tiny fixture.
- Default-compressed and level-0 CBQ both pass.
- The old CBQ-to-FASTQ adapter remains available only as a compatibility test
  oracle or fallback harness, clearly named as such.

## Phase 4: Smoke And Regression Tests

Add or update tests in STAR-suite:

- `tests/run_cbq_chromap_atac_e2e_smoke.sh`
  - synthetic genome;
  - synthetic ATAC R1/R2 plus barcode input;
  - FASTQ baseline;
  - CBQ in-memory provider path;
  - compare fragment tuples, barcode payloads, and summary metrics.
- Extend `tests/run_cbq_e2e_module_regression.sh` to call the new ATAC E2E
  smoke instead of the compatibility materialization smoke once stable.
- Keep the existing compatibility adapter smoke under a name that makes its
  status clear, for example `run_cbq_chromap_fastq_compat_smoke.sh`.
- Register the new ATAC CBQ E2E case in
  `tests/production_module_regression_manifest.tsv`.

Add or update tests in Chromap-suite:

- CLI ATAC vs libchromap ATAC parity on the PBMC 3K 100K fixture.
- FASTQ provider vs legacy loader parity while the legacy path still exists.
- Provider-fed synthetic ATAC fixture with barcode correction enabled.
- Negative tests for mismatched R1/R2/barcode counts.

Production-scale validation:

- Run the PBMC 3K 100K multiome fixture first.
- Then run the existing JAX multiome one-lane smoke.
- Only after those pass, run a full production sample.

## Phase 5: Documentation And MCP/Recipes

Update STAR-suite docs:

- `docs/EXPERIMENTAL_BINSEQ_INPUT.md`
- `docs/RUNBOOK_BINSEQ_INPUT_CONTRACT.md`
- `docs/LIBCHROMAP_CONTRACT.md`
- `docs/RUNBOOK_MULTIOME_MEX_MUDATA_20260516.md`
- `tests/ARTIFACTS.md`

Update Morphic recipes and MCP workflows after the code path is real:

- add CBQ ATAC input parameters alongside existing FASTQ ATAC parameters;
- mark Y/noY FASTQ emission as unsupported for CBQ unless implemented later;
- document that CBQ ATAC production means in-memory provider, not FASTQ
  materialization;
- expose smoke/regression commands for agents.

## Risks And Guardrails

- Do not implement a second CBQ reader in Chromap-suite.
- Do not use temporary FASTQ as the default production path.
- Do not expand scope to non-ATAC modalities while doing this integration.
- Do not break existing FASTQ ATAC production behavior.
- Watch read-id stability: Chromap duplicate removal, sorting, and output
  ordering rely on deterministic read IDs.
- Watch barcode range semantics: ATAC barcode extraction can use offset and
  reverse-complement behavior, so provider-filled batches must pass through the
  same effective range logic as FASTQ input.
- Watch lifetime: decoded CBQ spans may be borrowed from block buffers; either
  keep buffers alive through Chromap consumption or copy only into Chromap's
  existing batch storage.

## Suggested Implementation Order

1. In Chromap-suite, make ATAC CLI dispatch call `libchromap`.
2. Add FASTQ-backed ATAC provider with no output change.
3. Add an in-memory provider test harness inside Chromap-suite.
4. In STAR-suite, add the CBQ-to-provider adapter behind an explicit
   experimental flag or source enum.
5. Replace the current Chromap CBQ compatibility smoke with true in-memory
   ATAC parity.
6. Run PBMC 3K 100K and JAX one-lane multiome smokes.
7. Promote docs/MCP/recipes after test evidence is recorded.

## Completion Criteria

This runbook is complete when:

- ATAC FASTQ and ATAC CBQ in-memory paths share libchromap dispatch.
- STAR-suite can run Chromap ATAC from CBQ without materializing FASTQ.
- Synthetic and PBMC 3K 100K parity pass.
- JAX multiome one-lane smoke passes.
- Documentation and MCP/recipe surfaces clearly distinguish FASTQ production,
  CBQ production, and FASTQ compatibility shims.
