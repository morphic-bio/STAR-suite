# Hand-over to Multiomics Suite: STAR Suite 1.10.0

STAR Suite 1.10.0 no longer contains the Chromap integration. This document
lists every file and symbol removed from STAR, with the last STAR commit that
holds it, so that Multiomics Suite (step M1 of the single-binary design) can
copy them verbatim. Nothing listed here was rewritten on the way out.

- Removal commit on `dev-release-v1.10.0`: `826e019`
  (`826e01969ad0ec1154d6788baf9c4b28a042578f`). Its parent `90bb0e0`
  (`90bb0e0c150d5c02d68328ae19c34679d7c30f84`) is v1.9.5.a plus the two
  DOGMA-plex feature commits, and still holds every removed item.
- Host interface commit: `f8d2dae`; see `docs/HOST_API.md`, which also maps
  each STAR dependency of the moved code to the host API.
- Reference composition for byte identity (design decision 5): STAR
  `6e83853` (`6e8385375942d59b25b90b3aa9db0a87b9fefc7a`, branch
  `feat/atac-sidecar-without-bam`, on v1.9.5).

## 1. Which commit to copy from

Copy from **`6e83853`**. Six integration files are newer there than on the
1.10 line: they carry `--chromapAtacBarcodeSampleLimit` (`52c1ec8`, via merge
`fcfae2e`) and `--chromapAtacOutputFormat sidecar` (`6e83853`), which the
DOGMA-plex reference uses. The 1.10 line never received those two commits,
so on it the integration has 57 parameters, not 58. The two columns below give
both histories; "Same bytes" says whether the file is identical at `90bb0e0`
and `6e83853`.

```bash
STAR=/path/to/STAR-suite
git -C "$STAR" show 6e83853:core/features/libchromap_contract/src/star_chromap_contract.cpp > native/src/atac_contract.cpp
# or every removed path at once, preserving the tree:
git -C "$STAR" archive 6e83853 core/features/libchromap_contract \
    core/legacy/source/star_chromap_orchestration.cpp core/legacy/source/star_chromap_orchestration.h | tar -x -C staging/
```

The moved code includes STAR-internal headers (`Parameters.h`,
`ThreadControl.h`, `GlobalVariables.h`, `TimeFunctions.h`,
`SaturationPermitController.h`). Against STAR 1.10 these become
`host/StarHost.h` and `host/SaturationPermitController.h`, with the renames in
section 4. Keep the STAR-suite MIT notice on each moved file and record the
source commit in its header, as the design asks.

## 2. Removed files (45)

| Removed path | Last commit, 1.10 line | Last commit at `6e83853` | Same bytes | Lines | Destination (design 3.1) |
|---|---|---|---|---|---|
| `core/legacy/source/star_chromap_orchestration.cpp` | `012617d` | `6e83853` | **no** | 1697 | `native/src/atac_host.cpp` (re-pointed to the host API) |
| `core/legacy/source/star_chromap_orchestration.h` | `e2e68d4` | `e2e68d4` | yes | 42 | Multiomics-internal header |
| `core/legacy/source/star_chromap_orchestration_stub.cpp` | `3af52d5` | `3af52d5` | yes | 58 | none (the `WITH_CHROMAP=0` stub; not needed) |
| `core/features/libchromap_contract/.gitignore` | `9352462` | `9352462` | yes | 3 | reference only |
| `core/features/libchromap_contract/Makefile` | `6c5532b` | `bc84ea3` | **no** | 83 | reference for the CMake build |
| `core/features/libchromap_contract/include/multiome_atac_peak_mex.h` | `6c1107e` | `6c1107e` | yes | 60 | `native/` peak-step header |
| `core/features/libchromap_contract/include/star_chromap_contract.h` | `0a7a4f4` | `6e83853` | **no** | 159 | `native/` contract header |
| `core/features/libchromap_contract/src/multiome_atac_peak_mex.cpp` | `bc84ea3` | `bc84ea3` | yes | 1697 | `native/src/atac_peak_mex.cpp` |
| `core/features/libchromap_contract/src/star_chromap_contract.cpp` | `0a7a4f4` | `6e83853` | **no** | 509 | `native/src/atac_contract.cpp` |
| `core/features/libchromap_contract/tools/star_libchromap_contract_runner.cpp` | `0a7a4f4` | `6e83853` | **no** | 296 | subcommand `atac-contract-runner` |
| `core/features/libchromap_contract/tools/star_multiome_atac_peak_mex.cpp` | `6c1107e` | `6c1107e` | yes | 245 | subcommand `atac-peak-mex` |
| `core/features/libscrna/include/libscrna/AtacEvidenceFromPeaks.h` | `2a4f20f` | `2a4f20f` | yes | 27 | `native/` (AEV1 evidence reader) |
| `core/features/libscrna/src/AtacEvidenceFromPeaks.cc` | `2a4f20f` | `2a4f20f` | yes | 454 | `native/` (AEV1 evidence reader) |
| `core/features/libscrna/tools/scrna_build_atac_evidence.cpp` | `8de086a` | `8de086a` | yes | 342 | `native/` multiome tools |
| `core/features/libscrna/tools/scrna_multiome_combine.cpp` | `8de086a` | `8de086a` | yes | 615 | `native/` multiome tools |
| `core/features/libscrna/tools/scrna_compare_multiome_calls.cpp` | `8de086a` | `8de086a` | yes | 451 | `native/` multiome tools |
| `core/features/libscrna/tools/scrna_arc_barcode_table.cpp` | `8de086a` | `8de086a` | yes | 213 | `native/` multiome tools |
| `core/features/libscrna/tools/scrna_build_gex_evidence.cpp` | `8de086a` | `8de086a` | yes | 437 | `native/` multiome tools |
| `core/legacy/source/input/CbqChromapAdapter.h` | `4e00aea` | `4e00aea` | yes | 31 | Multiomics (CBQ support, section 5) |
| `core/legacy/source/input/CbqChromapAdapter.cpp` | `4e00aea` | `4e00aea` | yes | 225 | Multiomics (CBQ support, section 5) |
| `core/legacy/source/input/cbq_chromap_adapter_harness.cpp` | `4e00aea` | `4e00aea` | yes | 103 | Multiomics (CBQ support, section 5) |
| `tests/run_cbq_chromap_adapter_smoke.sh` | `4e00aea` | `4e00aea` | yes | 220 | Multiomics (CBQ support, section 5) |
| `tests/run_star_chromap_macs3_lowmem_smoke_100k.sh` | `7afa8a4` | `7afa8a4` | yes | 119 | Multiomics `tests/` |
| `tests/run_star_libchromap_cbq_contract_smoke.sh` | `2b01ea1` | `2b01ea1` | yes | 127 | Multiomics `tests/` |
| `tests/test_star_multiome_atac_peak_mex_signac_profile.sh` | `6c1107e` | `6c1107e` | yes | 66 | Multiomics `tests/` |
| `tests/test_catatac_trimodal_downsample_smoke.sh` | `81fb08e` | `81fb08e` | yes | 327 | Multiomics `tests/` |
| `tests/catatac_trimodal_downsample_verify.py` | `63545a7` | `63545a7` | yes | 267 | Multiomics `tests/` |
| `tests/multi_feature/test_hiv_dogma_four_arm_table_smoke.sh` | `e9589d7` | `e9589d7` | yes | 267 | Multiomics `tests/` |
| `scripts/build_atac_peak_matrix_from_fragments.py` | `847d4b9` | `847d4b9` | yes | 223 | Multiomics (section 6) |
| `scripts/build_multiome_mudata.py` | `6613482` | `6613482` | yes | 527 | Multiomics (section 6) |
| `scripts/convert_jax_multiome01_mt_adaptive_once.sh` | `8de086a` | `8de086a` | yes | 293 | Multiomics (section 6) |
| `scripts/normalize_multiome_atac_barcode_fastq.py` | `847d4b9` | `847d4b9` | yes | 113 | Multiomics (section 6) |
| `scripts/run_jax_multiome01_production.sh` | `8424098` | `8424098` | yes | 23 | Multiomics (section 6) |
| `scripts/run_multiome_cell_call_external_gex_from_arc.sh` | `8de086a` | `8de086a` | yes | 243 | Multiomics (section 6) |
| `scripts/run_multiome_cell_call_from_arc.sh` | `8de086a` | `8de086a` | yes | 240 | Multiomics (section 6) |
| `scripts/run_multiome_cell_call_harness_from_arc.sh` | `8de086a` | `8de086a` | yes | 164 | Multiomics (section 6) |
| `scripts/run_multiome_mudata_smoke.sh` | `8424098` | `8424098` | yes | 23 | Multiomics (section 6) |
| `scripts/run_remote_multiome_post_mex_rsync.sh` | `8424098` | `8424098` | yes | 23 | Multiomics (section 6) |
| `scripts/run_star_multiome_lane_smoke.sh` | `8424098` | `8424098` | yes | 23 | Multiomics (section 6) |
| `scripts/upload_jax_multiome01_large_files_globus.sh` | `8424098` | `8424098` | yes | 23 | Multiomics (section 6) |
| `mcp_server/workflows/morphic_multiome.yaml` | `89a008d` | `89a008d` | yes | 336 | Multiomics MCP |
| `docs/LIBCHROMAP_CONTRACT.md` | `6c5532b` | `52c1ec8` | **no** | 337 | Multiomics docs |
| `docs/RUNBOOK_CHROMAP_ATAC_CBQ_IN_MEMORY.md` | `2b01ea1` | `2b01ea1` | yes | 290 | Multiomics docs |
| `docs/HANDOFF_STAR_LIBCHROMAP_MACS3_INTEGRATION_20260425.md` | `0a7a4f4` | `0a7a4f4` | yes | 350 | Multiomics docs |
| `docs/RUNBOOK_MULTIOME_MEX_MUDATA_20260516.md` | `bc84ea3` | `bc84ea3` | yes | 1114 | Multiomics docs |

## 3. Code removed from files that stay in STAR

Line numbers are at `90bb0e0` / `6e83853`. The integration code in these
files is identical at both revisions except where noted.

| File | Lines | What | 1.10 replacement |
|---|---|---|---|
| `core/legacy/source/STAR.cpp` | 65 / 65 | `#include "star_chromap_orchestration.h"` | `host/StarHost.h` |
| | 1334-1346 / 1334-1346 | `preflightStarChromapAtacIfEnabled` (exit 102) and `StarChromapAtacAsyncRun chromapAtacAsyncRun; startStarChromapAtacIfEnabled` (exit 103), with their error texts | `hooks->preflight`, `hooks->start` |
| | 3580-3585 / 3580-3585 | `runStarChromapAtacIfEnabled` (exit 103) | `hooks->finish` |
| | 3604-3610 / 3604-3610 | `atacInUse=`/`atacRetainedLease=` log fields; `fixedPoolIncomplete` from `P.dynamicThreadAtacController == 2` | label from `externalLabel`; `hooks->requiresFullPoolAtExit` |
| `core/legacy/source/mapThreadsSpawn.cpp` | 14 / 14 | `#include "SaturationPermitController.h"` | not needed in STAR |
| | 261-263 / 261-263 | Flex fused guard: `P.chromapAtac.enabled != 0 \|\| !unsetToken(P.multiomeAtacPeakMex.inlineMode)` rejects with "Chromap/multiome ATAC output is enabled" | `hooks->externalActive` ("host external work is enabled") |
| | 925-933 / 925-933 | pool = `runThreadN + max(0, chromapAtac.threads)` when `chromapAtac.enabled == 1` | `hooks->extraPermitThreads` |
| | 939 / 939 | `configuredAtacFloor = max(0, P.dynamicThreadAtacFloor)` | `hooks->initialFloors` (floors[2]) |
| | 949-973 / 949-973 | mode-2 initial floors: `SaturationPermitController` initial decision when `interfaceEnabled && chromapAtac.enabled == 1 && dynamicThreadAtacController == 2`, fatal (exit 1) below 3 permits with a FEATURE estimate | `hooks->initialFloors` |
| | 1007 / 1007 | `atacController=` field of the "Dynamic thread interface enabled" line | the host logs its controller mode itself |
| `core/legacy/source/Parameters.h` | 66, 73, 80, 87 | `dynamicThreadTelemetryIntervalSec`, `dynamicThreadAtacFloor`, `dynamicThreadAtacController`, `dynamicThreadAtacWorkEstimate` | host options |
| | 819-891 / 819-892 | `struct {...} chromapAtac` (6e83853 adds `barcodeSampleLimit`) and `struct {...} multiomeAtacPeakMex` | host options |
| `core/legacy/source/Parameters.cpp` | 167, 169, 171, 174 | registration of the four permit parameters | host parser |
| | 795-847 / 795-848 | registration of the `chromapAtac*` and `multiomeAtac*` parameters | host parser |
| | 1020-1231 / 1021-1236 | default resets (`if (p->nameString == "..." && p->inputLevel < 0)`), including the only definition of the `chromapAtacBarcodeTranslateFromFirst` default | host defaults |
| | 2029-2051 / 2034-2056 | checks: telemetry interval >= 0; floors >= 0 (the message named `--dynamicThread{Map,Atac,Feature}Floor`); controller in 0..2; mode 2 needs telemetry and a positive interval | host `preflight` (STAR keeps the Map/Feature floor check) |
| | 2067-2073 / 2072-2078 | BGZF hierarchy refused with `dynamicThreadAtacController != 0` or `chromapAtac.enabled == 1` | host `preflight`, `RunView::bgzfHierarchy` |
| | 2218-2228 / 2223-2233 | mode 2 refused with an applying `--dynamicThreadPfControllerMode` (`pfControllerAppliesUpdates`) | host `preflight`, `RunView::pfControllerMode` |
| `core/legacy/source/parametersDefault` (+ `.xxd`) | 51-60, 70-73, 78-91, 102-105 / same | help text of the four permit parameters; ATAC wording in the Map floor, work-estimate and BGZF hierarchy entries | host `--help` section |
| | 2406-2597 / 2412-2609 | the `chromapAtac*` and `multiomeAtac*` help block | host `--help` section |
| `core/legacy/source/ReadAlignChunk_processChunks.cpp` | 1423-1426 | comment on the wider `chromapAtac` pool | comment only |
| `core/legacy/source/Makefile` | 49-106 / 49-97 | `WITH_CHROMAP` switch, `CHROMAP_SUITE_DIR`, `CHROMAP_LINK_LINE`, HTSlib pairing | `HTSLIB=bundled\|external` |
| | 262 / 253 | `$(CHROMAP_ORCH_OBJ)` in `OBJECTS` | none |
| | 451-458, 761-762 / 442-449, 752-753 | `cbq-chromap-adapter-harness`, `CbqChromapAdapter.o` rules | section 5 |
| | 806-839 / 797-805 | `chromap-unsupported-target`; at `90bb0e0` also `check-core-dependencies` (HTSlib and Chromap checkout checks) | `check-external-htslib` |
| | 902-912 / 870-880 | FORCE rules building `libstar_chromap_contract.a`, `libchromap.a`, `librapidmacs` | Multiomics build |
| `Makefile`, `build/*.mk` | | `star-libchromap-contract` target, `LIBCHROMAP_CONTRACT_DIR`, `DEV_RELEASE_CHECK_SCRIPTS` entry, `WITH_CHROMAP=0` arguments | none |
| `core/features/libscrna/Makefile` | | `CC_SRCS := src/AtacEvidenceFromPeaks.cc` and the five multiome tool targets | Multiomics build |
| `tests/production_module_regression_manifest.tsv` | 16, 30 / 15, 29 | rows `multiome cbq-chromap-adapter` and `multiome chromap-macs3-100k` | Multiomics test manifest |
| `tests/run_cbq_e2e_module_regression.sh` | | `cbq_chromap_adapter` case and `RUN_CHROMAP_MAPPING_SMOKE` | section 5 |
| `tests/run_partial_make_regression.sh` | | `core-without-chromap` and `chromap-core` cases | none |
| `tests/test_htslib_build_discovery.py` | | Chromap checkout fixtures and two Chromap tests | rewritten for `HTSLIB=` |
| `tests/test_ambient_fdr_anndata_mudata.py` (now `test_ambient_fdr_anndata.py`), last commit `7da2e6a` | 126-161 | MuData half: builds `out.h5mu` with `scripts/build_multiome_mudata.py` from an ATAC peak MEX and checks guide ambient-FDR columns in `mdata.obs` | Multiomics test next to its MuData builder |
| `mcp_server/config.yaml` | 196-202 | `morphic_multiome` workflow entry | Multiomics MCP |
| `mcp_server/launchpad/static/app.js`, `index.html`, `mcp_server/README.md` | | `morphic_multiome` in the default workflow list and docs | Multiomics MCP |
| `docker/Dockerfile`, `debian/rules`, CI | | `STAR_WITH_CHROMAP`, `core-portable` | `make core` |

## 4. Symbols

Removed from STAR (copy with the files in section 2):

- Orchestration (`star_chromap_orchestration.{h,cpp}`): `StarChromapAtacAsyncRun`
  (and its `Impl`: worker thread, permit-telemetry sampler, mode-1 drain-time
  and mode-2 saturation controllers), `preflightStarChromapAtacIfEnabled`,
  `startStarChromapAtacIfEnabled`, `runStarChromapAtacIfEnabled`,
  `validateAndBuildConfig`, `runInlinePeakMexIfEnabled`,
  `runEvidenceFromPeaksIfEnabled`, `parameterInputLevel`, the permit shims
  `chromapStarDynamicPermitAcquire` / `chromapStarDynamicPermitRelease` with
  `kAtacPermitHookContext`, and helpers (`parsePeakCallMode`,
  `resolveInlinePeakCallMode`, `isConcurrentStartMode`,
  `deriveSecondaryOutputPath`, ...).
- Contract (`star_chromap_contract.h`, namespace `star::multiome`):
  `ChromapOutputFormat`, `Tn5ShiftMode`, `ChromapMacs3FragPeaksSource`,
  `ChromapMacs3FragThresholdMode`, `ChromapInputFormat`,
  `ChromapContractStatus`, `ChromapPermitAcquireFn`,
  `ChromapPermitReleaseFn`, `ChromapPermitHooks`, `ChromapAtacConfig`,
  `ChromapAtacResult`, `runChromapAtac`, `chromapContractStatusName`.
- Peak/matrix step (`multiome_atac_peak_mex.h`, `star::multiome`):
  `MultiomeAtacPeakMexThresholdMode`, `MultiomeAtacPeakCallMode`,
  `MultiomeAtacPeakMexArgs`, `RunMultiomeAtacPeakMex`.
- libscrna (`libscrna/AtacEvidenceFromPeaks.h`, `libscrna::atac`):
  `RunAtacEvidenceFromBinary`.
- CBQ adapter (`input/CbqChromapAdapter.h`, `star::input`):
  `CbqChromapAdapterOptions`, `CbqChromapFastqPaths`,
  `materialize_cbq_chromap_fastqs`.

Renamed in STAR, used by the moved code (re-point per `docs/HOST_API.md`):

| 1.9.5.a | 1.10.0 |
|---|---|
| `ThreadControl::PermitDomain::ATAC` | `ThreadControl::PermitDomain::EXTERNAL`; hosts use `star::host::Domain::External` |
| `ThreadControl::MapPermitSnapshot` (`atacDomain`) | `star::host::PermitSnapshot` (`externalDomain`); `ThreadControl::MapPermitSnapshot` is an alias |
| `SaturationPermitController.h`, `star::multiome::SaturationPermitController` | `host/SaturationPermitController.h`, `star::permits::SaturationPermitController` |
| controller `Domain::ATAC`, `Phase::PROBE_ATAC`, `Reason::ATAC_PROBE_COMPLETE`, `ATAC_ETA_LATE`, `ATAC_COMPLETE` | `EXTERNAL`, `PROBE_EXTERNAL`, `EXTERNAL_PROBE_COMPLETE`, `EXTERNAL_ETA_LATE`, `EXTERNAL_COMPLETE` |
| controller fields `WorkEstimates::atac`, `atacOccupancy`, `atacUnitsDelta`, `atacInUse`, `atacWaiters`, `atacEtaSec`, `atacEstimateComplete`, `atacFloor`, `atacSaturation`, `atacSaturationKnown` | `external`, `externalOccupancy`, ... (same suffixes) |
| names `"atac"`, `"probe-atac"`, `"atac-probe-complete"`, `"atac-eta-late"`, `"atac-complete"` | `"external"`, ...; `domainName/phaseName/reasonName(x, "atac")` print the old strings |
| log fields `floors(map/feature/atac)`, `atacState(...)`, `atacInUse=`, stall `domain=atac` | the host's `externalLabel`; with `externalLabel = "atac"` the lines read as in 1.9.5 |

The sampler's own log lines (`[ATAC permit telemetry]`, `[ATAC saturation
controller]`, `[ATAC drain-time controller]`, `ATAC permit telemetry:`) are in
the moved orchestration and are unchanged by the move.

## 5. CBQ support (author decision, 28 September)

`CbqChromapAdapter` is moved, not deleted. Copy these verbatim; CBQ through
this adapter is to be tested at a later date. Chromap Suite keeps its own
native CBQ reader.

| Item | Last commit | Notes |
|---|---|---|
| `core/legacy/source/input/CbqChromapAdapter.h`, `CbqChromapAdapter.cpp` | `4e00aea` | CBQ to synchronized R1/R2/barcode FASTQs for Chromap |
| `core/legacy/source/input/cbq_chromap_adapter_harness.cpp` | `4e00aea` | standalone harness |
| Makefile rules `cbq-chromap-adapter-harness` and `CbqChromapAdapter.o` | `90bb0e0` lines 451-458, 761-762 | builds with STAR's `CbqInputModule.o` (`input/CbqInputModule.{h,cpp}`, `input/InputContract.h`), which stay in STAR; `libstar_suite.a` contains `CbqInputModule.o` |
| `tests/run_cbq_chromap_adapter_smoke.sh` | `4e00aea` | payload parity; optional tiny mapping when `CHROMAP_BIN` exists |
| manifest row `multiome cbq-chromap-adapter`; `cbq_chromap_adapter` case of `tests/run_cbq_e2e_module_regression.sh` | `90bb0e0` | |

## 6. Scripts

Checked against Multiomics Suite `master` (`57b868f`) on 28 September. Eleven
of the twelve have a current copy there (identical, or changed later in
Multiomics); one has none and must be copied.

| STAR script (removed) | Last STAR commit | Multiomics copy | State |
|---|---|---|---|
| `scripts/build_atac_peak_matrix_from_fragments.py` | `847d4b9` | `scripts/experimental/build_atac_peak_matrix_from_fragments.py` | identical |
| `scripts/build_multiome_mudata.py` | `6613482` (2026-06-13) | `scripts/build_multiome_mudata.py` | Multiomics newer (`8df8ac9`, 2026-09-26) |
| `scripts/convert_jax_multiome01_mt_adaptive_once.sh` | `8de086a` | none | **copy** |
| `scripts/normalize_multiome_atac_barcode_fastq.py` | `847d4b9` | `scripts/normalize_multiome_atac_barcode_fastq.py` | identical |
| `scripts/run_jax_multiome01_production.sh` | `8424098` | `recipes/run_jax_multiome01_production.sh` | Multiomics newer; STAR's was a 23-line launcher |
| `scripts/run_multiome_cell_call_external_gex_from_arc.sh` | `8de086a` | `scripts/experimental/…` | identical |
| `scripts/run_multiome_cell_call_from_arc.sh` | `8de086a` | `scripts/experimental/…` | identical |
| `scripts/run_multiome_cell_call_harness_from_arc.sh` | `8de086a` | `scripts/experimental/…` | identical |
| `scripts/run_multiome_mudata_smoke.sh` | `8424098` | `recipes/run_multiome_mudata_smoke.sh` | Multiomics newer; STAR's was a launcher |
| `scripts/run_remote_multiome_post_mex_rsync.sh` | `8424098` | `recipes/run_remote_multiome_post_mex_rsync.sh` | Multiomics newer; STAR's was a launcher |
| `scripts/run_star_multiome_lane_smoke.sh` | `8424098` | `recipes/run_star_multiome_lane_smoke.sh` | Multiomics newer; STAR's was a launcher |
| `scripts/upload_jax_multiome01_large_files_globus.sh` | `8424098` | `recipes/upload_jax_multiome01_large_files_globus.sh` | Multiomics newer; STAR's was a launcher |

The multiome tests of section 2 (`run_star_chromap_macs3_lowmem_smoke_100k.sh`,
`run_star_libchromap_cbq_contract_smoke.sh`,
`test_star_multiome_atac_peak_mex_signac_profile.sh`,
`test_catatac_trimodal_downsample_smoke.sh` with
`catatac_trimodal_downsample_verify.py`,
`multi_feature/test_hiv_dogma_four_arm_table_smoke.sh`) have no Multiomics
copies; Multiomics `recipes/run_catatac_trimodal_downsample_smoke.sh`,
`recipes/run_catatac_trimodal_full_benchmark.sh`,
`recipes/run_hiv_dogma_four_arm_downsample_smoke.sh` and
`tests/test_downsample_smoke_recipes.py` call them in the STAR tree today and
must point at the copies (step M3) before they are used with STAR 1.10.

## 7. Parameters (58)

In STAR registration order at `6e83853`. Names, types and defaults are to stay
unchanged in Multiomics (design 3.4). Help text: `parametersDefault` lines in
section 3.

| # | Parameter | Type | Default | Default defined in | On 1.9.5.a base |
|---|---|---|---|---|---|
| 1 | `dynamicThreadTelemetryIntervalSec` | int | `10` | parametersDefault | yes |
| 2 | `dynamicThreadAtacFloor` | int | `0` | parametersDefault | yes |
| 3 | `dynamicThreadAtacController` | int | `0` | parametersDefault | yes |
| 4 | `dynamicThreadAtacWorkEstimate` | uint64 | `0` | parametersDefault | yes |
| 5 | `chromapAtacEnable` | int | `0` | parametersDefault | yes |
| 6 | `chromapAtacReferenceFasta` | string | `-` | parametersDefault | yes |
| 7 | `chromapAtacIndex` | string | `-` | parametersDefault | yes |
| 8 | `chromapAtacInputFormat` | string | `fastq` | parametersDefault | yes |
| 9 | `chromapAtacRead1` | string | `-` | parametersDefault | yes |
| 10 | `chromapAtacRead2` | string | `-` | parametersDefault | yes |
| 11 | `chromapAtacBarcode` | string | `-` | parametersDefault | yes |
| 12 | `chromapAtacReadPairCbq` | string | `-` | parametersDefault | yes |
| 13 | `chromapAtacBarcodeCbq` | string | `-` | parametersDefault | yes |
| 14 | `chromapAtacReadFormat` | string | `-` | parametersDefault | yes |
| 15 | `chromapAtacBarcodeWhitelist` | string | `-` | parametersDefault | yes |
| 16 | `chromapAtacBarcodeSampleLimit` | uint64 | `20000000` | parametersDefault | no (added by 52c1ec8) |
| 17 | `chromapAtacBarcodeTranslate` | string | `-` | parametersDefault | yes |
| 18 | `chromapAtacBarcodeTranslateFromFirst` | int | `0` | reset in Parameters.cpp | yes |
| 19 | `chromapAtacOutputFragments` | string | `-` | parametersDefault | yes |
| 20 | `chromapAtacSecondaryFragments` | string | `-` | parametersDefault | yes |
| 21 | `chromapAtacOutputFormat` | string | `BED` | parametersDefault | yes |
| 22 | `chromapAtacSummary` | string | `-` | parametersDefault | yes |
| 23 | `chromapAtacTempDir` | string | `-` | parametersDefault | yes |
| 24 | `chromapAtacThreads` | int | `1` | parametersDefault | yes |
| 25 | `chromapAtacHtsThreads` | int | `0` | parametersDefault | yes |
| 26 | `chromapAtacSortBam` | int | `0` | parametersDefault | yes |
| 27 | `chromapAtacWriteIndex` | int | `0` | parametersDefault | yes |
| 28 | `chromapAtacSortBamRam` | uint64 | `8589934592` | parametersDefault | yes |
| 29 | `chromapAtacEmitNoYBam` | int | `0` | parametersDefault | yes |
| 30 | `chromapAtacEmitYBam` | int | `0` | parametersDefault | yes |
| 31 | `chromapAtacNoYOutput` | string | `-` | parametersDefault | yes |
| 32 | `chromapAtacYOutput` | string | `-` | parametersDefault | yes |
| 33 | `chromapAtacLowMem` | int | `0` | parametersDefault | yes |
| 34 | `chromapAtacLowMemRam` | uint64 | `0` | parametersDefault | yes |
| 35 | `chromapAtacCallMacs3FragPeaks` | int | `0` | parametersDefault | yes |
| 36 | `chromapAtacMacs3FragPeaksOutput` | string | `-` | parametersDefault | yes |
| 37 | `chromapAtacMacs3FragSummitsOutput` | string | `-` | parametersDefault | yes |
| 38 | `chromapAtacMacs3FragKeepIntermediates` | string | `-` | parametersDefault | yes |
| 39 | `chromapAtacMacs3FragPvalue` | double | `1e-5` | parametersDefault | yes |
| 40 | `chromapAtacMacs3FragQvalue` | double | `0` | parametersDefault | yes |
| 41 | `chromapAtacMacs3FragMinLength` | int | `200` | parametersDefault | yes |
| 42 | `chromapAtacMacs3FragMaxGap` | int | `30` | parametersDefault | yes |
| 43 | `chromapAtacMacs3FragUint8Counts` | int | `1` | parametersDefault | yes |
| 44 | `chromapAtacMacs3FragLowMem` | int | `0` | parametersDefault | yes |
| 45 | `chromapAtacEvidenceFromPeaksOutput` | string | `-` | parametersDefault | yes |
| 46 | `chromapAtacTn5ShiftMode` | string | `classical` | parametersDefault | yes |
| 47 | `chromapAtacStartMode` | string | `postMapping` | parametersDefault | yes |
| 48 | `multiomeAtacPeakMexInline` | string | `no` | parametersDefault | yes |
| 49 | `multiomeAtacPeakBarcodeTranslate` | string | `-` | parametersDefault | yes |
| 50 | `multiomeAtacPeakBarcodeTranslateFromFirst` | string | `yes` | parametersDefault | yes |
| 51 | `multiomeAtacPeakMetricsTsv` | string | `-` | parametersDefault | yes |
| 52 | `multiomeAtacPeakMexOutDir` | string | `-` | parametersDefault | yes |
| 53 | `multiomeAtacPeakNarrowPeak` | string | `-` | parametersDefault | yes |
| 54 | `multiomeAtacPeakSummits` | string | `-` | parametersDefault | yes |
| 55 | `multiomeAtacPeakCallMode` | string | `frag` | parametersDefault | yes |
| 56 | `multiomeAtacPeakMacsProfile` | string | `-` | parametersDefault | yes |
| 57 | `multiomeAtacPeakThreads` | int | `0` | parametersDefault | yes |
| 58 | `multiomeAtacPeakMaxBarcodes` | uint64 | `0` | parametersDefault | yes |

Explicit-set tracking: the orchestration's `parameterInputLevel()` read STAR's
input level for `chromapAtacMacs3FragPvalue`, `chromapAtacMacs3FragQvalue`,
`chromapAtacMacs3FragMinLength`, `chromapAtacMacs3FragMaxGap` and
`multiomeAtacPeakCallMode`. The host's `parameter` callback receives that level
(2 for the command line, 5 and above for parameter files); a parameter never
delivered is at its default.

## 8. Not moved

- CAT-ATAC guide arm (`catatac_guide` layout in `PfMultiConfig`, ATAC barcode
  namespace in `PfMultiProcess.cpp`, `catatac_crispri_guide_capture.csv`,
  `test_catatac_split_read.c` and the `test_catatac_guide_arm*` smokes): a
  feature-arm feature, stays in STAR (author decision 9).
- Permit allocator (`ThreadControl`) and the saturation controller: stay in
  STAR behind the host interface.
- HTSlib header sites: re-keyed from `WITH_CHROMAP` to `STAR_EXTERNAL_HTSLIB`
  (`HTSLIB=external`).
- libscrna EmptyDrops, OrdMag, occupancy, C API and `scrna_simpleed`.
- `share/star-suite/` recipe snapshot: still pinned to STAR-suite-recipes
  `1b12325` (six multiome workflows). Re-pin when those recipes point at the
  Multiomics binary (design 6.4); not done in 1.10.0-rc1.
- `core/legacy/source/STAR.release`, a tracked 1.0.0-era binary used by the
  UCSF regression row, contains old Chromap strings; it is not built or linked.
