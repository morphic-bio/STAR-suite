# Cross-module performance regression harnesses

See [the runbook](../docs/RUNBOOK_CROSS_MODULE_PERFORMANCE_20260910.md) and
[measured results](../docs/benchmarks/CROSS_MODULE_PERFORMANCE_20260910.md).
These tests compare execution/storage changes without changing model parameters.
Use fresh artifact directories and the host's identical-execution policy.

- `pf_bgzf_harness.cpp` and `run_pf_bgzf_tests.py`: standalone PF decoding,
  pairing, assignment, shared permits and input failures. Build with
  `make -C core/features/process_features ../../../tests/pf_bgzf_harness`.
  The Python driver creates the small FASTQ/whitelist/feature fixture.
- `pf_direct_failure.cpp`: concurrent malformed batches (quality length, layout,
  extreme field length) must fail without final MEX and release permits before
  cleanup; a maximum-sized newline-terminated sequence with default qualities
  must succeed. Compile against `libprocess_features.a` and `libscrna.a`, using
  `-Icore/features/process_features/include -std=c++17 -fopenmp -lpthread -lz -lhts`.
  Arguments: the PF fixture's `inputs` directory and an existing fresh output
  directory. Requires `R1.gz`, `R2.gz`, `whitelist.txt`, `features.csv`.
- `emptydrops/cross_module_caller.cpp`: arguments are MEX directory, output TSV,
  MC workers, `split|strided`, and simulation count. Uses four logical bootstrap
  streams in both arms. Compile with the libscrna include/archive and pthread;
  `-DCALLER_BASELINE` builds against the frozen pre-change C API. The changed
  arm also checks wide-offset rejection and weighted scheduler exception cleanup.
- `slam/cross_module_solver.cpp`: serializes both solver results in hexadecimal
  floating point, including likelihood, iterations and convergence. Link
  `SlamSolver.cpp` and `slam/source/libem/slam_vb_overdisp.cpp`; compile both the
  frozen and cached implementations with identical optimization flags.
- `slam/test_slam_fit_reuse.cpp`: build with
  `make -C core/legacy/source slam-fit-reuse-tests`; run
  `core/legacy/source/slam_fit_reuse_tests`. Covers both models, one/five workers,
  repeated writers, parameter changes and histogram mutation/merge.
- `transcriptvb/test_component_execution.cpp`: compile with
  `-Icore/features/vbem/source/libem -fopenmp` and link its `libem.a`. Covers empty
  and unsupported input, one/multiple components, more workers than components,
  global iterations/convergence, effective-length callbacks and numeric outputs.
- OCM native materializer: build `ocm-multi-unit-tests`. The existing
  `tests/test_ocm_multi_unit.cpp` accepts `OCM_TEST_THREADS` and
  `OCM_TEST_MEMORY` for a copied fixture with no pool filtered MEX. Use
  `OCM_TEST_FIXTURE_ROOT` and `OCM_TEST_LOG`; execute the `materialize` subcommand.
  Compare decoded matrices, namespaces, calls and diagnostics between arms.

The measured-workflow scripts, exact commands, frozen executables, input hashes,
comparison outputs and failed attempts are retained at the artifact location in
`tests/ARTIFACTS.md`. Integrated acceptance includes whole perturb guide calling,
real SLAM auto-trim/per-file reopening, core batch transitions, FLEX and CBQ.
Native BGZF deliberately falls back for TranscriptVB online model learning;
its integrated component gate uses deterministic one-worker mapping, while the
full saved EC comparison exercises eight-worker fitting without re-alignment.
