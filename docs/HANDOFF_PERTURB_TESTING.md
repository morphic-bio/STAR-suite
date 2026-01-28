# Handoff: MCP Server + Perturb Testing Focus

## Repo Status

- Branch: `master` updated and pushed.
- MCP server implementation merged and live under `mcp_server/`.
- No tests were run during the merge.
- Untracked artifacts remain (local): `mcp_server/__pycache__/`, `tests/*`, `core/features/vbem/test/nonflex_mex_smoke_output/`.
  - Do not commit these; see `tests/ARTIFACTS.md`.

## MCP Server Summary

Implemented FastMCP server with:
- Discovery tools (`list_datasets`, `list_test_suites`, `find_docs`, `find_tests`).
- Config reload (`reload_config`) with executor settings update.
- Preflight validation (trusted roots, fixtures, binaries, disk space, working_dir).
- Script execution queue with true concurrency (ThreadPoolExecutor).
- Output collection via `output_roots` (includes logs and out_dir/default).
- Docker + docker-compose + `config.docker.yaml`.
- Unit/integration tests under `mcp_server/tests/`.

Key locations:
- `mcp_server/app.py` (tools + auth)
- `mcp_server/config.yaml` (default config, public discovery)
- `mcp_server/tools/preflight.py`
- `mcp_server/tools/executor.py`
- `mcp_server/README.md` (usage)

Known issues / follow-ups:
- Build commands run via bare `"make"/"cmake"/"ninja"` (PATH) even though preflight
  validated a trusted binary path. Consider using resolved binary path.
- `cancel()` only terminates the parent process, not the full process group.
  Consider killing the process group (mirroring timeout logic).
- Default output dir (`/tmp/mcp_runs/<run_id>`) is recorded but not created if
  `out_dir` is omitted; consider creating it.

## Perturb Module Testing Focus

Primary tests to run:
- `tests/test_cr_compat_crispr_calling.sh`
- `tests/run_a375_gex_features_cr_parity_genefull.sh`
- `core/features/process_features/tests/test_assignbarcodes_regression.sh`
- `tests/run_cbub_regression_test.sh`

Datasets:
- A375 data: `/storage/A375` (see `tests/ARTIFACTS.md` for specific runs).

Docs to reference:
- `docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md`
- `docs/feature_barcodes.md`
- `docs/TODO_crispr_feature_calling.md`
- `docs/UNSORTED_BAM_CBUB_IMPLEMENTATION_SUMMARY.md`
- `tests/nbem_sceptre_parity_results.md` (NB‑EM/Poisson‑EM/SCEPTRE‑parity results + test sets)

Artifact hygiene:
- Write new outputs to `tests/*_output*/` or `/tmp`.
- Add new output locations to `tests/ARTIFACTS.md`.
