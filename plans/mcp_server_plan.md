# MCP Server Plan (Internal Agent Development)

Goal: provide a lightweight, stable MCP service for agents to discover data,
run preflight checks, and execute standard test/run scripts. HPC/Temporal
executor integration is deferred to a later stage.

## Stage 0: Design and Scope

- Define server scope: discovery + preflight + local execution only.
- Decide trusted roots:
  - Repo: `/mnt/pikachu/STAR-suite`
  - Data: `/storage/*` (explicit allowlist)
  - Temp: `/tmp/*`
- Establish minimal endpoints (see Stage 1).
- Decide auth: token or IP allowlist (internal-only).

## Stage 1: Minimal FastMCP Server (No Docker)

Deliverable: `mcp_server/` with a single `app.py` runnable in venv.

Endpoints:
- `list_datasets()` -> known dataset roots from config.
- `list_test_suites()` -> module-aligned test suites + required fixtures.
- `find_docs(topic)` -> returns doc paths (from `docs/`).
- `find_tests(tag)` -> returns test scripts (from `tests/`).
- `preflight(run_config)` -> validates paths, references, flags.
- `run_script(name, args)` -> runs a whitelisted script with args.
- `collect_outputs(run_id)` -> returns output paths and logs.

Config:
- `mcp_server/config.yaml` with:
  - repo_root
  - dataset_roots
  - test_suites (per module)
  - allowed_scripts
  - allowed_flags / schema for preflight
  - artifact_log_path

Logging:
- Append run summaries to `plans/artifacts/mcp_runs_YYYYMMDD/`.
- Record outputs in `tests/ARTIFACTS.md` when applicable.

## Module Test Suites (initial mapping)

Maintain a module->tests map so agents can discover a minimal, known-good
set of scripts quickly.

- Core / Solo (non-Flex):
  - `tests/run_cbub_regression_test.sh`
  - `tests/run_sorted_bam_cbub_nonflex_test.sh`
  - `tests/run_unsorted_cbub_smoke_test.sh`
- Flex:
  - `tests/run_flex_smoke.sh`
  - `tests/run_flex_cbub_validation_test.sh`
  - `tests/run_flex_multisample_test.sh`
- CR-compat / CRISPR:
  - `tests/test_cr_compat_crispr_calling.sh`
  - `tests/run_a375_gex_features_cr_parity*.sh`
- Feature tools:
  - `core/features/process_features/tests/test_assignbarcodes_regression.sh`
  - `core/features/process_features/tests/test_pf_api.c` (build+run)
- SLAM:
  - `slam/tests/...` (add once confirmed)

Each suite should declare required fixtures and output locations; the MCP
server can expose this as `list_test_suites()`.

## Stage 2: Stabilize + Docs

- Add `mcp_server/README.md` (how to run, auth, examples).
- Add a short link in `README.md` and `AGENTS.md`.
- Add minimal tests for schema validation and preflight errors.

## Stage 3: Containerize (Recommended for Shared Server)

Deliverable: `mcp_server/Dockerfile` + `docker-compose.yml`.

- Bind mounts:
  - repo root
  - `/storage` (read-only if desired)
  - `/tmp`
- Environment-driven config path.
- Basic token auth or IP allowlist.

## Stage 4 (Later): Temporal / SLURM Executor

Deferred. Requirements:
- `submit_hpc_job(run_config)` -> Temporal workflow
- `job_status(job_id)` / `collect_outputs(job_id)`
- Audit trail for all submissions

## Notes / Open Questions

- What dataset roots should be publicly exposed?
- Which scripts are safe for external contributors?
- Are we OK with read-only access to `/storage`?
- Do we need a public vs internal server split?
