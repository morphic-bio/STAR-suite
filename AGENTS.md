# STAR-suite Agent Guide

This document is for coding agents working in this repo. Keep it short and
actionable; link to deeper docs rather than copying them.

## Scope and Goals

- STAR-suite reorganizes STAR into modules while keeping `core/legacy` as the
  single source of truth for the STAR core.
- **Note on naming**: `core/legacy/` is the *current* canonical STAR core tree.
  “legacy” refers to the upstream STAR directory layout we preserve for
  compatibility and maintainability; it does not mean deprecated or “vanilla-only”.
- Changes should preserve STAR CLI compatibility unless explicitly planned.
- Large data and test outputs stay untracked; document their locations in
  `tests/ARTIFACTS.md`.

## Local Agent Override

- If `AGENTS.local.md` exists at the repo root, read it before doing any work
  that depends on local datasets, references, installation paths, or host
  policy.
- `AGENTS.local.md` is intentionally untracked and host-specific. Use
  `AGENTS.local.md.example` as the template.

## Repository Map

- `core/legacy/` - upstream STAR layout (single source of truth).
- `core/features/` - shared feature overlays (vbem, yremove, bamsort, etc.).
- `core/features/process_features/` - vendored `process_features` toolchain.
- `core/features/feature_barcodes/` - standalone tools (`assignBarcodes`, etc.).
- `core/features/libscrna/` - shared EmptyDrops/OrdMag/Occupancy implementations.
- `flex/` - STAR-Flex code and tools.
- `slam/` - SLAM-seq code and tools.
- `scripts/` - suite-level helper scripts (preflight, QC, parity, downstream).
  See `scripts/README.md` for a full catalogue.
- `tests/` - smoke/regression scripts; keep outputs untracked.
- `docs/` - technical summaries and runbooks.

## Build and Smoke Tests

- Core build: `make core` (binary: `core/legacy/source/STAR`).
- Flex tools: `make flex` or `make flex-tools`.
- Feature tools: `make feature-barcodes-tools`.
- CB/UB regression: `tests/run_cbub_regression_test.sh`.
- CRISPR calling test: `tests/test_cr_compat_crispr_calling.sh`.
- Flex smoke: `tests/run_flex_smoke.sh` (if fixtures available).

## Benchmark Hygiene

- Do not run benchmark jobs in parallel with other benchmark jobs on the same
  machine or from the same checkout.
- Serialize benchmarking runs to avoid memory pressure, disk contention, and
  misleading performance results.
- **Benchmark outdirs:** use a **fresh** output directory for each timed run so
  completion artifacts are not confused with an earlier attempt.
- **Completion signals (wrapper-driven runs):** treat **wrapper-written**
  completion artifacts as authoritative when present (for example
  `BENCHMARK_SUMMARY.txt` written only after the wrapped binary exits). If no
  wrapper summary exists, use `Log.final.out` together with a successful wrapper
  exit status or another explicit finished marker—**not** ambiguous partial logs.
- **Do not** infer “still running” from `Log.out`, `Log.progress.out`, or
  `pgrep` alone; those are weak or misleading signals for long Solo jobs.

## Build Hygiene

- **MANDATORY before debugging any crash or regression**: always do a clean
  rebuild first. Stale `.o` files from a previous branch, commit, or partial
  build are a proven source of false segfaults, conflict-guard trips, and
  read-id corruption that disappear after a clean build. Do not spend time
  investigating a failure until you have confirmed it reproduces on a clean
  binary.
- Clean rebuild command:
  `make -C core/legacy/source clean && make -C core/legacy/source -j8 STAR`.
- This applies after switching branches/commits, cherry-picking, reverting
  files, or any operation that changes source without rebuilding all objects.
- This is especially important for Flex/Solo debugging; stale objects can
  produce false segfault/regression signals.

## Data and Artifacts (Do Not Commit)

- External datasets live under `/storage/` (e.g., `/storage/A375/...`).
- Test outputs should be in `tests/*_output*/` or `/tmp/`.
- Add new artifact locations to `tests/ARTIFACTS.md`.

## Key Technical Defaults and Recent Changes

### Unsorted BAM CB/UB injection

- Unsorted BAM CB/UB tags are injected via the same path as sorted BAM.
- Buffered mode uses `g_unsortedTagBuffer` and a `SamtoolsSorter` "noSort"
  mode; output order is deterministic (spill files first, then in-memory).
- Legacy two-pass `--soloAddTagsToUnsorted` path is removed.
- Summary: `docs/UNSORTED_BAM_CBUB_IMPLEMENTATION_SUMMARY.md`.

### CB/UB Independence (Flex)

- CB and UB are handled independently in Flex.
- `status==2` means CB valid, UMI invalid; CB tags still emitted.
- Summary: `docs/CB_UB_INDEPENDENCE_IMPLEMENTATION_SUMMARY.md`.

### CRISPR Feature Calling (CR-compat)

- CR-compat mode runs GMM-based CRISPR calling automatically when Guide Capture
  features are present.
- `--crMinUmi` code default is 3; the CRISPR Guide bundle overrides it to 10. Use 2-3 for lineage barcodes.
- Summary: `docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md`.

### Heatmaps (process_features)

- Cairo PNGs removed; Plotly HTML+JSON outputs only.
- Outputs: `Feature_counts_heatmap.html` + `.json`, `Feature_types_heatmap.html`
  + `.json`.
- Summary: `docs/HEATMAP_REFACTOR_SUMMARY.md`.

### Feature Offset Detection (assignBarcodes / pf_api)

- Auto-detects a global offset from the `pattern` column.
- Errors if multiple offsets are detected (heterogeneity threshold 5%).
- Use `--feature_constant_offset N` for a fixed global offset.
- Use `--force-individual-offsets` for per-feature offsets.
- Docs: `docs/feature_barcodes.md`.

### CR-Compat Benchmark Threading (CRITICAL)

- **Always** use `--dynamicThreadInterface 1` for parallel phases.
- **Always** use `--crAssignConsumerThreads -1` (auto-size). **NEVER** hardcode
  to a small number like 4; this starves the feature assignment phase.
- **Always** use `--crAssignSearchThreads 1` to prevent oversubscription.
- Use `--soloFeatures GeneFull` only (skip Gene/Velocyto) for benchmarks.
- Use `--outSAMtype None` unless BAM output is specifically needed.
- Full reference command: `docs/feature_barcodes.md` § "Optimized Benchmark
  Parameters".

## Flex Integration Notes

- Flex now uses `libscrna` for EmptyDrops/OrdMag/Occupancy (no duplicate
  implementations).
- Ensure Flex builds link `libscrna` and include `core/features/libscrna/include`.
- Plan: `plans/refactor_flex_plan.md`.

## CR-compat GEX Parity Notes

- Use `GeneFull` for GEX parity (CR includes introns since v7).
- `--soloCrGexFeature` controls which GEX MEX is merged in CR-compat mode.
- CB/UB tags are independent (CB can be present without UB and vice versa).

### Poly-G Trimming

- `--clip3pPolyG yes|no|auto` (default `auto`) trims NovaSeq/NextSeq poly-G
  artifacts in CellRanger4 mode. Without this, poly-G reads inflate LINC00486
  and destroy gene-level Pearson correlation.
- Summary: `docs/HANDOFF_UCSF_FULL_GENE_PEARSON_ANOMALY_20260220.md`.

### OrdMag Cell Calling (Bootstrap)

- Non-Flex libscrna path uses CR9-style bootstrapped `recovered_cells`
  estimation (100 samples) instead of hardcoded `nExpectedCells=3000`.
- Flex path keeps fixed `nExpectedCells=3000` (CR 7.1 defaults).
- EmptyDrops MC simulations: 100K for non-Flex, 10K for Flex.
- BH FDR correction enabled for non-Flex; raw p-value for Flex.
- Validated: UCSF iPSC2 full sample Jaccard 0.99, gene Pearson 0.997.

## MCP Server (Agent Tooling)

An MCP server is available for automated agent workflows:

```bash
pip install -r mcp_server/requirements.txt
export MCP_AUTH_TOKEN="your-token"
python3 -m mcp_server.app --transport http --host 0.0.0.0 --port 8765
```

Convenience wrapper for humans (start/stop Launchpad + MCP on one port):

```bash
bash scripts/launchpad_server.sh up
bash scripts/launchpad_server.sh down
```

### Discovery tools (no auth required when `public_discovery: true`)

| Tool | Purpose |
|------|---------|
| `list_datasets` | Enumerate configured datasets with paths and metadata |
| `list_test_suites` | List test suites, scripts, and runnability status |
| `find_docs` | Search documentation by keyword |
| `find_tests` | Search test scripts by keyword |
| `list_workflows` | List all registered workflows with summaries |
| `describe_workflow` | Full workflow details: stages, parameter groups, caveats |
| `get_workflow_scripts` | Scripts composing a workflow: entry + helpers (provenance detail requires auth) |
| `get_workflow_parameter_schema` | Machine-readable parameter definitions with types, defaults, constraints |
| `scaffold_workflow_schema` | Parse a shell script and generate a draft workflow schema YAML |
| `validate_draft_workflow_schema` | Validate a draft schema against the WorkflowSchema model |

**Workflow visibility**: Each workflow in `config.yaml` can set `visibility: public`
(default) or `visibility: private`. Private workflows are omitted from
`list_workflows` and return "unknown workflow" from `describe_workflow` /
`get_workflow_scripts` / `get_workflow_parameter_schema` unless the caller
provides a valid `auth_token`. This is content-level filtering — the same tool
returns different results based on auth status.

**Discovery detail levels**: `get_workflow_scripts` applies response redaction
based on auth status:
- **Public / unauthenticated**: returns safe structural information only —
  `workflow_id`, `title`, `entry_script`, per-script `role`, relative `path`,
  `description`, `language`, `exists`, and `workflow_schema` in provenance.
  Host-specific fields (`absolute_path`, `repo_root`, `git_commit`,
  `git_remote`) are omitted.
- **Authenticated / trusted-local**: returns full detail including
  `absolute_path` on each script and full provenance (`repo_root`,
  `git_commit`, `git_remote`). This level is intended for same-host agents
  and trusted local tooling.

### Authenticated tools (require `auth_token`)

| Tool | Purpose |
|------|---------|
| `validate_workflow_parameters` | Validate parameter values (types, paths, constraints) |
| `render_workflow_command` | Render a validated parameter set into a shell command (`argv` + env) |
| `preflight` | Pre-run checks (disk space, binaries, script permissions) |
| `run_script` | Execute an allowlisted script with timeout and logging |
| `collect_outputs` | Retrieve run status, logs, and output file inventory |
| `reload_config` | Hot-reload config and workflow schemas |

### Workflow parameter pipeline

The recommended agent workflow for running scripts:

1. `list_workflows` / `describe_workflow` — discover available workflows
2. `get_workflow_scripts` — get entry + helper scripts with provenance (for encoders)
3. `get_workflow_parameter_schema` — get parameter types, defaults, constraints
4. `validate_workflow_parameters` — validate user-supplied parameters
5. `render_workflow_command` — get the exact `argv` and env overrides
6. `preflight` — verify disk space, binaries, permissions
7. `run_script` — execute; get `run_id` for tracking
8. `collect_outputs` — retrieve results when done

### Schema authoring (propose-only)

To add new workflows, agents can use `scaffold_workflow_schema` to parse an
existing shell script and generate a draft YAML, then `validate_draft_workflow_schema`
to check it. The draft must be committed to the repo and loaded via `reload_config`
to become active — these tools never modify the running config.

### Key features

- All paths validated against trusted roots
- Job queue (1 concurrent, 10 queued)
- Timeout handling with process group cleanup
- Logs stored in `plans/artifacts/mcp_runs_YYYYMMDD/`

Documentation: `mcp_server/README.md`  
Configuration: `mcp_server/config.yaml`  
Workflow schemas: `mcp_server/workflows/`

### MCP Usage for New Agents

- Prefer MCP tools for discovery/preflight/run rather than ad-hoc shell scripts.
- Endpoints differ by client:
  - **Codex** (streamable-HTTP): `POST /`
  - **Cursor** (SSE): `GET /sse` + `POST /messages`
- If running on a remote host, local SSH tunnels are not visible to agents.
  Use a public tunnel (e.g., `cloudflared tunnel --url http://localhost:8765`)
  so agents can reach the MCP server.

## Docs to Check First

- `docs/UNSORTED_BAM_CBUB_IMPLEMENTATION_SUMMARY.md`
- `docs/CRISPR_FEATURE_CALLING_IMPLEMENTATION_SUMMARY.md`
- `docs/HEATMAP_REFACTOR_SUMMARY.md`
- `docs/Github-actions.md`
- `docs/Star-binary-distribution.md`
- `docs/feature_barcodes.md`
- `docs/todos`
- `tests/ARTIFACTS.md`
- `mcp_server/README.md`

## Commits

- Never add `Co-Authored-By` trailers to commit messages. All commits are
  authored by the human developer.
- The default branch is `master`. Push to `master`, not `main`.

## Branching and Merges

- Feature branches merge into `perturb`, then merge into `master`.
- Do not squash-merge branches that touch shared core files; preserve the DAG
  with `git merge --no-ff` so later integrations keep a usable merge base.
- Keep large binaries and datasets untracked; update `.gitignore` if needed.

### Safe Merge Policy

Large shared-file merges have caused silent regressions by dropping critical
code when two branches modify the same file. Follow these rules:

- **Never squash-merge branches that touch `assignBarcodes.c` or other large
  shared files.** Use `git merge --no-ff` instead; this preserves the DAG so
  git can produce proper three-way conflict markers.
- **Pre-merge diff audit**: before any squash merge, verify no expected
  functions were dropped:
  ```bash
  # After staging the squash but BEFORE committing:
  git diff HEAD -- core/features/process_features/src/assignBarcodes.c \
    | grep '^-.*pf_search_hash_offsets\|^-.*pf_single_offset_hash_search' \
    && echo "FATAL: tiered search removed" && exit 1
  ```
- **Post-merge symbol check**: after every merge that touches
  `process_features`, confirm critical fast-path symbols are present:
  ```bash
  grep -c 'pf_search_hash_offsets\|pf_single_offset_hash_search' \
    core/features/process_features/src/assignBarcodes.c
  # Must be >= 4 (2 definitions + 2 call sites)
  ```
- **Prefer rebase-merge** (`git rebase <target> && git merge --no-ff`) for
  branches touching core files. This gives clean linear history while
  preserving merge-base information.

## GitHub Actions CI/CD Policy

- Do not publish images or release artifacts on every push.
- Pull requests: run fast checks only (build sanity + Tier A smoke); no publish.
- Push to `dev`: run integration checks and optionally publish `dev-<sha>` images.
- Push to `master`: run required checks and publish multi-arch images (`amd64`,
  `arm64`) to stable tags (`latest`, `master-<sha>`).
- Tags `v*`: run release pipeline (multi-arch publish + GitHub Release artifacts
  + source package upload for PPA build).
- CI path filters are enabled for `ci-pr.yml`, `ci-dev.yml`, and
  `ci-master.yml`; these workflows run only when build/test/release infra paths
  change (see exact globs in `docs/Github-actions.md`).
- `release.yml` remains tag-triggered on `v*` and is intentionally not path-scoped.
- Fixture-heavy tests (Tier B) should be run only when fixtures are available
  (self-hosted/scheduled/manual), not on all PRs.
- Details: `docs/Github-actions.md`.

## Output Hygiene

- Do not commit generated binaries or test outputs.
- If a new test creates outputs, add its location to `tests/ARTIFACTS.md`.
