# UCSF STAR-suite Parameter Schema Audit

> Date: 2026-04-08
> Scope: completeness audit of `mcp_server/workflows/ucsf_star_suite_production.yaml`

## Method

Traced every parameter consumed by each script in the UCSF production
workflow chain:

1. `scripts/paper/run_ucsf_corrected_production_workflow.sh` (entry, 21 CLI flags)
2. `scripts/run_ucsf_perturb_yremove_batch.sh` (batch runner, env vars + forwarded flags)
3. `scripts/run_scrna_downstream_gene_full_velocyto.sh` (downstream, env vars)
4. `scripts/run_remote_cellbender_rsync.sh` (remote CellBender)
5. `scripts/run_remote_cellbender_scan.sh` (watcher)
6. Six Python helper scripts

Classified each parameter using the framework from the input surface
classification runbook (canonical workflow input / source customization /
runtime-infrastructure / internal).

## Previously Exposed (21 params)

All entry-script CLI flags were already in the schema: `samples`,
`all_samples`, `dataset_root`, `feature_ref`, `genome_dir`,
`solo_cb_whitelist`, `cr_whitelist`, `star_bin`, `threads`,
`cellbender_cpu_cores`, `downsample_reads`, `downsample_seed`,
`cellbender_gpu`, `trim_qc`, `star_only`, `dry_run`,
`globus_src_endpoint`, `globus_dst_endpoint`, `globus_dst_root`,
`globus_poll_seconds`, `out_root`.

## Newly Exposed (6 params, total now 27)

| Parameter | Type | Source script | Delivery | Default | Why exposed |
|-----------|------|---------------|----------|---------|-------------|
| `sample_map` | file | batch_runner | env `SAMPLE_MAP` | null | Overrides FASTQ lookup; affects which reads are processed |
| `trim_qc_max_reads` | int | batch_runner | env `STAR_TRIM_QC_MAX_READS` | 250000 | Controls trim-QC sampling depth |
| `min_genes` | int | downstream | env `MIN_GENES` | 200 | QC filter: minimum genes per cell |
| `max_genes` | int | downstream | env `MAX_GENES` | 2500 | QC filter: base max genes (adaptive overrides) |
| `mt_pct_cutoff` | float | downstream | env `MT_PCT_CUTOFF` | 5 | QC filter: mitochondrial percentage cutoff |
| `n_mad` | float | downstream | env `N_MAD` | 3 | Adaptive filtering: MADs for max_genes threshold |

All 6 use `env_var` mapping with `skip_when_default: true`, so they only
appear in `render_workflow_command` output's `env_overrides` when the user
sets a non-default value.

## Not Exposed (intentionally)

### Locked-on: `adaptive_filter`

The entry script hardcodes `--adaptive-filter` to the batch runner when
`star_only=false`. There is no CLI path to disable it through the top-level
wrapper. Documented in schema caveats but not a settable parameter.

### Internal / forwarded-only

| Parameter | Script | Reason not exposed |
|-----------|--------|--------------------|
| `cr_config` | batch_runner | Auto-derived from sample directory layout |
| `cr_sample_id` | batch_runner | Auto-derived; same as sample_id |
| `downstream_output_name` | batch_runner | Derived from sample_id |
| `CELLBENDER_LAYER` | downstream | Fixed layer name (internal convention) |
| `CELLBENDER_FLAGS` | downstream | Expert override; not user-facing |
| Container image vars | downstream | Runtime/infrastructure, not workflow semantics |

### Runtime / infrastructure

| Parameter | Script | Reason |
|-----------|--------|--------|
| Remote host / SSH config | remote_cellbender | Infrastructure routing |
| rsync paths | remote_cellbender | Infrastructure routing |
| Scan polling interval | remote_scan | Infrastructure timing |

## Bug Fix: `skip_when_default` for env_var params

While writing tests, discovered that `render_workflow_command` was not
applying `skip_when_default` to `env_var` params — they always appeared in
`env_overrides` even at default values. Fixed in `tools/workflows.py` to
apply the same default-check logic for env_var rendering.

## Schema Changes Summary

- **Parameters**: 21 → 27
- **Parameter groups**: 6 → 7 (added `downstream_qc`)
- **Constraints**: 5 → 7 (added `trim_qc_max_reads` dependency, extended positive checks)
- **Caveats**: 3 → 5 (adaptive filtering always-on, env_var consumption)
- **Flag order**: updated to include all 6 new params

## Test Coverage

- 16 new tests in `test_ucsf_workflow_e2e.py` (`TestUCSFNewlyExposedParams`)
- Updated scaffold consistency test to exclude sub-script params from
  scaffolder coverage check (scaffolder only parses entry-script getopts)
- Full suite: 430 tests, all green

## Remaining Gaps

- Multi-sample bridge support (bridge currently requires exactly one sample)
- `CELLBENDER_FLAGS` could be exposed as an expert/advanced parameter in a
  future pass if users need to tune CellBender options
- No schema for the Python helper scripts' own arguments (they are called
  internally by the bash orchestration, not user-facing)
