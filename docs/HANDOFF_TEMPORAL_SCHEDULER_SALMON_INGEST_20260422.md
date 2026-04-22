# Temporal Scheduler Ingest Contract for `simple_salmon`

This note is the short answer for anyone implementing the scheduler-side JSON
ingest path for the salmon walking-skeleton example.

## Answer

For the salmon example, the scheduler should ingest the **scheduler payload**
that is posted to `POST /start_workflow`.

Use this concrete golden example:

- `/mnt/pikachu/bwb-nextflow-utils/tests/fixtures/simple_salmon/submission_payload.json`

Do **not** treat either of these as the public scheduler ingest contract:

- `/mnt/pikachu/bwb-nextflow-utils/tests/fixtures/simple_salmon/ir.json`
- `/mnt/pikachu/bwb-nextflow-utils/tests/fixtures/simple_salmon/bwb.json`

Those are earlier pipeline stages:

- `ir.json` = canonical Workflow IR
- `bwb.json` = legacy/backend BWB lowering
- `submission_payload.json` = scheduler-ready payload

## Canonical Spec

The current cross-repo boundary registry says the scheduler consumes the
**Scheduler HTTP Payload**, stewarded by `bwb-nextflow-utils`:

- `/mnt/pikachu/agentic-repo-coordinator/registry/boundaries.yaml`
  - boundary id: `temporal-scheduler-http`

The canonical payload spec is here:

- `/mnt/pikachu/bwb-nextflow-utils/plans/modular_architecture_contracts.md`
  - `Contract 5: Scheduler HTTP Payload`
  - `POST /start_workflow body`
- `/mnt/pikachu/bwb-nextflow-utils/tools/workbench_temporal_adapter.py`
  - `bundle_to_workflow_def(...)`
  - `build_submission_payload(...)`

## Exact Files to Read

### Use These

- Golden salmon payload:
  - `/mnt/pikachu/bwb-nextflow-utils/tests/fixtures/simple_salmon/submission_payload.json`
- Boundary registry:
  - `/mnt/pikachu/agentic-repo-coordinator/registry/boundaries.yaml`
- Contract doc:
  - `/mnt/pikachu/bwb-nextflow-utils/plans/modular_architecture_contracts.md`
- Adapter implementation:
  - `/mnt/pikachu/bwb-nextflow-utils/tools/workbench_temporal_adapter.py`

### Context Only

- Workflow IR schema:
  - `/mnt/pikachu/bwb-nextflow-utils/schema/workflow_ir_v0.schema.json`
- Salmon IR example:
  - `/mnt/pikachu/bwb-nextflow-utils/tests/fixtures/simple_salmon/ir.json`
- Salmon BWB example:
  - `/mnt/pikachu/bwb-nextflow-utils/tests/fixtures/simple_salmon/bwb.json`

## Practical Guidance for the Scheduler Teammate

If you are implementing the scheduler or SLURM-side ingest:

1. Parse the shape used in `submission_payload.json`.
2. Expect the top-level request body used by `POST /start_workflow`.
3. Treat `workflow_def` as the scheduler-facing graph contract.
4. Do not require callers to send raw `ir.json`.
5. Do not require callers to send raw `bwb.json`.

## Why This Is the Right Layer

The coordinator RFC is explicit that the scheduler payload is downstream of the
agent execution brief and downstream of Workflow IR:

- `/mnt/pikachu/agentic-repo-coordinator/rfcs/agent-handoff-v0.md`
  - "It is **not** the Workflow IR."
  - "It is **not** the scheduler payload. The scheduler payload is produced by
    the temporal adapter after bundle building and approval."

So the scheduler boundary is:

`bundle -> temporal adapter -> submission payload -> scheduler`

not:

`IR -> scheduler`

and not:

`BWB JSON -> scheduler`
