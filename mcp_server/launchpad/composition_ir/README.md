# Bulk composition Workflow IR (STAR-suite)

These JSON documents are the **Launchpad composition** source of truth for
`star_genome_generate` and `star_bulk_pe_batch` when building the
`star_index_bulk_pe_then_bulk_de` archetype. They are **not** read from
`bwb-nextflow-utils` test fixtures.

Keep them aligned with the operational semantics of:

- `mcp_server/workflows/star_genome_generate.yaml`
- `mcp_server/workflows/star_bulk_pe_batch.yaml`

If those workflows change in ways that affect constructor handoffs (ports,
types, required inputs), update the matching IR here and run
`pytest mcp_server/tests/test_launchpad.py` (Composition tests).
