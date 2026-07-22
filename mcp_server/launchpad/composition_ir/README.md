# Bulk composition Workflow IR (STAR-suite)

These JSON documents are the **Launchpad composition** source of truth for
`star_genome_generate` and `star_bulk_pe_batch` when building the
`star_index_bulk_pe_then_bulk_de` archetype. They are **not** read from
`bwb-nextflow-utils` test fixtures.

Keep them aligned with the operational semantics of:

- `mcp_server/workflows/star_genome_generate.yaml`
- `mcp_server/workflows/star_bulk_pe_batch.yaml`

The SLAM PE smoke/production and SLAM DESeq2 entries are shell workflow
recipes exposed by Launchpad and MCP. They do not currently have composition IR
documents because they orchestrate host-local datasets, Globus transfer, and
container build actions rather than a portable bulk STAR workflow graph.

If those workflows change in ways that affect constructor handoffs (ports,
types, required inputs), update the matching IR here and run
`pytest mcp_server/tests/test_launchpad.py` (Composition tests).
