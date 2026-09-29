# STAR Suite 1.9.5.b Launchpad portability

User request: backport Multiome Launchpad support to 1.9.5.b and remove
workstation paths from Launchpad. Work is isolated on `dev-release-v1.9.5.b`
in `/mnt/pikachu/STAR-suite-v195b-launchpad-20260929`, based on master `4824548`.
Existing tags and the parallel Multiomics integration worktree remain unchanged.

Scope: portable shipped config/schema defaults, bundled Multiome recipe routing,
explicit input/runtime validation, tracked local execution, UI status/cancel,
installation instructions, packaging, and fast automated gates. Reference and
data paths must come from users or their site config. Historical benchmark
scripts and biological results are not rewritten or rerun.

The existing 1.9.5 portable release builds omit Chromap. The UI must distinguish
those builds from a Multiome-capable runtime and fail clearly before processing.
The new Multiomics ownership architecture remains a separate integration.

## Implemented

- Shipped config/workflow forms use explicit user inputs and portable roots.
- Both Multiome recipe IDs use the bundled official snapshot. FASTQ roles, lane
  checks, runtime capability checks and executable hashes are shared by UI/API.
- Local jobs have fresh outputs, status/logs, cancellation and browser reconnect.
- Source/package launcher and relocatable payload are wired into tarballs/.deb.
- Alpine 3.14.3 is bundled with its MIT license so the browser needs no CDN.
- Release artifact builds depend on the reusable Launchpad test workflow.

## Validation

Evidence: `tests/launchpad_portability_output/` (ignored; see `tests/ARTIFACTS.md`).

- Clean `nice -n 10 make -j8 core-portable`: passed. `STAR --version` is
  `1.9.5.b`; `--build-features` reports `chromap_atac:false`.
- Complete MCP test suite: **607 passed**. After browser reconnect refinement,
  final affected Launchpad/config/render/validation suite: **144 passed**.
- Browser test stages the real install payload into a relocated prefix with
  spaces, starts the installed launcher from another directory, blocks external
  browser requests, fills the RNA/ATAC form, launches a synthetic process, then
  reloads and verifies the recipe, inputs and successful run status.
- Real bundled recipe dry-run checks RNA mate order, native ATAC barcode format
  and clipping, including runtime overrides. No aligner is invoked in that test.
- Official snapshots validator: **11 recipes, 10 evidence records passed**.
- JavaScript/Python/shell syntax, workflow YAML and whitespace checks passed.

This is a local release candidate. No tag/release was created and no biological
job was run. A fresh Chromap-enabled engine and full tarball/Debian container
builds were not run in this change; release CI retains those packaging checks.
